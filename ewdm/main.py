#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

import numpy as np
import xarray as xr
import logging

from tqdm import tqdm

from .wavelets import cwt, xwt
from .density import estimate_directional_distribution, estimate_radial_distribution
from .helpers import get_sampling_frequency
from .sources import SpotterBuoysDataSource, CDIPDataSourceRealTime
from .parameters import VARIABLE_NAMES


# logging.basicConfig(level=logging.INFO)
logging.basicConfig(level=logging.WARN)
logger = logging.getLogger("main")


RADTODEG = 180. / np.pi
DEGTORAD = np.pi / 180.
GRAV = 9.8


class _BaseClass(object):
    """Base class to different estimation methods"""
    def __init__(
            self,
            dataset: xr.Dataset,
            fs: float = None,
            interpolate: bool = True,
            max_nan_ratio: float = 0.1,
            max_time_gap: str = "10s",
            normalise: bool = True
        ) -> xr.Dataset:
        """Initialise class"""
        
        self.dataset = dataset
        if fs is None:
            try:
                self.fs = self.dataset.sampling_rate
                logger.info(f"Sampling frequency from dataset: {self.fs} Hz")
            except AttributeError:
                self.fs = get_sampling_frequency(self.dataset["time"])
                logger.info(f"Sampling frequency estimated: {self.fs} Hz")
        else:
            self.fs = fs

        self.interpolate = interpolate
        self.max_nan_ratio = max_nan_ratio
        self.max_time_gap = max_time_gap
        self.normalise = normalise


    def interpolate_dataset(self, dataset, max_nan_ratio, max_time_gap):
        """Interpolate dataset if it contanins invalid values
        
        Arguments:
            dataset (xr.Dataset): It should contain surface_elevation
            max_nan_ratio (float): Maximum threshold for invalid to
                valid data ratio. If dataset invalid values supasses
                this threshold, the function return a dataset full of nan.
            max_time_gap (str): Maximum tolerable time gap to interpolated.

        Retunrs:
            xr.Dataset interpolated
        """

        ntime = len(dataset["time"])
        nan_ratio = (dataset["surface_elevation"].isnull().sum() / ntime)
        if nan_ratio > max_nan_ratio:
            logger.warning(
                f"The invalid data ratio is {nan_ratio} and is "
                f"greater than the allowed value {max_nan_ratio}. "
                f"Returning dataset full of nans."
            )
            return dataset * np.nan
        else:
            logger.info("The dataset will be linearly interpolated.")
            return dataset.interpolate_na(
                dim="time", method="quadratic", max_gap=max_time_gap,
                fill_value="extrapolate"
            )



class Arrays(_BaseClass):
    """Perform EWDM for spatial arrays of surface eleavtion data.

    Arguments:
        dataset (xr.Dataset): Dataset containing input data. Typical wave staff
        measurements are characterised by sea surface elevation data at
        different spatial location. ADCP along-beam echo-based surface elevation
        can also be considered as a form of spatial arrays. The convenction
        followed for variable names is:

        .. code-block:: python

           <xarray.Dataset>
           Dimensions:            (time: 4096, element: 6)
           Coordinates:
             * time               (time) float64
             * element            (element) int64 0 1 2 3 4 5
           Data variables:
               surface_elevation  (time, element) float64
               position_x         (element) float64
               position_y         (element) float64
           Attributes:
               sampling_rate:  4

    Returns:
        xr.Dataset: Dataset containing produced directional spectra
        and directional spreading function.
    """

    def __init__(
            self,
            dataset: xr.Dataset,
            fs: float = None,
            interpolate: bool = True,
            max_nan_ratio: float = 0.1,
            max_time_gap: str = "10s",
            normalise: bool = True
        ):

        # initialise child class with parent args
        super().__init__(
            dataset, fs, interpolate, max_nan_ratio, max_time_gap, normalise
        )

        # number of elements, unique possible equations and pairs
        self.npoints = len(self.dataset["element"])
        self.neqs = self.npoints * (self.npoints-1) // 2
        self.npairs = (self.neqs * (self.neqs - 1)) // 2


    @classmethod
    def from_numpy(
        cls,
        time: np.ndarray,
        surface_elevation: np.ndarray,
        position_x: np.ndarray,
        position_y: np.ndarray,
        **kwargs
    ):
        """
        Create an instance of Arrays from numpy arrays.

        Arguments:
            time: Time array
            surface_elevation: Surface elevation array. The size of this object
                should be consistent with the number of elements in the array.
            position_x: Easting coordinate of array elements in metres
            position_y: Northing coordinate of array elements in metres
        """

        # determine number of elements
        if len(position_x) == len(position_y):
            # element numbering is always assumed from 0 to number of elements
            elements = np.arange(len(position_x))
        else:
            msg = "Number of elements in `x` and `y` are not consistent"
            logger.exception(msg)
            raise Exception(msg)

        # TODO: check surface_elevation size
        
        # create dataset from numpy arrays
        dataset = xr.Dataset(
            data_vars = {
                "surface_elevation": (["time", "element"], surface_elevation),
                "position_x": ("element", position_x),
                "position_y": ("element", position_y)
            },
            coords = {"time": time, "element": elements},
        )

        # store sampling_rate as a global attribute
        if 'fs' in kwargs:
            dataset.attrs = {'sampling_rate': kwargs["fs"]}

        return cls(dataset, **kwargs)

    
    
    def wavelet_coefficients(self, dataset: xr.Dataset) -> xr.Dataset:
        """Estimate wavelet coefficients

        This method takes the dataset containing the sea surface
        elevation data at different elements of the array and returns
        the corresponding wavelet complex coefficients for each
        element.
        """

        if "surface_elevation" in dataset:
            coeffs = xr.Dataset()
            for element in dataset["element"]:
                cwt_result = cwt(
                    dataset["surface_elevation"].sel(element=element),
                    freqs=self.freqs, fs=self.fs
                )
                coeffs[element.item()] = xr.DataArray(
                    cwt_result,
                    dims=["frequency", "time"],
                    coords={"frequency": self.freqs, "time": dataset["time"]},
                    attrs={"sampling_rate": self.fs}
                )
            return coeffs.to_array(dim='element')
            

        else:
            raise Exception(
                "Local wavelet power cannot be computed because "
                "`surface_elevation` data is not available in "
                "the dataset."
            )

    
    def array_geometry(self, dataset: xr.Dataset) -> np.ndarray:
        """Compute array of spatial differences for each point.
            
        This method takes the dataset containing the position (x
        and y) of each element of the spatial array and calculate
        the position difference vector.
        """

        try:
            x = dataset["position_x"].data
            y = dataset["position_y"].data
        except KeyError:
            raise Exception(
                "Local wavelet power cannot be computed because "
                "`position_x` and `position_y` data is not available "
                "in the dataset."
            )

        logger.info(f"Computing spatial distances in array.")

        dx = np.zeros((self.neqs, 2))
        index = 0
        for m in range(self.npoints):
            for n in range(m + 1, self.npoints):
                dx[index, 0] = x[m] - x[n]
                dx[index, 1] = y[m] - y[n]
                index += 1

        return dx
    
    
    def compute_angle(self, complex_data):
        """Compute and wrap angle in radians from complex argument."""

        return np.arctan2(complex_data.imag, complex_data.real)
        # return (angle - np.pi) % (2 * np.pi) - np.pi


    def phase_differences(
        self, coeffs: xr.Dataset, cross_wavelet=False
    ) -> np.ndarray:
        """Compute array with phase differences between array pairs"""

        # get freq nd time lenghts
        ntimes = len(coeffs["time"])
        nfreqs = len(coeffs["frequency"])

        # initialise delta phi array
        dphi = np.zeros([self.neqs, nfreqs, ntimes])

        # loop for each pair of points
        index = 0
        for m in range(self.npoints):
            for n in range(m + 1, self.npoints):
                logger.info(f"Processing for pair {index}: {m},{n}")
                if cross_wavelet:
                    dphi[index, :, :] = (
                        self.compute_angle(
                            coeffs.isel(element=m) * 
                            coeffs.isel(element=n).conj()
                        ).data
                    )
                else:
                    dphi[index, :, :] = (
                        self.compute_angle(coeffs.isel(element=m)) - 
                        self.compute_angle(coeffs.isel(element=n))
                    ).data
                index += 1

        return dphi

    def compute_wavenumbers(self, dx, dphi, solver="lstsq"):
        """Compute wavenumber vector from phase differences and distances.
        
        This method take the point distances and the phase differences 
        arrays and returns the wavenumber vector components and residuals.
        """

        # assuming that dphi array is neqs, nfreqs, ntimes
        logger.info(f"phase array size is {dphi.shape}")
        _, nfreqs, ntimes = dphi.shape

        # initialise delta kx, ky and epsilon
        kx = np.zeros([nfreqs, ntimes])
        ky = np.zeros([nfreqs, ntimes])
        epsilon = np.zeros([nfreqs, ntimes])

        if solver in ["lstsq", "least-squares", 1]:

            logger.info("Finding least-square solution")
            
            # loop for each frequency
            for ifrq in range(nfreqs):
                kk, residuals, _, _ = np.linalg.lstsq(
                    dx, dphi[:,ifrq,:], rcond=self.rcond
                )
                kx[ifrq,:] = kk[0,:]
                ky[ifrq,:] = kk[1,:]
                epsilon[ifrq,:] = residuals
                
            return kx, ky, epsilon

        # solve each pair of elemens separately
        elif solver in ["pair-wise", "pairwise", 2]:

            logger.info("Finding pair-wise solution")
            
            # loop for each frequency
            for ifrq in range(nfreqs):
                kk_set = []
                for i in range(self.neqs-1):
                    for j in range(i+1, self.neqs):
                        dot_ij = (
                            (np.dot(dx[i], dx[j])) /
                            (np.linalg.norm(dx[i]) * np.linalg.norm(dx[j]))
                        )
                        if np.abs(dot_ij) < 0.5:
                            dx_ij = np.array([dx[i], dx[j]])
                            dphi_ij = np.array(
                                [dphi[i,ifrq,:], dphi[j,ifrq,:]]
                            )
                            kk = np.linalg.solve(dx_ij, dphi_ij)
                            kk_set.append(kk)

                kk_set = np.array(kk_set)
                kx[ifrq,:] = np.mean(kk_set, axis=0)[0]
                ky[ifrq,:] = np.mean(kk_set, axis=0)[1]
                # epsilon[ifrq,:] = np.std(kk_set, axis=0)

            return kx, ky, epsilon


        else:
            raise Exception(
                "Solver argument only accepts the following options:\n"
                "1: `least-squares` `lstsq`\n"
                "2: `pair-wise` `pairwise`\n"
            )

    def estimate_directional_distribution(self, dataset) -> xr.Dataset:
        """Estimate the directional distribution of wave energy.

        This method calculates the wavelet power and local wave direction using
        the local estimates of wavenumber vector and then estimates the
        directional distribution function and directional spectra.

        Returns:
            xr.Dataset: Dataset containing the directional spectrum, directional
            distribution, and frequency spectrum.

        Raises:
            Exception: If `use` is not one of `displacements`, `velocities`,
            `accelerations` or `slopes`.
        """

        # if self.interpolate:
            # _dataset = self.interpolate_dataset(
                # dataset, self.max_nan_ratio, self.max_time_gap
            # )
        # else:
            # _dataset = dataset.copy()
        
        _dataset = dataset.copy()

        # first obtaing wavelet coefficients
        coeffs = self.wavelet_coefficients(_dataset)

        # array geometry
        dx = self.array_geometry(_dataset)

        # Largest baseline (maximum pairwise separation) of the array, in
        # metres. The coarsest wavenumber the array can resolve is of order
        # 2*pi / L_baseline. This is exposed as `self._L_baseline` so a
        # caller can, if they wish, choose an explicit `bandwidth_floor`
        # informed by the array geometry; it is not applied automatically,
        # because an aperture wavelength is far wider than a sensible KDE
        # kernel and would over-smooth the spectrum.
        if dx.size:
            self._L_baseline = float(np.max(np.hypot(dx[:, 0], dx[:, 1])))
        else:
            self._L_baseline = np.nan

        # TODO: add the following paramaters as class attrs
        # xmin = np.min(np.abs(dx[:, 0] + 1j*dx[:, 1]))
        # xmax = np.max(np.abs(dx[:, 0] + 1j*dx[:, 1]))

        # kmin = GRAV * (1 / (self.fs * xmax))**2  # wave_speed < sampling_speed
        # kmax = np.pi / xmin  # when k*dx = pi

        # fmin = np.sqrt(GRAV * kmin) / (2*np.pi)
        # fmax = np.sqrt(GRAV * kmax) / (2*np.pi)

        # period_bounds = [1/fmax, 1/fmin]
        # wavelength_bounds = [2*np.pi/kmax, 2*np.pi/kmin]

        # phase differences
        dphi = self.phase_differences(coeffs, cross_wavelet=self.cross_wavelet)
        dphi = (dphi - np.pi)  % (2* np.pi) - np.pi

        # limit = np.pi
        # min_phase = 0
        # dphi[dphi == 0] = min_phase
        # dphi[dphi >  limit] = min_phase
        # dphi[dphi < -limit] = min_phase

        # compute wave number components and resiudal
        # TODO: restrict to residuals less than certain threshold
        kx, ky, residuals = self.compute_wavenumbers(
            dx, dphi, solver=self.solver
        )
        
        # compute local wave direction
        theta =  xr.DataArray(
            RADTODEG * np.arctan2(ky, kx),
            dims = ["frequency", "time"],
            coords={"frequency": self.freqs, "time": _dataset["time"]},
        )

        # local wavenumber magnitude k = |(kx, ky)| in rad/m. This is
        # solved from the spatial phase gradients and is therefore an
        # observed quantity, independent of any dispersion assumption.
        wavenumber = xr.DataArray(
            np.hypot(kx, ky),
            dims = ["frequency", "time"],
            coords={"frequency": self.freqs, "time": _dataset["time"]},
        )

        # local inverse phase speed nu = k / omega = k / (2 pi f),
        # following Björkqvist et al. (2019), in s/m.
        omega = 2 * np.pi * self.freqs
        nu = wavenumber / xr.DataArray(
            omega, dims=["frequency"], coords={"frequency": self.freqs}
        )

        # compute wavelet power from wavelet coefficients
        data_std = _dataset["surface_elevation"].std("time")
        power = np.abs(coeffs) ** 2
        if self.normalise:
            # rescale so that the integral of the spectrum equals the
            # surface-elevation variance: int <|W|^2>_t df -> sigma^2.
            # The factor is variance / wavelet-energy (sigma^2 / E), not
            # sigma / sqrt(E), which would only be exact when E == sigma^2.
            wavelet_energy = power.mean("time").integrate("frequency")
            mean_power = (power * data_std**2 / wavelet_energy).mean("element")
        else:
            mean_power = power.mean("element")

        # frequency-direction spectrum (default, unchanged behaviour)
        if self.coordinate in ("frequency", None):
            return estimate_directional_distribution(
                mean_power, theta, dd=self.dd, kappa=self.kappa
            )

        # wavenumber-direction or inverse-phase-speed-direction spectrum.
        # The radial-KDE bandwidth is left to Silverman's data-driven rule
        # (see `density._gaussian_radial_kde`), which together with the
        # log-spaced default grid yields F(k)/Q(nu) spectra about as smooth
        # as S(f) while keeping the spectral peak unbiased. An explicit
        # `self.bandwidth_floor` (default None) is forwarded for callers who
        # want to impose a physical lower bound on the kernel width; when
        # None the kernel is only floored at the local bin spacing to guard
        # against degenerate (zero-spread) radial samples.
        elif self.coordinate == "wavenumber":
            return estimate_radial_distribution(
                mean_power, theta, wavenumber, "wavenumber",
                bins_radial=self.bins_radial,
                dd=self.dd, kappa=self.kappa, bandwidth=self.bandwidth,
                bandwidth_floor=self.bandwidth_floor
            )

        elif self.coordinate == "nu":
            return estimate_radial_distribution(
                mean_power, theta, nu, "nu",
                bins_radial=self.bins_radial,
                dd=self.dd, kappa=self.kappa, bandwidth=self.bandwidth,
                bandwidth_floor=self.bandwidth_floor
            )

        else:
            raise Exception(
                "`coordinate` should be either `frequency`, "
                "`wavenumber` or `nu`."
            )


    def compute(
            self,
            omin: float = -5,
            omax: float = None,
            nvoice: float = 16,
            cross_wavelet: bool = False,
            solver: str = "lstsq",
            rcond: float = None,
            dd: float = 5.0,
            kappa: float = 36.0,
            use: str = "displacements",
            block_size: str = "30min",
            coordinate: str = "frequency",
            bins_radial: np.ndarray = None,
            bins_per_octave: int = 24,
            bandwidth: str = "silverman",
            bandwidth_floor: float = None,
        ) -> xr.Dataset:
        """Perform computation using specified parameters.

        Args:
            omin (float, optional): Minimum octave (default is -5).
            omax (float, optional): Maximum octave. If None, it is
                automatically determined based sampling frequency. The final
                frequency array is logaritmically distributed from
                `2**omin` to `2**omax`.
            nvoice (float, optional): Number of voices for the computation
                (default is 16).
            cross_wavelet (bool, optional): Whether or not use cross-wavelet
                product to compute phase differences between array pairs. If
                False then arithmetic substraction is used (default).
            solver (str, option): Wavenumber solver method. Two options are
                available: 1 or 'lstsq' or 'least-square': Least-square
                solution of all possible pair. combination. 2 or 'pairwise'
                or 'pair-wise': Pair-wise solution and averaging...
            rcond (float, optional): If solver is `lstsq` this argument is passed
                to `np.linalg.lstsq` function. It is a cut-off ration for small
                singular values of matrix `dphi`.
            dd (float, optional): Directional resolution in degrees
                (default is 5 degrees).
            kappa (float, optional): Smoothness parameter for Kernel
                Density Estimation. Small values of `kappa` produce oversmooth
                results (default is 36.0).
            use (str, optional): Type of data to perform estimation.
                It should be should be either `displacements`, `velocities`
                or `accelerations`." (default is "displacements").
            block_size (str): If dataset contains more than one hour of data,
                split dataset into blocks of `block_size` and perform
                computation over each block. The resulting output will have a
                time dimension. It is advisable to choose values of no more than
                half-hour. Default `block_size="30min"`.
            coordinate (str, optional): Independent radial coordinate of the
                output spectrum. One of `frequency` (default, the usual
                frequency-direction spectrum), `wavenumber`, or `nu` (the
                inverse phase speed nu = k / omega, yielding the Q(nu, theta)
                spectrum of Björkqvist et al., 2019). The latter two are
                available only for `Arrays` because they rely on the
                wavenumber vector solved from the spatial phase gradients.
            bins_radial (np.ndarray, optional): Bin centres for the radial
                axis when `coordinate` is `wavenumber` (rad/m) or `nu`
                (s/m). If None, a sensible range is derived from the data.
                Ignored when `coordinate="frequency"`.
            bins_per_octave (int, optional): Number of bins per octave for
                the default log-spaced radial grid used when `bins_radial`
                is None (default 24). Logarithmic spacing matches the
                wavelet frequency grid and yields smoother F(k)/Q(nu)
                spectra. Ignored when `bins_radial` is supplied or when
                `coordinate="frequency"`.
            bandwidth (str or float, optional): Bandwidth for the Gaussian
                kernel used along the radial axis. Either a number or the
                string `silverman` (default). Ignored when
                `coordinate="frequency"`.
            bandwidth_floor (float, optional): Physical lower bound on the
                radial-kernel bandwidth, in the units of the radial
                coordinate (rad/m for `wavenumber`, s/m for `nu`). Default
                None, in which case Silverman's data-driven bandwidth is
                used and only floored at the local bin spacing to guard
                degenerate samples. Provide a value to impose a known
                resolution limit (note the array baseline length is
                available as `self._L_baseline` after `compute`). Ignored
                when `coordinate="frequency"`.

        Returns:
            xr.Dataset: Dataset containing the directional spectrum, directional
                distribution, and frequency spectrum.
        """

        # the maximum frequency is given by the nyquist frequency
        if omax is None:
            omax = int(np.log2(self.fs / 2))

        self.omin = omin
        self.omax = omax
        self.nvoice = nvoice

        # fourier equivalent frequencies
        self.freqs = 2. ** np.linspace(omin, omax, nvoice*abs(omin-omax)+1)

        # determine phasediff and wavenumber solver 
        self.cross_wavelet = cross_wavelet
        self.solver = solver
        self.rcond = rcond

        # directional resolution
        self.dd = dd
        self.kappa = kappa

        # data used for estimation
        self.use = use
        self.block_size = block_size

        # radial coordinate of the output spectrum
        self.coordinate = coordinate
        self.bandwidth = bandwidth
        self.bandwidth_floor = bandwidth_floor

        # number of bins per octave for the default (log-spaced) radial
        # grids. Logarithmic spacing mirrors the wavelet frequency grid
        # `self.freqs` (also log-uniform), placing bins where the spectral
        # structure is and rendering smoothly on the usual log-log axes.
        # `nperoct` only controls how finely the radial KDE is *sampled*;
        # the physical reliability of the result is governed by the kernel
        # bandwidth and its `bandwidth_floor` (see `compute_spectrum`), not
        # by the bin count.
        self.bins_per_octave = bins_per_octave

        # default radial bins if none are provided. We keep the
        # Björkqvist et al. (2019) span — wavenumber over roughly
        # [1/12, 9] rad/m and inverse phase speed nu = k / omega bounded by
        # the slowest resolved wave — but distribute the bins
        # logarithmically rather than linearly so the omnidirectional F(k)
        # and Q(nu) spectra look as smooth as the frequency spectrum S(f).
        # The lower bound stays one Björkqvist bin-width above zero so the
        # polar directional spectrum (which carries a 1/r Jacobian) has no
        # singular bin at the origin.
        def _logbins(lo, hi, nperoct):
            noct = np.log2(hi / lo)
            npts = int(np.ceil(nperoct * noct)) + 1
            return 2.0 ** np.linspace(np.log2(lo), np.log2(hi), npts)

        if bins_radial is None and coordinate == "wavenumber":
            k_lo = 1.0 / 12.0
            k_hi = 9.0
            self.bins_radial = _logbins(k_lo, k_hi, self.bins_per_octave)
        elif bins_radial is None and coordinate == "nu":
            nu_lo = 1.0 / 50.0
            nu_hi = 9.0 / (2 * np.pi * self.freqs[0])
            nu_hi = min(nu_hi, 2.0)  # cap at a physically reasonable slowness
            self.bins_radial = _logbins(nu_lo, nu_hi, self.bins_per_octave)
        else:
            self.bins_radial = bins_radial

        # determine length of time series
        # if dataset contains more than one hour of data, it will be splitted
        # into `block_size` and the output will be time-dependent
        time_length = (
            (self.dataset["time"][-1] - self.dataset["time"][0]).item() / 3600e9
        )
        if time_length > 1.0:
            logger.warning(
                f"Length of time series is {time_length:.2f} hours."
                f"I understand that you want spectra every {self.block_size}. "
                f"Time series are being splitted."
            )
            groups = self.dataset.resample(time=self.block_size)
            results = (
                self.estimate_directional_distribution(subset)
                .compute(use=self.use)
                .expand_dims({"time": [time]})
                for time, subset in tqdm(groups, desc="Processing:")
                if len(subset["time"]) > 1
            )
            return xr.concat(results, dim="time")
        else:
            return self.estimate_directional_distribution(self.dataset)


class Triplets(_BaseClass):
    """Perform EWDM for triplet-based data such as wave buoys or ADCPs.

    Arguments:
        dataset (xr.Dataset): Dataset containing input data. Users may have
        different input variables depending on the kind of devide. For example,
        Typical GPS buoys deliver horizontal displacements or velocities. Other
        buoys only provide horizontal acceleration. ADPCs provide two
        dimensional components of horizonal velocities and echo-based sea
        surface elevation. Hence, the dataset should contain either
        displacements, velocities or accelerations. Sea surface elevation should
        be provided in all cases. The convention followed for variable names is:

            .. code-block:: python

                <xarray.Dataset>
                Dimensions:                 (time)
                Coordinates:
                  * time                    (time) datetime64[ns]
                Data variables:
                    eastward_displacement   (time) float32
                    northward_displacement  (time) float32
                    surface_elevation       (time) float32
                    eastward_velocity       (time) float32
                    northward_velocity      (time) float32
                    eastward_acceleration   (time) float32
                    northward_acceleration  (time) float32
                    eastward_slope          (time) float32
                    northward_slope         (time) float32
                Attributes: (1/1)
                    sampling_rate:           2.5

    Returns:
        xr.Dataset: Dataset containing produced directional spectra
        and directional spreading function.
    """


    @classmethod
    def from_numpy(
        cls,
        time: np.ndarray,
        surface_elevation: np.ndarray,
        eastward_displacement: np.ndarray = None,
        northward_displacement: np.ndarray = None,
        eastward_velocity: np.ndarray = None,
        northward_velocity: np.ndarray = None,
        eastward_acceleration: np.ndarray = None,
        northward_acceleration: np.ndarray = None,
        **kwargs
    ):
        """
        Create an instance of Triplets from numpy arrays.

        Arguments:
            time: time
            surface_elevation: Surface elevation array
            eastward_displacement: Eastward displacements
            northward_displacement: Northward displacements
            eastward_velocities: Eastward velocities
            northward_velocities: Northward velocities
            eastward_acceleration: Eastward accelerations
            northward_acceleration: Northward accelerations
            time: Time values.
        """
        
        # create dataset from numpy arrays
        data_vars = {
            'time': time,
            'surface_elevation': ("time", surface_elevation)
        }

        if eastward_displacement is not None:
            data_vars['eastward_displacement'] = ("time", eastward_displacement)
        
        if northward_displacement is not None:
            data_vars['northward_displacement'] = ("time", northward_displacement)
        
        if eastward_velocity is not None:
            data_vars['eastward_velocity'] = ("time", eastward_velocity)
        
        if northward_velocity is not None:
            data_vars['northward_velocity'] = ("time", northward_velocity)
        
        if eastward_acceleration is not None:
            data_vars['eastward_acceleration'] = ("time", eastward_acceleration)

        if northward_acceleration is not None:
            data_vars['northward_acceleration'] = ("time", northward_acceleration)
        
        # store sampling_rate as a global attribute
        attrs = {}
        if 'fs' in kwargs:
            attrs = {'sampling_rate': kwargs["fs"]}

        dataset = xr.Dataset(data_vars=data_vars, attrs=attrs)

        return cls(dataset, **kwargs)


    def compute_velocities(self):
        """Compute velocity componentes from displacements"""
        try:
            self.dataset["eastward_velocity"] = (
                self.dataset["eastward_displacement"]
                .differentiate("time", datetime_unit="s")
            )
            self.dataset["northward_velocity"] = (
                self.dataset["northward_displacement"]
                .differentiate("time", datetime_unit="s")
            )
        except KeyError:
            raise Exception(
                "`eastward_displacement` and `northward_displacement` "
                "are required to calculate velocity components."
            )

    def compute_accelerations(self):
        """Compute acceleration componentes from velocities"""
        try:
            self.dataset["eastward_acceleration"] = (
                self.dataset["eastward_velocity"]
                .differentiate("time", datetime_unit="s")
            )
            self.dataset["northward_acceleration"] = (
                self.dataset["northward_velocity"]
                .differentiate("time", datetime_unit="s")
            )
        except KeyError:
            raise Exception(
                "`eastward_velocity` and `northward_velocity` "
                "are required to calculate acceleration components. "
                "consider runing `self.compute_velocities()` first."
            )


    def estimate_wavelet_power(self, dataset) -> xr.DataArray:
        """Estimate the wavelet power of the surface elevation data.

        This method computes the continuous wavelet transform (CWT) of the
        surface elevation data in the dataset to estimate the wavelet power.

        Returns:
            np.ndarray: The wavelet power of the surface elevation data.

        Raises:
            Exception: If the 'surface_elevation' data is not available
            in the dataset.
        """
        if "surface_elevation" in dataset:
            data_std = dataset["surface_elevation"].std().item()
            Wzz = cwt(
                dataset["surface_elevation"],
                freqs=self.freqs, fs=self.fs,
            )
            power = np.abs(Wzz)**2
            if self.normalise:
                # rescale so that the integral of the spectrum equals the
                # surface-elevation variance: int <|W|^2>_t df -> sigma^2.
                # The factor is variance / wavelet-energy (sigma^2 / E), not
                # sigma / sqrt(E), which would only be exact when E == sigma^2.
                wavelet_energy = power.mean("time").integrate("frequency")
                return power * data_std**2 / wavelet_energy.item()
            else:
                return power

        else:
            raise Exception(
                "Local wavelet power cannot be computed because "
                "`surface_elevation` data is not available in "
                "the dataset."
            )


    def theta_from_displacements(self, dataset) -> xr.DataArray:
        """Compute local wave direction from wave displacements.

        This method calculates the local wave direction using the eastward and
        northward displacements along with the surface elevation from the
        dataset. This method is based on Peláez-Zapata et al (2024).

        Returns:
            xr.DataArray: Local wave direction in degrees.

        Raises:
            Exception: If `eastward_displacement` and `northward_displacement`
            are not available in the dataset.
        """
        try:
            Wxz = xwt(
                dataset["eastward_displacement"],
                dataset["surface_elevation"],
                freqs=self.freqs, fs=self.fs
            )
            Wyz = xwt(
                dataset["northward_displacement"],
                dataset["surface_elevation"],
                freqs=self.freqs, fs=self.fs
            )
            return RADTODEG * np.arctan2((1j*Wyz).real, (1j*Wxz).real)
        except KeyError:
            logger.exception(
                "Local wave direction cannot be computed from wave "
                "displacements because `eastward_displacement` and "
                "`northward_displacement` are not available in the "
                "dataset."
            )


    def theta_from_velocities(self, dataset) -> xr.DataArray:
        """Compute local wave direction from wave velocities.

        This method calculates the local wave direction using the eastward and
        northward velocities along with the surface elevation from the dataset.
        This method is based on Peláez-Zapata et al (2024).

        Returns:
            xr.DataArray: Local wave direction in degrees.

        Raises:
            Exception: If `eastward_displacement` and `northward_displacement`
            are not available in the dataset.
        """
        # if velocities dont exist in the dataset then
        try:
            Wxz = xwt(
                dataset["eastward_velocity"],
                dataset["surface_elevation"],
                freqs=self.freqs, fs=self.fs
            )
            Wyz = xwt(
                dataset["northward_velocity"],
                dataset["surface_elevation"],
                freqs=self.freqs, fs=self.fs
            )
            return RADTODEG * np.arctan2(Wyz.real, Wxz.real)
        except KeyError:
            logger.exception(
                "Local wave direction cannot be computed from wave "
                "velocities because `eastward_velocity` and "
                "`northward_velocity` are not available in the "
                "dataset."
            )


    def theta_from_accelerations(self, dataset) -> xr.DataArray:
        """Compute local wave direction from wave accelerations.

        This method calculates the local wave direction using the eastward and
        northward velocities along with the surface elevation from the dataset.
        This method is based on Peláez-Zapata et al (2024).

        Returns:
            xr.DataArray: Local wave direction in degrees.

        Raises:
            Exception: If `eastward_acceleration` and `northward_acceleration`
            are not available in the dataset.
        """
        # if accelerations exist in the dataset then
        try:
            Wxz = xwt(
                dataset["eastward_acceleration"],
                dataset["surface_elevation"],
                freqs=self.freqs, fs=self.fs
            )
            Wyz = xwt(
                dataset["northward_acceleration"],
                dataset["surface_elevation"],
                freqs=self.freqs, fs=self.fs
            )
            return RADTODEG * np.arctan2(Wyz.imag, Wxz.imag)
        except KeyError:
            logger.exception(
                "Local wave direction cannot be computed from wave "
                "accelerations because required variables are not "
                "available in the dataset."
            )


    def theta_from_slopes(self, dataset):
        """Compute local wave direction from wave slopes.

        This method calculates the local wave direction using the eastward and
        northward slopes, also known as roll and pitch, respectively, along
        with the surface elevation from the dataset. This method is based on
        Krogstad et al. (2005).

        Returns:
            xr.DataArray: Local wave direction in degrees.

        Raises:
            Exception: If `eastward_slope` and `northward_slope`
            are not available in the dataset.
        """
        # if accelerations exist in the dataset then
        try:
            Wxz = xwt(
                dataset["eastward_slope"],
                dataset["surface_elevation"],
                freqs=self.freqs, fs=self.fs
            )
            Wyz = xwt(
                dataset["northward_slope"],
                dataset["surface_elevation"],
                freqs=self.freqs, fs=self.fs
            )
            return RADTODEG * np.arctan2(Wyz.imag, Wxz.imag)
        except KeyError:
            raise Exception(
                "Local wave direction cannot be computed from wave "
                "slopes because required variables are not "
                "available in the dataset."
            )


    def estimate_directional_distribution(self, dataset) -> xr.Dataset:
        """Estimate the directional distribution of wave energy.

        This method calculates the wavelet power and local wave direction using
        the specified method (displacements, velocities, or accelerations) and
        then estimates the directional distribution function and directional
        spectra.

        Returns:
            xr.Dataset: Dataset containing the directional spectrum, directional
            distribution, and frequency spectrum.

        Raises:
            Exception: If `use` is not one of `displacements`, `velocities`,
            `accelerations` or `slopes`.
        """

        if self.interpolate:
            _dataset = self.interpolate_dataset(
                dataset, self.max_nan_ratio, self.max_time_gap
            )
        else:
            _dataset = dataset.copy()

        power = self.estimate_wavelet_power(_dataset)
        if self.use == "displacements":
            theta = self.theta_from_displacements(_dataset)
        elif self.use == "velocities":
            theta = self.theta_from_velocities(_dataset)
        elif self.use == "accelerations":
            theta = self.theta_from_accelerations(_dataset)
        elif self.use == "slopes":
            theta = self.theta_from_slopes(_dataset)
        else:
            raise Exception(
                "`use` should be either `displacements`, `velocities` "
                "`accelerations` or `slopes`."
            )
        return estimate_directional_distribution(
            power, theta, dd=self.dd, kappa=self.kappa
        )


    def compute(
            self,
            omin: float = -5,
            omax: float = None,
            nvoice: float = 16,
            dd: float = 5.0,
            kappa: float = 36.0,
            use: str = "displacements",
            block_size: str = "30min",
        ) -> xr.Dataset:
        """Perform computation using specified parameters.

        Args:
            omin (float, optional): Minimum octave (default is -5).
            omax (float, optional): Maximum octave. If None, it is
                automatically determined based sampling frequency. The final
                frequency array is logaritmically distributed from
                `2**omin` to `2**omax`.
            nvoice (float, optional): Number of voices for the computation
                (default is 16).
            dd (float, optional): Directional resolution in degrees
                (default is 5 degrees).
            kappa (float, optional): Smoothness parameter for Kernel
                Density Estimation. Small values of `kappa` produce oversmooth
                results (default is 36.0).
            use (str, optional): Type of data to perform estimation.
                It should be should be either `displacements`, `velocities`
                or `accelerations`." (default is "displacements").
            block_size (str): If dataset contains more than one hour of data,
                split dataset into blocks of `block_size` and perform
                computation over each block. The resulting output will have a
                time dimension. It is advisable to choose values of no more than
                half-hour. Default `block_size="30min"`.

        Returns:
            xr.Dataset: Dataset containing the directional spectrum, directional
                distribution, and frequency spectrum.
        """

        # the maximum frequency is given by the nyquist frequency
        if omax is None:
            omax = int(np.log2(self.fs / 2))

        self.omin = omin
        self.omax = omax
        self.nvoice = nvoice

        # fourier equivalent frequencies
        self.freqs = 2. ** np.linspace(omin, omax, nvoice*abs(omin-omax)+1)

        # directional resolution
        self.dd = dd
        self.kappa = kappa

        # data used for estimation
        self.use = use
        self.block_size = block_size

        # determine length of time series
        # if dataset contains more than one hour of data, it will be splitted
        # into `block_size` and the output will be time-dependent
        time_length = (
            (self.dataset["time"][-1] - self.dataset["time"][0]).item() / 3600e9
        )
        if time_length > 1.0:
            logger.warning(
                f"Length of time series is {time_length:.2f} hours."
                f"I understand that you want spectra every {self.block_size}. "
                f"Time series are being splitted."
            )
            groups = self.dataset.resample(time=self.block_size)
            results = (
                self.estimate_directional_distribution(subset)
                .compute(use=self.use)
                .expand_dims({"time": [time]})
                for time, subset in tqdm(groups, desc="Processing:")
                if len(subset["time"]) > 1
            )
            return xr.concat(results, dim="time")
        else:
            return self.estimate_directional_distribution(self.dataset)

# }}}



if __name__ == "__main__":
    pass

# --- end of file ---

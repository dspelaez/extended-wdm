#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

import numpy as np
import xarray as xr
import logging

from tqdm import tqdm

from .wavelets import cwt, xwt
from .density import estimate_directional_distribution
from .helpers import get_sampling_frequency
from .sources import SpotterBuoysDataSource, CDIPDataSourceRealTime
from .parameters import VARIABLE_NAMES


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
            except AttributeError:
                self.fs = get_sampling_frequency(self.dataset["time"])
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
            return dataset * np.nan
        else:
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
            Dimensions:                 (x, y, time)
            Coordinates:
              * time                    (time) datetime64[ns]
              * x                       (time) datetime64[ns]
              * y                       (time) datetime64[ns]
            Data variables:
                surface_elevation       (x, y, time) float32
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
        position_x: np.ndarray,
        position_y: np.ndarray
        ):
        """
        Create an instance of Arrays from numpy arrays.

        Arguments:
            surface_elevation: Surface elevation array
            eastward_displacement: Eastward displacements
            northward_displacement: Northward displacements
            eastward_velocities: Eastward velocities
            northward_velocities: Northward velocities
            eastward_acceleration: Eastward accelerations
            northward_acceleration: Northward accelerations
            time: Time values.
        """
        pass

    
    def wavelet_coefficients(self, dataset: xr.Dataset) -> xr.Dataset:
        ...

    def array_geometry(self, dataset: xr.Dataset) -> np.ndarray:
        ...
    
    def compute_phase(self, dataset: xr.Dataset) -> xr.Dataset:
        ...

    def compute_power(self, dataset: xr.Dataset) -> xr.Dataset:
        ...




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
    ):
        """
        Create an instance of Arrays from numpy arrays.

        Arguments:
            surface_elevation: Surface elevation array
            eastward_displacement: Eastward displacements
            northward_displacement: Northward displacements
            eastward_velocities: Eastward velocities
            northward_velocities: Northward velocities
            eastward_acceleration: Eastward accelerations
            northward_acceleration: Northward accelerations
            time: Time values.
        """
        pass


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
                wavelet_energy = power.mean("time").integrate("frequency")**0.5
                return power * data_std / wavelet_energy.item()
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
        if all(
            var in dataset for var in
                ["eastward_displacement", "northward_displacement"]
            ):
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
        else:
            raise Exception(
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
            raise Exception(
                "Local wave direction cannot be computed from wave "
                "velocities because `eastward_velocity` and "
                "`northward_velocity` are not available in the "
                "dataset.\n"
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
            raise Exception(
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
            groups = self.dataset.resample(time=self.block_size)
            results = (
                self.estimate_directional_distribution(subset)
                .compute(use=self.use)
                .expand_dims({"time": [time]})
                for time, subset in tqdm(groups)
                if len(subset["time"]) > 1
            )
            return xr.concat(results, dim="time")
        else:
            return self.estimate_directional_distribution(self.dataset)

# }}}



if __name__ == "__main__":
    pass

# --- end of file ---

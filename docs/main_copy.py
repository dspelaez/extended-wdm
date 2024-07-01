#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

import numpy as np
import xarray as xr
import logging

from .wavelets import cwt, xwt
from .helpers import get_sampling_frequency
from .sources import SpotterBuoysDataSource, CDIPDataSourceRealTime
from .parameters import VARIABLE_NAMES


RADTODEG = 180. / np.pi
DEGTORAD = np.pi / 180.


# von mises kernel desnsity estimation
def _vonmises_kde(arr, bins, kappa):
    """Return Kernel-Density Estimation using von Mises distribution"""

    # define input parameters
    x = np.radians(bins[:,None] - arr[None,:])

    # integrate vonmises kernels
    kde = (
        np.exp(kappa * np.cos(x)).sum(axis=1) / (2 * np.pi * np.i0(kappa))
    )
    kde /= np.trapz(kde, x=bins)
    return kde


# function to get histogram along freq axis
def _get_density(arr, bins, kappa):
    if np.isnan(arr).all():
        return np.zeros_like(bins, dtype="float") * np.nan
    else:
        if kappa is not None:
            return _vonmises_kde(arr, bins=bins, kappa=kappa)
        else:
            bins_edges = np.r_[bins, bins[-1]+np.diff(bins)[0]]
            return np.histogram(arr, bins=bins_edges, density=True)[0]


# function to actually estimate the spectrum
def _estimate_directional_distribution(power, theta, dd, kappa):
    """Construct directional distribution function from local wave directions"""

    # array of directiontions where dd is the directional resolution
    bins = np.arange(-180, 180, dd)

    # directional distribution function
    D = np.apply_along_axis(
        _get_density, arr=theta, bins=bins, kappa=kappa, axis=1
    )

    # determine average wavelet power
    S = power.mean("time").data

    # array containing directional wave spectra
    E = S[:,None] * D

    # return dataset
    output_ds = xr.Dataset(
        data_vars = {
            "directional_spectrum": (["frequency", "direction"], E),
            "directional_distribution": (["frequency", "direction"], D),
            "frequency_spectrum": (["frequency"], S)
        },
        coords = {
            "frequency": power["frequency"].data,
            "direction": bins,
        }
    )
    output_ds["frequency"].attrs = VARIABLE_NAMES["frequency"]
    output_ds["direction"].attrs = VARIABLE_NAMES["direction"]
    for var in output_ds:
        output_ds[var].attrs = VARIABLE_NAMES[var]

    return output_ds


class ExtendedWDM(object):
    ...

class WaveStaffsExtendedWDM(object):
    ...

class PitchRollBouysExtendedWDM(object):
    ...

class AccelBouysExtendedWDM(object):
    ...

class GPSBouysExtendedWDM(object):
    ...

class AdcpExtendedWDM(object):
    ...

# class to compute directional spectrum from wave buoy data
class BuoysEWDM(object):
    """Perform EWDM for wave buoys.

    Args:
        dataset (xr.Dataset): Dataset containing buoy data. Users may have different
            input variables depending on the kind of buoy. Typical GPS buoys deliver
            horizontal displacements or velocities. Other buoys only provide acceleration.
            Hence, the dataset should contain either displacements, velocities or
            accelerations. Sea surface elevation should be provided in all cases.
            The convention followed for variable names is:

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
                Attributes: (1/1)
                    sampling_rate:           2.5

    Returns:
        xr.Dataset: Dataset containing produced directional spectra and directional
        spreading function.
    """

    def __init__(
            self,
            dataset: xr.Dataset,
            fs: float = None,
            clean_dataset: bool = True,
            max_nan_ratio: float = 0.1,
            max_time_gap: str = "10s",
        ) -> xr.Dataset:
        """Initialise BuoysEWDM class"""

        # TODO: Pitch-roll has not been tested yet.

        # interpolate dataset to remove nan variables
        ntime = len(dataset["time"])
        if clean_dataset:
            nan_ratio = (dataset["surface_elevation"].isnull().sum() / ntime)
            if nan_ratio > max_nan_ratio:
                self.dataset = dataset * np.nan
            else:
                self.dataset = dataset.interpolate_na(
                    dim="time", method="quadratic", max_gap=max_time_gap,
                    fill_value="extrapolate"
                )
        else:
            self.dataset = dataset

        if fs is None:
            try:
                self.fs = self.dataset.sampling_rate
            except AttributeError:
                self.fs = get_sampling_frequency(self.dataset["time"])
        else:
            self.fs = fs


    @classmethod
    def from_numpy(cls,
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
        Create an instance of BuoyEWDM from numpy arrays.

        Args:
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



    def estimate_wavelet_power(self):
        """Estimate the wavelet power of the surface elevation data.

        This method computes the continuous wavelet transform (CWT) of the surface
        elevation data in the dataset to estimate the wavelet power.

        Returns:
            np.ndarray: The wavelet power of the surface elevation data.

        Raises:
            Exception: If the 'surface_elevation' data is not available in the dataset.
        """
        if "surface_elevation" in self.dataset:
            return cwt(
                self.dataset["surface_elevation"],
                freqs=self.freqs, fs=self.fs
            )
        else:
            raise Exception(
                "Local wavelet power cannot be computed because "
                "`surface_elevation` data is not available in "
                "the dataset."
            )


    def theta_from_displacements(self) -> xr.DataArray:
        """Compute local wave direction from wave displacements.

        This method calculates the local wave direction using the eastward and northward
        displacements along with the surface elevation from the dataset.

        Returns:
            xr.DataArray: Local wave direction in degrees.

        Raises:
            Exception: If `eastward_displacement` and `northward_displacement` are not
            available in the dataset.
        """
        if all(
            var in self.dataset for var in
                ["eastward_displacement", "northward_displacement"]
            ):
            Wxz = xwt(
                self.dataset["eastward_displacement"],
                self.dataset["surface_elevation"],
                freqs=self.freqs, fs=self.fs
            )
            Wyz = xwt(
                self.dataset["northward_displacement"],
                self.dataset["surface_elevation"],
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


    def theta_from_velocities(self) -> xr.DataArray:
        """Compute local wave direction from wave velocities.

        This method calculates the local wave direction using the eastward and northward
        velocities along with the surface elevation from the dataset. If velocities are
        not available in the dataset, it tries to compute them from displacements.

        Returns:
            xr.DataArray: Local wave direction in degrees.

        Raises:
            Exception: If `eastward_displacement` and `northward_displacement` are not
            available in the dataset.
        """
        # if velocities dont exist in the dataset then
        if not all(
            var in self.dataset for var in
                ["eastward_velocity", "northward_velocity"]
            ):

            # try to compute velocities from displacements
            try:
                self.dataset["eastward_velocity"] = (
                    self.dataset["eastward_displacement"]
                    .differentiate("time", datetime_unit="s")
                )
                self.dataset["northward_velocity"] = (
                    self.dataset["northward_displacement"]
                    .differentiate("time", datetime_unit="s")
                )
            # if not displacements available then raise error
            except KeyError:
                raise Exception(
                    "Local wave direction cannot be computed from wave "
                    "velocities because `eastward_displacement` and "
                    "`northward_displacement` are not available in the "
                    "dataset."
                )

        # if we do have velocities then compute local wave direction
        Wxz = xwt(
            self.dataset["eastward_velocity"],
            self.dataset["surface_elevation"],
            freqs=self.freqs, fs=self.fs
        )
        Wyz = xwt(
            self.dataset["northward_velocity"],
            self.dataset["surface_elevation"],
            freqs=self.freqs, fs=self.fs
        )

        return RADTODEG * np.arctan2(Wyz.real, Wxz.real)


    def theta_from_accelerations(self) -> xr.DataArray:
        """Compute local wave direction from wave accelerations.

        This method calculates the local wave direction using the eastward and northward
        velocities along with the surface elevation from the dataset. If accelerations
        or are not available in the dataset, it tries to compute them from displacements
        or velocities.

        Returns:
            xr.DataArray: Local wave direction in degrees.

        Raises:
            Exception: If `eastward_displacement` and `northward_displacement` are not
            available in the dataset.
        """
        # if accelerations dont exist in the dataset then
        if not all(
            var in self.dataset for var in
                ["eastward_acceleration", "northward_acceleration"]
            ):
            #
            # try to compute accelerations from velocities
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
                pass

            # if velocities not found, compute them from displacements
            try:
                # compute velocities first
                self.dataset["eastward_velocity"] = (
                    self.dataset["eastward_displacement"]
                    .differentiate("time", datetime_unit="s")
                )
                self.dataset["northward_velocity"] = (
                    self.dataset["northward_displacement"]
                    .differentiate("time", datetime_unit="s")
                )
                # then compute accelerations
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
                    "Local wave direction cannot be computed from wave "
                    "accelerations because required variables are not "
                    "available in the dataset."
                )

        Wxz = xwt(
            self.dataset["eastward_acceleration"],
            self.dataset["surface_elevation"],
            freqs=self.freqs, fs=self.fs
        )
        Wyz = xwt(
            self.dataset["northward_acceleration"],
            self.dataset["surface_elevation"],
            freqs=self.freqs, fs=self.fs
        )
        return RADTODEG * np.arctan2(Wyz.imag, Wxz.imag)


    def estimate_directional_distribution(self) -> xr.Dataset:
        """Estimate the directional distribution of wave energy.

        This method calculates the wavelet power and local wave direction using the
        specified method (displacements, velocities, or accelerations) and then
        estimates the directional distribution function and directional spectra.

        Returns:
            xr.Dataset: Dataset containing the directional spectrum, directional
            distribution, and frequency spectrum.

        Raises:
            Exception: If `use` is not one of `displacements`, `velocities`,
            or `accelerations`.
        """
        self.power = np.abs(self.estimate_wavelet_power()) ** 2
        if self.use == "displacements":
            self.theta = self.theta_from_displacements()
        elif self.use == "velocities":
            self.theta = self.theta_from_velocities()
        elif self.use == "accelerations":
            self.theta = self.theta_from_accelerations()
        else:
            raise Exception(
                "`use` should be either `displacements`, `velocities` or `accelerations`."
            )
        return _estimate_directional_distribution(
            self.power, self.theta, dd=self.dd, kappa=self.kappa
        )


    def compute(
            self,
            omin: float = -5,
            omax: float = None,
            nvoice: float = 16,
            dd: float = 5.0,
            kappa: float = 36.0,
            use: str = "displacements",
        ) -> xr.Dataset:
        """Perform computation using specified parameters.

        Args:
            omin (float, optional): Minimum octave (default is -5).
            omax (float, optional): Maximum octave. If None, it is automatically determined
                based sampling frequency. The final frequency array is logaritmically
                distributed from `2**omin` to `2**omax`.
            nvoice (float, optional): Number of voices for the computation (default is 16).
            dd (float, optional): Directional resolution in degrees (default is 5 degrees).
            kappa (float, optional): Smoothness parameter for Kernel Density Estimation.
                Small values of `kappa` produce oversmooth results (default is 36.0).
            use (str, optional): Type of data to perform estimation. It should be should be
                either `displacements`, `velocities` or `accelerations`." (default is "displacements").

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


        # determine length of time series
        # if dataset contains more than one hour of data, it will be splitted
        # into `block_size` and the output will be time-dependent
        time_length = (
            (self.dataset["time"][-1] - self.dataset["time"][0]).item() / 3600e9
        )
        if time_length > 1.0:
            print("Warning")

        return self.estimate_directional_distribution()

# }}}



if __name__ == "__main__":
    pass

# --- end of file ---

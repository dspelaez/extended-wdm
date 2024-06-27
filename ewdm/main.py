#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

"""
@author: Daniel Peláez-Zapata
@github: http://github.com/dspelaez
@created: 2024-03-12
"""

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


# class to compute directional spectrum from wave buoy data
class BuoysEWDM(object):

    def __init__(
            self,
            dataset: xr.Dataset,
            fs: float = None,
            clean_dataset: bool = True,
            max_nan_ratio: float = 0.1,
            max_time_gap: str = "10s",
        ) -> xr.Dataset:
        """Perform EWDM for wave buoys

        Arguments:
        - dataset (xr.Dataset): Dataset contaning buoy data. Users may have
          different input variables depending on the kind of buoy. Typical GPS
          buoys deliver horizontal displacements or velocities. Other buoys only
          provide acceleration. Hence, the dataset should containe either
          dispacements, velocities or accelerations. Sea surface elevation
          should be provided in all cases. The convention followed for variable
          names is:

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
        - output (xr.Dataset): Dataset contaning produced directional spectra
          and directional spreading function.
        """

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
        surface_elevation: np.ndarray,
        eastward_displacement: np.ndarray = None,
        northward_displacement: np.ndarray = None,
        eastward_velocity: np.ndarray = None,
        northward_velocity: np.ndarray = None,
        eastward_acceleration: np.ndarray = None,
        northward_acceleration: np.ndarray = None,
        time: np.ndarray = None
    ):
        """
        Create an instance of BuoyEWDM from numpy arrays.

        Args:
        - surface_elevation (np.ndarray): Surface elevation array
        - eastward_displacement (np.ndarray, optionl): Eastward displacements
        - northward_displacement (np.ndarray, optional): Northward displacements
        - eastward_velocities (np.ndarray, optional): Eastward velocities
        - northward_velocities (np.ndarray, optional): Northward velocities
        - eastward_acceleration (np.ndarray, optional): Eastward accelerations
        - northward_acceleration (np.ndarray, optional): Northward accelerations
        - time (np.ndarray, optional): Time values (optional).

        """
        pass

    
    
    def estimatee_wavelet_power(self):
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


    def theta_from_displacements(self):
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


    def theta_from_velocities(self):
        
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


    def theta_from_accelerations(self):

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


    def estimate_directional_distribution(self):

        self.power = np.abs(self.estimate_wavelet_power())**2
        if self.use == "displacements":
            self.theta = self.theta_from_displacements()
        elif self.use == "velocities":
            self.theta = self.theta_from_velocities()
        elif use == "accelerations":
            self.theta = self.theta_from_accelerations()
        else:
            raise Exception(
                "`estimate_from` should be either `displacements`, "
                "`velocities` or `accelerations`"
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
            **kwargs
        ):
        """Compute"""

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

    from plots import plot_directional_spectrum
    from matplotlib import pyplot as plt
    plt.ion()
    dataset = (
        CDIPDataSourceRealTime(189)
        .read_dataset(time_start='2024-06-24T08:00')
    )

    # dataset = (
        # CDIPDataSourceRealTime(188)
        # .read_dataset(time_start='2024-06-22T06:00')
    # )

    # spotter = SpotterBuoysDataSource("../data/displacement.csv")
    # dataset = spotter.read_dataset()

    # groups = dataset.resample(time="60min")
    # results = (
        # BuoysEWDM(subset)
        # .compute(use="displacements")
        # .expand_dims({"time": [time]})
        # for time, subset in groups if len(subset["time"]) > 1
    # )
    # a1 = xr.concat(results, dim="time")

    # for i in range(0,24,3):
        # print(i)
        # plot_directional_spectrum(
            # a1.isel(time=i).directional_distribution, dirs="direction", frqs="frequency",
            # levels=None, colorbar=True, axes_kw={"rmax": 0.5, "is_period": True}
        # )
        
    spec = BuoysEWDM(dataset)
    a1=spec.compute(use="displacements")
    # a2=spec.compute(use="velocities")
    # a3=spec.compute(use="accelerations")

    plot_directional_spectrum(
        a1.directional_distribution, dirs="direction", frqs="frequency",
        levels=None, colorbar=True, axes_kw={"rmax": 0.5, "is_period": True}
    )

    # plot_directional_spectrum(
        # a2.directional_distribution, dirs="direction", frqs="frequency",
        # levels=None, colorbar=True, axes_kw={"rmax": 0.5, "is_period": True}
    # )

    # plot_directional_spectrum(
        # a3.directional_distribution, dirs="direction", frqs="frequency",
        # levels=None, colorbar=True, axes_kw={"rmax": 0.5, "is_period": True}
    # )

# --- end of file ---

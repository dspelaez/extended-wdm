#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8


import numpy as np
import xarray as xr
from scipy.signal import welch
from matplotlib import pyplot as plt

import ewdm
from ewdm.sources import CDIPDataSourceRealTime, SpotterBuoysDataSource
from ewdm.plots import plot_directional_spectrum
plt.ion()


if True:
    spotter = SpotterBuoysDataSource("../data/displacements.csv")
    dataset = spotter.read_dataset()
    spec = ewdm.Triplets(dataset)
    output = spec.compute(
        omin=-5, omax=0, nvoice=16, dd=5, kappa=36, use="displacements"
    )
    output = spec.compute()

    fig, ax = plt.subplots()
    plot_directional_spectrum(
        output.directional_spectrum, ax=ax,
        levels=None, colorbar=False,
        axes_kw={"rmax": 1, "as_period": False}
    )
    plt.savefig("spotter-example-directional-spectrum.png")


if False:
    cdip =  CDIPDataSourceRealTime(166)
    dataset = cdip.read_dataset(time_start='2024-06-09T08:30')
    spec = ewdm.Triplets(dataset)
    output = spec.compute(
        omin=-5, omax=0, nvoice=16, dd=5, kappa=36, use="displacements"
    )

    fig, ax = plt.subplots()
    plot_directional_spectrum(
        output.directional_spectrum, ax=ax,
        levels=None, colorbar=False, contours=3,
    )
    plt.savefig("cdip-example-directional-spectrum.png")

    fourier_freq, fourier_power = welch(
        dataset["surface_elevation"], fs=dataset.sampling_rate, nperseg=512
    )
    fig, ax = plt.subplots()
    ax.loglog(fourier_freq, fourier_power, label="Fourier-based")
    output.frequency_spectrum.plot(ax=ax, ls="--", label="Wavelet-based")
    ax.legend()


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

# a2=spec.compute(use="velocities")
# a3=spec.compute(use="accelerations")

# plot_directional_spectrum(
    # a2.directional_distribution, dirs="direction", frqs="frequency",
    # levels=None, colorbar=True, axes_kw={"rmax": 0.5, "is_period": True}
# )

# plot_directional_spectrum(
    # a3.directional_distribution, dirs="direction", frqs="frequency",
    # levels=None, colorbar=True, axes_kw={"rmax": 0.5, "is_period": True}
# )

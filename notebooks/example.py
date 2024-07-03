#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8


import numpy as np
import xarray as xr
import os

from scipy.signal import welch
from matplotlib import pyplot as plt

import ewdm
from ewdm.sources import CDIPDataSourceRealTime, SpotterBuoysDataSource
from ewdm.plots import plot_directional_spectrum

plt.ion()

# HERE = os.path.dirname(os.path.abspath(__file__))
# SPOTTER_FILE = os.path.join(HERE, "../data/displacement.csv")

if False:
    spotter = SpotterBuoysDataSource("../data/displacement.csv")
    dataset = spotter.read_dataset()
    spec = ewdm.Triplets(dataset)
    spec.compute_velocities()
    output = spec.compute(
        omin=-5, omax=0, nvoice=16, dd=5, kappa=36,
        block_size="60min", use="displacements"
    )
    output = spec.compute()

    fig, axs = plt.subplots(2, 3, figsize=(9,6), layout="tight")
    axs = axs.ravel()
    for ax, tt in zip(axs, output.time):
        plot_directional_spectrum(
            output.sel(time=tt).directional_distribution,
            ax=ax, vmin=0.003, vmax=0.015, levels=None, colorbar=False,
            axes_kw={"rmax": 0.6, "as_period": False}
        )
        ax.set_title(tt.dt.strftime("%Y-%m-%d %H:%M:%S").item())
    plt.savefig("spotter-example-directional-spectrum.png")


if True:
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
    output.frequency_spectrum.plot(ax=ax, ls="--", lw=2, label="Wavelet-based")
    ax.legend()
    plt.savefig("cdip-fourier-vs-wavelet-avg-power.png")

#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

"""
Fourier vs. wavelet power density
=================================

A quick comparison between azimutally-integrated wavelet power and Fourier-based
power density spectrum.
"""

import numpy as np
import xarray as xr

from scipy.signal import welch
from matplotlib import pyplot as plt

import ewdm
from ewdm.sources import SpotterBuoysDataSource

spotter = SpotterBuoysDataSource("../../data/displacement.csv")
dataset = spotter.read_dataset()

subset = dataset.sel(time=slice("2020-08-20 10:00", "2020-08-20 10:30"))

spec = ewdm.Triplets(subset)
output = spec.compute()

fourier_freq, fourier_power = welch(
    subset["surface_elevation"], fs=dataset.sampling_rate, nperseg=512
)
fig, ax = plt.subplots()
ax.loglog(fourier_freq, fourier_power, label="Fourier-based")
output.frequency_spectrum.plot(ax=ax, ls="--", lw=2, label="Wavelet-based")
ax.legend()

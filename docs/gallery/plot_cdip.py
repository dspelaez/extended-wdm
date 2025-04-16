#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

"""
Using real time CDIP data
=========================

This example shows how to compute the directional wave spectrum from waverider wave buoys data obtained from CDIP dataset in real time.
"""

# import the dependencies
import numpy as np
import xarray as xr

from matplotlib import pyplot as plt

import ewdm
from ewdm.sources import CDIPDataSourceRealTime
from ewdm.plots import plot_directional_spectrum


# %%
# Let's download 30-min time series from Ocean Station Papa Buoy - 166p1. In
# this case, the buoy captures the wave-induced hortizonal $x(t)$, $y(t)$, 
# and vertical displacements $\eta(t)$.
cdip =  CDIPDataSourceRealTime(166)
dataset = cdip.read_dataset(
    time_start="2024-06-09T08:30", time_end="2024-06-09T09:00"
)
print(dataset)

# plot time series of the first 5 minutes
fig, (ax1, ax2, ax3) = plt.subplots(
    3,1, figsize=(6,6), sharex=True, sharey=True
)
subset = dataset.sel(time=slice("2024-06-09T08:30", "2024-06-09T08:35"))
subset["surface_elevation"].plot(ax=ax1)
subset["eastward_displacement"].plot(ax=ax2)
subset["northward_displacement"].plot(ax=ax3)
_ = ax1.set(xlabel="", ylabel="$\\eta(t)$ [m]")
_ = ax2.set(xlabel="", ylabel="$x(t)$ [m]")
_ = ax3.set(xlabel="", ylabel="$y(t)$ [m]")

# %%
# Now, let's compute the directional spectra using the `ewdm.Triplet` class. The
# input parameters `omin`, `omax` and `nvoice` control the frequency resolution,
# it means that a frequency array ranging from 2^-5 to 2^0 with 16 subvisions
# per octave will be created in this case. The parameters `dd` and `kappa`
# control the directional resolution and smoothness, respectively.
spec = ewdm.Triplets(dataset)
output = spec.compute(
    omin=-5, omax=0, nvoice=16, dd=5, kappa=36, use="displacements"
)
print(output)

# plot the results
fig, ax = plt.subplots(figsize=(6,5))
plot_directional_spectrum(
    output.directional_spectrum, ax=ax, levels=None, colorbar=True,
    cbar_kw={"label": "$E(f,\\theta)\\,\\mathrm{[m^2 Hz^{-1} deg^{-1}]}$"}
)

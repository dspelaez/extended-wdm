#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

"""
Using numpy arrays as input
===========================

You don't have a `xarray.Dataset` ready to use as input? No worries, you can use
the `.from_numpy()` class method to initialise the `Triplets` and `Arrays` classes from `numpy.ndarray` objects
"""

import numpy as np
import xarray as xr

from scipy.io import loadmat
from matplotlib import pyplot as plt

import ewdm
from ewdm.sources import SpotterBuoysDataSource
from ewdm.plots import plot_directional_spectrum


# %%
# Testing with ewdm.Triplets
# --------------------------
#
# Let's see how this feature work with the ewdm.Triplets class

spotter = SpotterBuoysDataSource("../../data/displacement.csv")
dataset = spotter.read_dataset()

subset = dataset.sel(time=slice("2020-08-20 10:00", "2020-08-20 10:30"))

# extract numpy arrays from dataset
time = subset["time"].data
surface_elevation = subset["surface_elevation"].data
eastward_displacement = subset["eastward_displacement"].data
northward_displacement = subset["northward_displacement"].data
sampling_rate = subset.attrs["sampling_rate"]


# %%
# Now, we can initialise the `ewdm.Triplets` class with the numpy arrays
spec = ewdm.Triplets.from_numpy(
    time = time,
    surface_elevation = surface_elevation,
    eastward_displacement = eastward_displacement,
    northward_displacement = northward_displacement,
    fs=sampling_rate,
    max_time_gap="30s"
)
print("Sampling rate is: ", spec.fs)
print("Maximum allowed time gap is: ", spec.max_time_gap)

# %%
# Let's compute the spectra and print the output
output_triplets = spec.compute()
print(output_triplets)


# %%
# Testing with ewdm.Arrays
# ------------------------
#
# It works for `ewdm.Arrays` class too.

# loading matlab file
mat_fname = "../../data/donelan_run62.mat"
mat_data = loadmat(mat_fname, simplify_cells=True)

sampling_rate = mat_data["fs"]
time = np.arange(len(mat_data["eta"])) / sampling_rate
elements = np.arange(len(mat_data["x"]))

print("Surface elevation shape: ", mat_data["eta"].shape)
print("Time array shape: ", time.shape)
print("Position x shape : ", mat_data["x"].shape)
print("Position y shape : ", mat_data["y"].shape)


# %%
# Now let's create spec object

spec = ewdm.Arrays.from_numpy(
    time = time,
    surface_elevation = mat_data["eta"],
    position_x = mat_data["x"],
    position_y = mat_data["y"],
    fs = sampling_rate,
    max_nan_ratio = 0.05,
)

print("Sampling rate is: ", spec.fs)
print("Maximum allowed nan ratio: ", spec.max_nan_ratio)

# %%
# Let's compute the spectra and print the output
output_array = spec.compute()
print(output_array)


# %%
# Let's plot the results
# ----------------------
#
# It's always nice to see some specral plots, isn't it?

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,4), layout="constrained")

plot_directional_spectrum(
    output_array.directional_spectrum, ax=ax1, levels=None, colorbar=True,
    axes_kw={"rmin": 0.1, "rmax": 1.2, "rstep": 0.2, "angle": 135},
    cbar_kw={"label": ""}
)

plot_directional_spectrum(
    output_triplets.directional_spectrum, ax=ax2, levels=None, colorbar=True,
    axes_kw={"rmin": 0.1, "rmax": 1.2, "rstep": 0.2, "angle": 135},
    cbar_kw={"label": "$E(f,\\theta)$"}
)

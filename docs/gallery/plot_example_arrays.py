#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

"""
Donelan `run82` test case
=========================

This example shows `ewdm.Arrays` in action. We use the run62 as test case. More
information about this test case can be found in Donelan et al. 2015
"""

import numpy as np
import xarray as xr

from scipy.io import loadmat
from matplotlib import pyplot as plt

import ewdm
from ewdm.plots import plot_directional_spectrum

# loading matlab file
mat_fname = "../../data/donelan_run82.mat"
mat_data = loadmat(mat_fname, simplify_cells=True)

sampling_rate = mat_data["fs"]
time = np.arange(len(mat_data["eta"])) / sampling_rate
elements = np.arange(len(mat_data["x"]))

# creating dataset
dataset = xr.Dataset(
    data_vars = {
        "surface_elevation": (["time", "element"], mat_data["eta"]),
        "position_x": ("element", mat_data["x"]),
        "position_y": ("element", mat_data["y"])
    },
    coords = {"time": time, "element": elements},
    attrs = {"sampling_rate": sampling_rate}
)

# compute directional spectra
spec = ewdm.Arrays(dataset)
output = spec.compute(
    cross_wavelet=True, solver="least-squares", kappa=36
)

fig, ax = plt.subplots()
plot_directional_spectrum(
    output.directional_spectrum, ax=ax, levels=None, colorbar=True,
    axes_kw={"rmin": 0.1, "rmax": 1.2, "rstep": 0.2, "angle": 135},
    cbar_kw={"label": "$E(f,\\theta)$"}
)

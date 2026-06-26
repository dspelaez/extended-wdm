#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

"""
Multi-aperture spectra from stereo-video
========================================
The companion stereo-imaging example picks a handful of pixels from the field of
view and feeds their elevation time series to `ewdm.Arrays`, just as one would
with a real wave-staff array. Here we instead hand the *whole* dense elevation
field eta(y, x, t) to `ewdm.MultiApertureArrays`.

Rather than a single fixed array, the estimator seeds virtual staffs into a
ladder of nested apertures: coarse, widely spaced clusters for the long waves and
tight clusters for the short waves. Each aperture is only trusted within its own
anti-alias band and the per-aperture spectra are stitched together. From a single
field this yields not only the frequency-direction spectrum but also the
wavenumber spectrum F(k) and the inverse-phase-speed spectrum Q(nu), measured
directly from the spatial phase gradients without assuming a dispersion relation.

This example uses the same `Guimaraes et al (2020)`_ Black Sea data as the
stereo-imaging example, and compares the result against `ewdm.Arrays` on a
sparse pentagon of pixels.

.. _Guimaraes et al (2020): https://doi.org/10.1038/s41597-020-0492-9
"""

import os
import subprocess

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

import ewdm
from ewdm.plots import plot_directional_spectrum

# define paths and filenames (shared with the stereo-imaging example)
URL = "https://data-dataref.ifremer.fr/stereo/BS_2013/"
VIDEO_URL = URL + "2013-09-22_13-00-01_10Hz/nc/Surfaces_20130922_130001_short.nc"
CACHE_DIR = "../../data/"
LOCAL_VIDEO_FILE = os.path.join(CACHE_DIR, "stereo-video.nc")

# Black Sea acquisition parameters
SAMPLING_RATE = 10.0        # Hz
DEPTH = 30.0                # m (deep water for these short wind waves)
NFRAMES = 4096              # subset of frames


# %%
# Downloading the stereo-video dataset
# ------------------------------------
#
# We download the netCDF4 file with the stereo-imaging surfaces and place it in
# the data folder. This might take a few minutes the first time.

if not os.path.exists(LOCAL_VIDEO_FILE):
    print("Downloading the stereo-video file. It might take a few minutes.")
    subprocess.call(f"wget {VIDEO_URL} -O {LOCAL_VIDEO_FILE}", shell=True)
else:
    print("File already exists. Skipping download.")

nc_obj = nc.Dataset(LOCAL_VIDEO_FILE)


# %%
# Building the elevation field
# ----------------------------
#
# The dataset stores the gridded sea surface elevation `Z(time, y, x)` together
# with the pixel coordinates `X` and `Y` in metres. We read a subset of the
# frames and arrange the field as (ny, nx, time), the layout expected by
# `MultiApertureArrays.from_field`.

X = np.asarray(nc_obj["X"][:], float)
Y = np.asarray(nc_obj["Y"][:], float)
dx = (np.nanmax(X) - np.nanmin(X)) / (X.shape[1] - 1)

Z = np.asarray(nc_obj["Z"][:NFRAMES, :, :], float)   # (time, ny, nx)
field = np.transpose(Z, (1, 2, 0))                   # (ny, nx, time)
del Z

print(f"field (ny, nx, time) = {field.shape}")
print(f"dx = {dx:.3f} m   field of view = {field.shape[1] * dx:.1f} m")


# %%
# Computing the multi-aperture spectra
# ------------------------------------
#
# The field of view is only a few wavelengths wide, so this is the small-aperture
# regime the default ladder is tuned for. We keep the anti-alias gate off (the
# default) since here the dominant wave is near the field-of-view scale.

spec = ewdm.MultiApertureArrays.from_field(field, dx, depth=DEPTH,
                                           fs=SAMPLING_RATE)
output = spec.compute(n_staff=16, seed=20, return_apertures=True)
print(output)

Hs = 4 * np.sqrt(float(output["var_eta"]))
f = output["frequency"].values
fp = f[np.nanargmax(output["frequency_spectrum"].values)]
print(f"Hs = {Hs:.2f} m   Tp = {1 / fp:.1f} s")


# %%
# Comparing against a sparse array
# --------------------------------
#
# As a reference we run `ewdm.Arrays` on a pentagon of pixels, exactly as in the
# stereo-imaging example. The two frequency spectra peak at the same period; the
# dense field additionally gives a smoother spectrum and a fuller tail.

indices = [(95, 95), (95, 5), (171, 62), (138, 161), (55, 161), (22, 62)]
px = np.array([X[i, j] for i, j in indices])
py = np.array([Y[i, j] for i, j in indices])
eta = np.array([field[i, j, :] for i, j in indices]).T
time = np.arange(NFRAMES) / SAMPLING_RATE

arr = ewdm.Arrays.from_numpy(
    time=time, surface_elevation=eta, position_x=px, position_y=py,
    fs=SAMPLING_RATE
)
pentagon = arr.compute(cross_wavelet=True, solver="least-squares", kappa=36)

fig, ax = plt.subplots()
ax.loglog(output["frequency"], output["frequency_spectrum"],
          label="multiaperture (dense field)")
ax.loglog(pentagon["frequency"], pentagon["frequency_spectrum"], "--",
          label="Arrays (pentagon)")
ax.set(xlabel="f [Hz]", ylabel="$S(f)$ [m$^2$/Hz]", title="Frequency spectrum")
ax.legend()


# %%
# The wavenumber and inverse-phase-speed spectra
# ----------------------------------------------
#
# The single fixed array can also be
# asked for a wavenumber spectrum (`coordinate="wavenumber"`), but it only
# resolves the narrow band around its own baseline and falls off a cliff at
# higher k. The multi-aperture estimator stitches the nested apertures together
# and carries a clean F(k) tail across the short waves the fixed array is blind
# to. The inverse-phase-speed spectrum Q(nu) of Björkqvist et al. (2019) tells
# the same story.

pentagon_k = arr.compute(cross_wavelet=True, solver="least-squares", kappa=36,
                         coordinate="wavenumber")
pentagon_nu = arr.compute(cross_wavelet=True, solver="least-squares", kappa=36,
                          coordinate="nu")

fig, (axk, axn) = plt.subplots(1, 2, figsize=(10, 4))
axk.loglog(output["wavenumber"], output["wavenumber_spectrum"],
           label="multiaperture (dense field)")
axk.loglog(pentagon_k["wavenumber"], pentagon_k["wavenumber_spectrum"], "--",
           label="Arrays (pentagon)")
axk.set(xlabel="k [rad/m]", ylabel="$F(k)$ [m$^3$]",
        title="Wavenumber spectrum")
top = np.nanmax(output["wavenumber_spectrum"].values)
axk.set_ylim(top * 1e-4, top * 3)
axk.legend()

axn.loglog(output["nu"], output["nu_spectrum"],
           label="multiaperture (dense field)")
axn.loglog(pentagon_nu["nu"], pentagon_nu["nu_spectrum"], "--",
           label="Arrays (pentagon)")
axn.set(xlabel=r"$\nu = k / \omega$ [s/m]", ylabel="$Q(\\nu)$ [m$^3$/s]",
        title="Inverse-phase-speed spectrum")
top = np.nanmax(output["nu_spectrum"].values)
axn.set_ylim(top * 1e-4, top * 3)
axn.legend()
fig.tight_layout()


# %%
# Finally, the polar wavenumber-direction spectrum. The dominant wave system
# appears as a single lobe at its wavenumber and direction (degrees clockwise
# from North).

fig, ax = plt.subplots(figsize=(6, 5))
plot_directional_spectrum(
    output["directional_spectrum_k"], frqs="wavenumber", ax=ax, levels=None,
    colorbar=True,
    axes_kw={"rmin": 0.2, "rmax": 2.0, "rstep": 0.4},
    cbar_kw={"label": "$\\Psi(k,\\theta)$"}
)

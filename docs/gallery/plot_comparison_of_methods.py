#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

"""
Comparison of different methods
===============================

"""

import numpy as np
import xarray as xr
import scipy.signal as signal

from matplotlib import pyplot as plt

import ewdm
from ewdm.sources import SpotterBuoysDataSource
from ewdm.plots import plot_directional_spectrum

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

# %% 
# Load sample dataset
# ---------------------------------------
#
spotter = SpotterBuoysDataSource("../../data/displacement.csv")
dataset = (
    spotter
    .read_dataset()
    .sel(time=slice("2020-08-20 10:00", "2020-08-20 10:30"))
)
print(dataset)

# plot time series of the first 5 minutes
fig, (ax1, ax2, ax3) = plt.subplots(
    3,1, figsize=(6,6), sharex=True, sharey=True
)
subset = dataset.sel(time=slice("2020-08-20T10:00", "2020-08-20T10:05"))
subset["surface_elevation"].plot(ax=ax1)
subset["eastward_displacement"].plot(ax=ax2)
subset["northward_displacement"].plot(ax=ax3)
_ = ax1.set(xlabel="", ylabel="$\\eta(t)$ [m]")
_ = ax2.set(xlabel="", ylabel="$x(t)$ [m]")
_ = ax3.set(xlabel="", ylabel="$y(t)$ [m]")


# %% 
# Classic directional spectral analysis
# -------------------------------------
#
class ClassicSpectralAnalysis(object):
    """This class implements the classic directional spectral analysis"""

    def __init__(self, dataset: xr.Dataset, fs: float=1.0, nperseg=256):
        self.dataset = dataset # time series
        self.fs = fs # sampling frequency
        self.nperseg = nperseg # number of fourier components

    def cross_spectral_matrix(self) -> xr.Dataset:

        # extract variables from dataset
        x, y, z = (
            self.dataset["eastward_displacement"],
            self.dataset["northward_displacement"],
            self.dataset["surface_elevation"]
        )

        # constants
        fft_args = {
            "fs": self.fs,
            "detrend": "constant",
            "nperseg": min(self.nperseg, len(self.dataset["time"])),
        }

        # auto-spectra
        frq, Sxx = signal.welch(x, **fft_args)
        frq, Syy = signal.welch(y, **fft_args)
        frq, Szz = signal.welch(z, **fft_args)

        # cross-spectra
        frq, Sxz = signal.csd(x, z, **fft_args)
        frq, Syz = signal.csd(y, z, **fft_args)
        frq, Sxy = signal.csd(x, y, **fft_args)

        return  xr.Dataset(
            coords = {
                "frequency": ("frequency", frq),
            },
            data_vars = {
                "Sxx": ("frequency", Sxx),
                "Syy": ("frequency", Syy),
                "Szz": ("frequency", Szz),
                "Sxz": ("frequency", Sxz),
                "Syz": ("frequency", Syz),
                "Sxy": ("frequency", Sxy),
            }
        )

    def directional_moments(self, Phi: xr.Dataset) -> xr.Dataset:
        Exx, Eyy, Ezz, Cxy, Qxz, Qyz = (
            np.real(Phi["Sxx"]), np.real(Phi["Syy"]), np.real(Phi["Szz"]),
            np.real(Phi["Sxy"]), np.imag(Phi["Sxz"]), np.imag(Phi["Syz"])
        )

        return  xr.Dataset(
            {
                "a0": Phi["Szz"],
                "a1": Qxz / np.sqrt(Ezz * (Exx + Eyy)),
                "b1": Qyz / np.sqrt(Ezz * (Exx + Eyy)),
                "a2": (Exx - Eyy) / (Exx + Eyy),
                "b2": 2 * Cxy / (Exx + Eyy)
            }
        )


csp = ClassicSpectralAnalysis(dataset, fs=dataset.sampling_rate)

# computing the cross-spectral matrix
Phi = csp.cross_spectral_matrix()
print(Phi)

# computing the directional moments, aka first-five Fourier coefficients
moments = csp.directional_moments(Phi)
print(moments)


# %% 
# Truncated Fourier Series
# ------------------------
def tfs_distribution(moments, smoothing=32):
    dirs =  xr.Variable(dims=("direction"), data=np.arange(-180,180,5))
    D = (
        1/2 +
        moments["a1"] * np.cos(np.radians(dirs)) +
        moments["b1"] * np.sin(np.radians(dirs)) +
        moments["a2"] * np.cos(2*np.radians(dirs)) +
        moments["b2"] * np.sin(2*np.radians(dirs))
    )
    D.coords["direction"] = dirs
    D = D.where(D > 0, 0)
    return (
        (D / D.integrate("direction"))
        .rolling(frequency=smoothing, center=True)
        .median()
    )

D_tfs = tfs_distribution(moments, smoothing=2)
print(D_tfs)


# %% 
# Maximum Entropy Method
# ----------------------
def mem_distribution(moments, smoothing=32):
    dirs =  xr.Variable(dims=("direction"), data=np.arange(-180,180,5))

    c1 = moments["a1"] + 1j*moments["b1"]
    c2 = moments["a2"] + 1j*moments["b2"]

    phi1 = (c1 - c2 * c1.conj()) / (1 - np.abs(c1)**2)
    phi2 = c2 - c1.conj() * phi1

    sigma_e = 1 - phi1 * c1.conj() - phi2 * c2.conj()

    D = (1/(2*np.pi)) * np.real(
        sigma_e.expand_dims({"direction": dirs}) /
        np.abs(
            1 - phi1.expand_dims({"direction": dirs}) * np.exp(-1j*dirs*np.pi/180)
              - phi2.expand_dims({"direction": dirs}) * np.exp(-2j*dirs*np.pi/180)
        )**2
    )

    return (
        (D.T / D.integrate("direction"))
        .rolling(frequency=smoothing, center=True)
        .median()
    )

D_mem = mem_distribution(moments, smoothing=2)
print(D_mem)


# %% 
# Wavelet-based directional wave spectrum
# ---------------------------------------
#
spec = ewdm.Triplets(dataset)
output = spec.compute()
D_ewdm = output["directional_distribution"]
print(D_ewdm)


# %% 
# Direactional distribution function
# ----------------------------------
#
vmin, vmax = 0, 0.05
axes_kw={"rmax": 0.8, "rstep": 0.2, "as_period": True}

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(8,8/3))
plot_directional_spectrum(
    D_tfs, ax=ax1, levels=None, colorbar=False,
    axes_kw=axes_kw, vmin=vmin, vmax=vmax
)
plot_directional_spectrum(
    D_mem, ax=ax2, levels=None, colorbar=False,
    axes_kw=axes_kw, vmin=vmin, vmax=vmax
)
plot_directional_spectrum(
    D_ewdm, ax=ax3, levels=None, colorbar=True, 
    cbar_kw={"label": "$D(f,\\theta)$"},
    axes_kw=axes_kw, vmin=vmin, vmax=vmax
)
_ = ax1.set(xlabel="", title="TFS")
_ = ax2.set(xlabel="", ylabel="", title="MEM")
_ = ax3.set(xlabel="", ylabel="", title="EWDM")

# %% 
# Directional wave spectrum
# -------------------------
#
E_tfs = moments["a0"] * D_tfs
E_mem = moments["a0"] * D_mem
E_ewdm = output["directional_spectrum"]


vmin, vmax = -3, 0
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(8,8/3), layout="constrained")
plot_directional_spectrum(
    E_tfs.pipe(np.log10), ax=ax1, levels=None, colorbar=False,
    axes_kw=axes_kw, vmin=vmin, vmax=vmax
)
plot_directional_spectrum(
    E_mem.pipe(np.log10), ax=ax2, levels=None, colorbar=False,
    axes_kw=axes_kw, vmin=vmin, vmax=vmax
)
plot_directional_spectrum(
    E_ewdm.pipe(np.log10), ax=ax3, levels=None, colorbar=True,
    axes_kw=axes_kw, vmin=vmin, vmax=vmax
)

# %% 
# Azimutally-integrated power spectrum
# ------------------------------------
# A quick comparison between azimutally-integrated wavelet power and Fourier-based
# power density spectrum.
#
fig, ax = plt.subplots()
ax.loglog(moments["frequency"], moments["a0"], label="Fourier-based")
output.frequency_spectrum.plot(ax=ax, ls="--", lw=2, label="Wavelet-based")
ax.legend()

#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

"""
Comparison of different methods
===============================

It is a common practice to express the directional wave spectrum,
:math:`E(f,\\theta)`, as the product of two functions:

    .. math:: E(f,\\theta) = S(f) D(f,\\theta)

where :math:`S(f)` is the frequency spectrum and :math:`D(f,\\theta)` is the
directional distribution function. The problem basically consist of
finding :math:`D(f,\\theta)`.

This section compares different methods for estimating the directional
distribution :math:`D(f,\\theta)` and the the directional wave
spectra :math:`E(f,\\theta)`. We compare the wavelet-based directional method
implemented in EWDM with two of the conventional methods based on Fourier
spectral analysis, namely Truncated Fourier Series (TSF) and Maximum Entropy
Method (MEM).

The results are visualised and compared to understand the performance and
differences between these methods.
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
# For this example, were use data from a Spotter buoy deployed off
# the west coast of Ireland in August 2020. In particular, we are using
# a sample of 30 minutes.
#
# The dataset contains the two horizontal wave-induced displacements and
# the vertical sea surface elevation.
#
spotter = SpotterBuoysDataSource("../../data/displacement.csv")
dataset = (
    spotter
    .read_dataset()
    .sel(time=slice("2020-08-20 10:00", "2020-08-20 10:30"))
)
print(dataset)

# %%
# Let's plot a subsample of 5 minutes of these time series to have
# an idea on what they look like. As it can be seen, the amplitude of
# the wave-induced buoy movement, in both the horizontal and vertical
# coordinates, is about 3-4 metres.

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
# The following class implements the conventional cross-spectral analysis
# using Fourier transform. This class takes in a dataset containing eastward
# and northward displacements along with the surface elevation data. The
# input parameters are sampling frequency of the time series and number of
# Fourier components for analysis.
#
# Considering the three time series, :math:`x(t)`, :math:`y(t)` and
# :math:`z(t)` (We normally use :math:`\eta(t)` for the vertical but it is
# changed here to :math:`z(t)` for simplicity and convenience), the
# cross-spectral matrix is computed as:
#
# .. math:: \Phi(f) = \begin{bmatrix}
#                        S_{xx} & S_{xy} & S_{xz} \\
#                        S_{yx} & S_{yy} & S_{yz} \\
#                        S_{zx} & S_{zy} & S_{zz}
#                     \end{bmatrix}
#
# where :math:`S_{xy}(f)` is the Fourier cross-spectrum between
# :math:`x(t)` and :math:`y(t)`, respectively.
#
# Each cross-spectrum can be written in terms of a real (co-spectrum)
# and an imaginary (quad-spectrum) component:
#
# .. math:: S_{xy}(f) = C_{xy} + i Q_{xy}
#
# The auto-spectrum is the cross-spectrum of the same signal, and can be
# written as :math:`E_{xx}(f) = S_{xx} S_{xx}^*`, where :math:`*`
# denotes a complex conjugate.
#
# For typical buoy recording of wave-induced displacements, the circular
# moments can be written in terms of these auto-, co-, and quad-spectra, like:
#
# .. math:: a_0 = S_{zz}(f)
# .. math:: a_1 = \frac{Q_{xz}}{\sqrt{E_{zz} (E_{xx} + E_{yy})}}
# .. math:: a_2 = \frac{Q_{yz}}{\sqrt{E_{zz} (E_{xx} + E_{yy})}}
# .. math:: b_1 = \frac{E_{xx} - E_{yy}}{E_{xx} + E_{yy}}
# .. math:: b_2 = \frac{2 C_{xy}}{E_{xx} + E_{yy}}
#
# See `Kuik et al. (1988)`_ and Appendix C in `Peláez-Zapata et al. (2024)`_ for further details.
#
# .. _Kuik et al. (1988): https://doi.org/10.1175/1520-0485(1988)018<1020:AMFTRA>2.0.CO;2
# .. _Peláez-Zapata et al. (2024): https://theses.fr/2024UPASM004

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
                "a0": Ezz,
                "a1": Qxz / np.sqrt(Ezz * (Exx + Eyy)),
                "b1": Qyz / np.sqrt(Ezz * (Exx + Eyy)),
                "a2": (Exx - Eyy) / (Exx + Eyy),
                "b2": 2 * Cxy / (Exx + Eyy)
            }
        )


csp = ClassicSpectralAnalysis(dataset, fs=dataset.sampling_rate)

# computing the cross-spectral matrix
Phi = csp.cross_spectral_matrix()
print("\nThe cross-spectra matrix is:\n" + 28*"-")
print(Phi)

# computing the directional moments, aka first-five Fourier coefficients
moments = csp.directional_moments(Phi)
print("\nThe circular moments are:\n" + 25*"-")
print(moments)


# %%
# Truncated Fourier Series
# ------------------------
#
# The most straightforward way of obtaining the directional distribution
# function :math:`D(f,\theta)` is by expanding in Fourier series:
#
# .. math:: D(f,\theta) = \frac{1}{2\pi} \left[ 1 + \sum_{n=1}^{N}
#                a_n(f) \cos n\theta + b_n(f) \sin n\theta
#           \right]
#
# This Fourier series is truncated to :math:`N=2` because it is the maximum
# number of components that can be extracted from buoy measurements.
#
def tfs_distribution(moments, smoothing=32):
    """Implementation of the Truncated Fourier Series method"""

    dirs =  xr.Variable(dims=("direction"), data=np.arange(-180,180,5))
    D = (1/(2*np.pi)) * (
        1 +
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
#
# A more sophisticated way of obtaining the directional distribution
# function :math:`D(f,\theta)` is by using a maximum entropy estimator.
#
# Following `Lygre and Krogstad (1983)`_ and `Alves and Melo (1999)`_,
# the form of the directional distribution function can be written as:
#
# .. math:: D(f,\theta) = \frac{1}{2\pi} \left[
#                   \frac{1 - \phi_1 c_1^* - \phi_2 c_2^*}
#                        { |1 - \phi_1 e^{-i\theta} - \phi_2 e^{-i2\theta}|^2 }
#               \right]
#
# where :math:`c_1` and :math:`c_2` are the complex representation of the
# Fourier coefficients, i.e.,
#
# .. math:: c_1(f) = a_1(f) + i b_1(f)
# .. math:: c_2(f) = a_2(f) + i b_2(f)
#
# and
#
# .. math:: \phi_1 = \frac{c_1 - c_2 c_1^*}{1 - |c_1|^2}
# .. math:: \phi_2  = c_2 - c_1^* \phi_1
#
# It is worth noting that this is just one of the possible implementations
# of MEM. There are other variations that might potentially produce better
# results. For more details, see `Christie (2024)`_ and
# `Simanesew et al. (2018)`_.
#
# .. _Lygre and Krogstad (1983): https://doi.org/10.1175/1520-0485(1986)016<2052:MEEOTD>2.0.CO;2
# .. _Alves and Melo (1999): https://doi.org/10.1016/S0141-1187(99)00019-X
# .. _Christie (2024): https://www.sciencedirect.com/science/article/pii/S0141118723003711?via%3Dihub
# .. _Simanesew et al. (2018): https://doi.org/10.1175/JTECH-D-17-0007.1
#
def mem_distribution(moments, smoothing=32):
    """Implementation of the Maximum Entropy Method"""

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
# Here, we use the `ewdm.Triplets` module to get an estimation of the
# directional distribution function :math:`D(f,\theta)` using the
# wavelet-based method applied over the buoy data.
#
spec = ewdm.Triplets(dataset)
output = spec.compute()
D_ewdm = output["directional_distribution"]
print(D_ewdm)


# %%
# Directional distribution function
# ----------------------------------
#
# The results for the directional distribution function :math:`D(f,\theta)`
# obtained from the three different methods evaluated are shown side-by-side.
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
# Now, we show the corresponding directional wave spectra
# :math:`E(f,\theta)` for each method. Recalling that this quantity is
# constructed as:
#
# .. math:: E(f,\theta) = S(f) D(f,\theta)
#
# For the Fourier-based methods, the frequency spectrum :math:`S(f)` is
# simply the Fourier transform of the surface elevation signal, i.e.,
# :math:`E_{zz}(f)`, the auto-spectrum of :math:`z(t)`. For the wavelet
# method, :math:`S(f)` is obtained time-averaging the squared
# wavelet amplitudes. See `Directional distribution and directional spectrum <https://extended-wdm.readthedocs.io/en/latest/maths.html#directional-distribution-and-directional-spectrum>`_.
#

# compute the directional wave spectrum
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
# A quick comparison between azimutally-integrated wavelet power and
# Fourier-based power density spectrum, :math:`S(f)`, is presented.
#
# This comparison reveals similarities. The wavelet spectrum inherently
# offers lower frequency resolution compared to the Fourier spectrum,
# thus resulting in a naturally smoother representation of the energy
# distribution. This characteristic of the continuous wavelet transform
# has been extensively discussed. See e.g., `Torrence and Compo (1998)`_.
#
# .. _Torrence and Compo (1998): http://journals.ametsoc.org/doi/10.1175/1520-0477(1998)079<0061:APGTWA>2.0.CO;2
#
fig, ax = plt.subplots()
ax.loglog(moments["frequency"], moments["a0"], label="Fourier-based")
output.frequency_spectrum.plot(ax=ax, ls="--", lw=2, label="Wavelet-based")
ax.set(xlabel="$f$ [Hz]", ylabel="$S(f)$ [m$^2$/Hz]")
ax.legend()

# %%
# Discussion and conclusion
# -------------------------
#
# The comparison of three methods for estimating the directional wave
# spectrum reveals distinct characteristics:
#
# - TFS is too broad, i.e., it  spreads too much the energy across the
#   directions, which is inconsistent with previous observations.
#
# - MEM seems to be too narrow and shows spurious peaks. In addition, it
#   tends to produce inconsistent bimodal distributions. For example, it
#   identifies waves travelling south which appear to e unrealistic.
#
# - EWDM results appear consistent showing a smooth transition of wave
#   energy across frequencies and directions.
#
# However, objectively evaluating the quality of these estimations is
# challenging due to the lack of a ground truth reference. Consequently,
# further research is imperative to enhance the understanding and
# assessment of these methods. Despite this limitation, EWDM exhibits promise
# in providing more reliable estimates for the directional wave spectrum.

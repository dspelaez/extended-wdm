#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

"""
Multi-aperture directional wave spectra
=======================================

A single fixed array of wave staffs resolves only a narrow band of wavenumbers.
The phase difference between two staffs separated by a baseline :math:`b` wraps
once :math:`k\\,b > \\pi`, so an array is trustworthy only below
:math:`k_{max} = \\pi / b_{max}`. A 2 m pentagon therefore stops at
:math:`k_{max} \\approx 1.6` rad/m, far short of the wind sea.

`ewdm.MultiApertureArrays` removes this ceiling. Given a dense elevation field
:math:`\\eta(y, x, t)` it seeds virtual staffs into a *ladder of nested
apertures*: coarse, widely spaced clusters for the long waves and tight clusters
for the short waves. Each aperture is trusted only within its own anti-alias band
and the per-aperture spectra are stitched into one estimate that spans the swell
through the wind sea, all measured directly from the spatial phase gradients
without assuming a dispersion relation.

This example runs the method end to end on the `Guimaraes et al (2020)`_ Black
Sea stereo-video field (19 m field of view, 0.1 m pixels). Supporting functions
live in ``multiaperture_helper.py`` next to this script.

.. _Guimaraes et al (2020): https://doi.org/10.1038/s41597-020-0492-9
"""

# sphinx_gallery_thumbnail_number = 2

import os
import sys
import subprocess

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle

import ewdm
from ewdm import MultiApertureArrays
from ewdm.multiaperture import auto_apertures, k_dispersion, seed_aperture
from ewdm.parameters import GRAV

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))
                if "__file__" in globals() else os.getcwd())
from multiaperture_helper import (
    circular_tukey, pentagon_spectra, direct_spectra, spreading)

# %%
# Load the field and run the shared computation
# ---------------------------------------------
#
# The stereo-video file is downloaded once into the data folder. We then build
# the dense field, run the multi-aperture estimator (its default ladder is sized
# to the field of view), and pre-compute the reference spectra that the figures
# below compare against: two fixed pentagons (`ewdm.Arrays`) and the direct 3-D
# FFT spectra.

URL = ("https://data-dataref.ifremer.fr/stereo/BS_2013/"
       "2013-09-22_13-00-01_10Hz/nc/Surfaces_20130922_130001_short.nc")
LOCAL = os.path.join("../../data", "stereo-video.nc")
if not os.path.exists(LOCAL):
    print("Downloading the stereo-video file (a few minutes the first time).")
    subprocess.call(f"wget {URL} -O {LOCAL}", shell=True)

o = nc.Dataset(LOCAL)
X = np.asarray(o["X"][:], float)
Y = np.asarray(o["Y"][:], float)
dx = (np.nanmax(X) - np.nanmin(X)) / (X.shape[1] - 1)
fs, depth, NT = 10.0, 30.0, 2048
field = np.transpose(np.asarray(o["Z"][:NT], float), (1, 2, 0))   # (ny, nx, time)
ny, nx, T = field.shape
t = np.arange(NT) / fs
var_eta = float((field - field.mean(2, keepdims=True)).var(2).mean())
ci = cj = (nx - 1) // 2

freqs = np.logspace(np.log10(0.05), np.log10(2.0), 60)
k_grid = 2.0 ** np.linspace(np.log2(0.03), np.log2(np.pi / dx), 90)
nu_grid = 2.0 ** np.linspace(np.log2(0.02), np.log2(3.0), 80)

# multi-aperture estimate (data-aware default ladder, per-aperture diagnostics on)
M = MultiApertureArrays.from_field(
    field, dx, depth=depth, fs=fs, direction_convention="math").compute(
    freqs=freqs, k_grid=k_grid, nu_grid=nu_grid, n_staff=16, seed=20,
    nu_f_lim=None, nu_k_lim=None, lo_frac_broad=1.0, rel_bandwidth=0.10,
    return_apertures=True)

# fixed pentagons: a broad one and the 2 m ASIS pentagram
R = int(round(1.0 / dx))
broad = pentagon_spectra(field, X, Y, dx, fs, t,
                         [(95, 95), (95, 5), (171, 62), (138, 161), (55, 161), (22, 62)])
narrow_idx = [(ci, cj)] + [(int(round(ci - R * np.sin(np.radians(90 + 72 * i)))),
                            int(round(cj + R * np.cos(np.radians(90 + 72 * i)))))
                           for i in range(5)]
narrow = pentagon_spectra(field, X, Y, dx, fs, t, narrow_idx)

D = direct_spectra(field, dx, fs, var_eta, freqs, k_grid)          # direct 3-D FFT

# shared style / handy aliases
MULT, EQC, SATC, DIRC, DISPC, FIXC, ZTOP = "#e31a1c", "#7e3ff2", "C1", "k", "C2", "C0", 12
SF_M = M["frequency_spectrum"].values
kw = M["wavenumber"].values
fpk = freqs[np.nanargmax(SF_M)]
kpk = kw[np.nanargmax(M["wavenumber_spectrum"].values)]
topS, topK = np.nanmax(SF_M), np.nanmax(M["wavenumber_spectrum"].values)
print(f"Hs = {4*np.sqrt(var_eta):.2f} m   fp = {fpk:.2f} Hz   "
      f"apertures = {len(M['aperture_name'])}")

# %%
# The problem: a single fixed array
# ---------------------------------
#
# The 2 m pentagon follows the truth up to its baseline limit, then throws up a
# spurious hump and a sharp energy cliff at the anti-alias ceiling
# :math:`k_{max} = \pi / b_{max}`, blind to everything shorter.

k_n = narrow["Fk"]["wavenumber"].values
F_n = narrow["Fk"]["wavenumber_spectrum"].values * narrow["fac"]
fig, ax = plt.subplots(figsize=(7, 5))
ax.loglog(D["klin"], D["klin"] ** 3 * D["Fk_d"], "--", color="0.45", lw=1.8,
          label="true spectrum (3-D FFT)")
ax.loglog(k_n, k_n ** 3 * F_n, "-", color="C0", lw=2.6,
          label="single fixed array (2 m)")
ax.axvline(narrow["kmax"], color="red", ls=":", lw=1.6)
ax.text(narrow["kmax"] * 0.93, 2e-4, r"$k_{max}=\pi/b_{max}$", rotation=90,
        va="top", ha="right", color="red", fontsize="small")
ax.set(xscale="log", yscale="log", xlim=(1e-2, 3e1), ylim=(1e-5, 1e-1),
       xlabel=r"$k$ [rad m$^{-1}$]", ylabel=r"$k^{3}\,F(k)$ [rad]")
ax.grid(which="major", ls="-", lw=0.5, alpha=0.6)
ax.grid(which="minor", ls=":", lw=0.4, alpha=0.5)
ax.legend(loc="upper left")

# %%
# The aperture ladder
# -------------------
#
# The default ladder is sized to the field: coarse windows for the swell down to
# tight five-staff *plus* crosses for the wind sea (no hand-tuning). Four of the
# larger apertures are drawn over a field snapshot; the panel below plots every
# aperture's saturation spectrum :math:`k^{3} F(k)` (the four highlighted, the
# rest dotted) with the stitched composite and the direct 3-D FFT.

apname = [str(a) for a in M["aperture_name"].values]
apFk, klo, khi = M["aperture_Fk"].values, M["aperture_klo"].values, M["aperture_khi"].values
kg, Fk_M, extd = M["wavenumber"].values, M["wavenumber_spectrum"].values, dict(auto_apertures(ny, nx, dx))
bmax_all = np.pi / np.asarray(khi); bl = np.log(bmax_all)
cands = [i for i in range(len(bl)) if bmax_all[i] >= 2.0] or \
        sorted(range(len(bl)), key=lambda i: -bmax_all[i])[:4]
show_idx = []
for tgt in np.linspace(bl[cands].max(), bl[cands].min(), 4):
    pool = [i for i in cands if i not in show_idx]
    show_idx.append(pool[int(np.argmin(np.abs(bl[pool] - tgt)))])
palette = ['#5e3c99', '#1b9e77', '#e6ab02', '#1f78b4']
colof = {idx: palette[n] for n, idx in enumerate(show_idx)}

i_snap = 1500 if NT > 1500 else NT // 2
Zsnap = (field - field.mean(2, keepdims=True))[:, :, i_snap]
xg = (np.arange(nx) - (nx - 1) / 2.0) * dx
yg = (np.arange(ny) - (ny - 1) / 2.0) * dx
vZ = np.percentile(np.abs(Zsnap), 99)

fig = plt.figure(figsize=(8.5, 13), constrained_layout=True)
gs = fig.add_gridspec(3, 4, height_ratios=[1, 1, 1.7])
fax = []
for (r, c), idx in zip([(0, 0), (0, 2), (1, 0), (1, 2)], show_idx):
    nm, col, ext, w = apname[idx], colof[idx], extd[apname[idx]], float(bmax_all[idx])
    ax = fig.add_subplot(gs[r, c:c + 2]); fax.append(ax)
    im = ax.pcolormesh(xg, yg, Zsnap, cmap="coolwarm", vmin=-vZ, vmax=vZ,
                       shading="auto", rasterized=True); ax.set_aspect("equal")
    ax.add_patch(Circle((0, 0), w / 2, fill=False, edgecolor=col, lw=2.5))
    ii, jj, _, _, _ = seed_aperture(ny, nx, dx, ext, 16, 20 + idx)
    gx, gy = (jj - (nx - 1) / 2.0) * dx, (ii - (ny - 1) / 2.0) * dx
    ins = np.hypot(gx, gy) <= w / 2.0 + 1e-9
    ax.scatter(gx[ins], gy[ins], s=14, facecolor="white", edgecolor="black", lw=0.5, zorder=5)
    ax.text(0.94, 0.94, nm, transform=ax.transAxes, color=col, fontsize="small",
            fontweight="bold", ha="right", va="top", zorder=6,
            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.7))
    ax.set(xlim=(-10, 10), ylim=(-10, 10), xticks=[-10, -5, 0, 5, 10], yticks=[-10, -5, 0, 5, 10])
    ax.set_ylabel("y [m]") if c == 0 else ax.set_yticklabels([])
    ax.set_xlabel("x [m]") if r == 1 else ax.set_xticklabels([])
fig.colorbar(im, ax=fax, location="top", orientation="horizontal",
             label=r"$\eta$ [m]", shrink=0.6, aspect=40, pad=0.02)

ax = fig.add_subplot(gs[2, :]); k3 = kg ** 3; sel = set(show_idx); _lbl = True
for idx in range(len(apname)):
    if idx in sel:
        continue
    ax.loglog(kg, k3 * apFk[idx], ":", color="k", lw=0.7, alpha=0.45,
              label=("other apertures" if _lbl else None)); _lbl = False
for idx in show_idx:
    inb = (kg >= klo[idx]) & (kg <= khi[idx])
    ax.loglog(kg, k3 * apFk[idx], ":", color=colof[idx], lw=1.2, alpha=0.7)
    ax.loglog(kg[inb], (k3 * apFk[idx])[inb], "-", color=colof[idx], lw=2.0, label=apname[idx])
fin = np.isfinite(Fk_M)
ax.loglog(kg[fin], (k3 * Fk_M)[fin], "-", color=MULT, lw=2.6, label="multi-aperture (composite)")
ax.loglog(D["klin"], D["klin"] ** 3 * D["Fk_d"], "k--", lw=2.4, label="direct")
ax.set(xlim=(1e-2, 5e1), ylim=(1e-5, 1e-1),
       xlabel=r"k [rad m$^{-1}$]", ylabel=r"k$^3$ F(k) [rad]")
ax.grid(which="major", ls="-", lw=0.6); ax.grid(which="minor", ls=":", lw=0.5)
ax.legend(loc="upper left", ncol=2, fontsize="small")

# %%
# Omnidirectional spectra
# -----------------------
#
# The composite carries a clean :math:`F(k)` tail far past the reach of any
# single pentagon, out to :math:`k \sim 10` rad/m, in agreement with the direct
# FFT and the dispersion image of the frequency spectrum.

def sl(x, xr, yr, p):
    xx = x[x >= xr]
    return xx, yr * (xx / xr) ** p

def yat(x, y, xr):
    return float(np.interp(xr, x, y))

def f_of_k(k):
    return np.sqrt(GRAV * k * np.tanh(np.clip(k * depth, 1e-9, 50))) / (2 * np.pi)

floorS, floorF = lambda s: 2 * s ** 2 / fs, lambda s: s ** 2 * dx / np.pi
SIG_REP, SIG_FAR = 0.01, 0.03
kk_d = k_dispersion(freqs, depth)
Fk_disp = SF_M * np.abs(np.gradient(freqs, kk_d))

fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4.6), gridspec_kw=dict(wspace=0.24))
a.loglog(freqs, D["Sf_d"], color=DIRC, lw=1.2, alpha=.75, label="direct 3-D FFT")
xr_ = 2.2 * fpk
a.loglog(*sl(freqs, xr_, yat(freqs, SF_M, xr_), -4), ":", color=EQC, lw=2,
         label=r"equil. ($f^{-4}$, $k^{-5/2}$)")
a.loglog(*sl(freqs, xr_, 0.45 * yat(freqs, SF_M, xr_), -5), "-.", color=SATC, lw=2,
         label=r"sat. ($f^{-5}$, $k^{-3}$)")
a.loglog(freqs, SF_M, lw=2.2, color=MULT, zorder=ZTOP, label="multiaperture (field)")
a.axvline(f_of_k(broad["kmax"]), color="0.5", ls=":", lw=1)
a.axhspan(floorS(SIG_REP), floorS(SIG_FAR), color="0.6", alpha=0.18, zorder=0)
a.axhline(floorS(SIG_REP), color="0.45", lw=1.3); a.axhline(floorS(SIG_FAR), color="0.45", lw=0.7)
a.set(xlim=(0.05, 5), ylim=(min(topS * 3e-4, floorS(SIG_REP) / 3), topS * 4),
      xlabel="f [Hz]", ylabel=r"$S(f)$ [m$^2$/Hz]")
a.text(0.052, np.sqrt(floorS(SIG_REP) * floorS(SIG_FAR)), "WASS quant. floor (~1$\\to$3 cm)",
       fontsize="small", color="0.3", va="center", ha="left")
a.legend(loc="upper right", fontsize="small")

b.loglog(broad["Fk"]["wavenumber"], broad["Fk"]["wavenumber_spectrum"] * broad["fac"],
         "--", color=FIXC, label="broad pentagram")
b.loglog(narrow["Fk"]["wavenumber"], narrow["Fk"]["wavenumber_spectrum"] * narrow["fac"],
         ":", color="tab:blue", lw=1.8, label="narrow pentagram (2 m)")
b.loglog(D["klin"], D["Fk_d"], "o-", color=DIRC, ms=3, lw=1, alpha=.7, label="_nolegend_")
b.loglog(kk_d, Fk_disp, "--", color=DISPC, lw=1.5, label=r"disp. image of $S(f)$")
kr = 2 * kpk
b.loglog(*sl(kw, kr, yat(kw, M["wavenumber_spectrum"].values, kr), -2.5), ":", color=EQC, lw=2, label="_nolegend_")
b.loglog(*sl(kw, kr, 0.6 * yat(kw, M["wavenumber_spectrum"].values, kr), -3), "-.", color=SATC, lw=2, label="_nolegend_")
b.loglog(kw, M["wavenumber_spectrum"], lw=2.2, color=MULT, zorder=ZTOP, label="_nolegend_")
for FX in (broad, narrow):
    b.axvline(FX["kmax"], color="0.5", ls=":", lw=1.3)
    b.text(FX["kmax"], 1.02, r"$k_{max}$=%.2f" % FX["kmax"], transform=b.get_xaxis_transform(),
           color="0.3", va="bottom", ha="center", fontsize="small", clip_on=False)
b.axhspan(floorF(SIG_REP), floorF(SIG_FAR), color="0.6", alpha=0.18, zorder=0)
b.axhline(floorF(SIG_REP), color="0.45", lw=1.3); b.axhline(floorF(SIG_FAR), color="0.45", lw=0.7)
b.set(xlim=(0.03, np.pi / dx), ylim=(min(topK * 3e-4, floorF(SIG_REP) / 3), topK * 4),
      xlabel="k [rad/m]", ylabel=r"$F(k)$ [m$^3$]")
b.legend(loc="upper right", fontsize="small")

# %%
# Directional spectra and spreading
# ---------------------------------
#
# Because each aperture resolves direction as well as wavenumber, the stitched
# estimate is a full directional spectrum. We compare it against the
# wavelet-direction (WDM) estimate of the fixed 2 m pentagon and the direct FFT,
# as the directional spectrum :math:`E(X,\theta)` and, row-normalised, the
# directional spreading :math:`D(X,\theta)`.

cols = ["fixed pentagram (2 m)", "multiaperture", "direct 3-D FFT"]
daf = {"fixed pentagram (2 m)": narrow["Sf"]["directional_spectrum"],
       "multiaperture": M["directional_spectrum_f"], "direct 3-D FFT": D["Dfft"]}
dak = {"fixed pentagram (2 m)": narrow["Fk"]["directional_spectrum"],
       "multiaperture": M["directional_spectrum_k"], "direct 3-D FFT": D["Dfkd"]}
FLIM, KLIM = (0.1, 1.5), (0.05, 4.0)
L_fov = nx * dx
f_fov = np.sqrt(GRAV * (2 * np.pi / L_fov)) / (2 * np.pi)
k_fov = 2 * np.pi / L_fov
SFd = {c: spreading(daf[c], "frequency", *FLIM) for c in cols}
th, fv, Dd = SFd["direct 3-D FFT"]; Dd = Dd.copy(); Dd[fv < f_fov] = np.nan
SFd["direct 3-D FFT"] = (th, fv, Dd)
SKd = {c: spreading(dak[c], "wavenumber", *KLIM,
                    blank_lowE=(c not in ("direct 3-D FFT", "multiaperture"))) for c in cols}
th, kv, Dd = SKd["direct 3-D FFT"]; Dd = Dd.copy(); Dd[kv < k_fov] = np.nan
SKd["direct 3-D FFT"] = (th, kv, Dd)


def two_wide(panels_f, panels_k, cbl_f, cbl_k, cmap, kw_f, kw_k, hline=None):
    nrow = len(cols)
    fig, axS = plt.subplots(nrow, 2, figsize=(9.2, 3.1 * nrow), gridspec_kw=dict(wspace=0.34))
    for i, c in enumerate(cols):
        thh, xv, Z = panels_f[c]; axS[i, 0].pcolormesh(thh, xv, Z, cmap=cmap, shading="auto", **kw_f)
        axS[i, 0].set(ylabel=r"$f$ [Hz]", yscale="log", ylim=FLIM, xlim=(-180, 180))
        thh, xv, Z = panels_k[c]; axS[i, 1].pcolormesh(thh, xv, Z, cmap=cmap, shading="auto", **kw_k)
        axS[i, 1].set(ylabel=r"$k$ [rad/m]", yscale="log", ylim=KLIM, xlim=(-180, 180))
        # field-of-view lower limit of the direct 3-D FFT (one wavelength across
        # the field of view)
        if hline is not None:
            axS[i, 0].axhline(f_fov, ls="--", color=hline, lw=1.3, zorder=4)
            axS[i, 1].axhline(k_fov, ls="--", color=hline, lw=1.3, zorder=4)
        for ax in (axS[i, 0], axS[i, 1]):
            ax.set_xticks(np.arange(-180, 181, 45))
        axS[i, 1].text(1.14, 0.5, c, transform=axS[i, 1].transAxes, rotation=-90,
                       va="center", ha="left", fontsize="medium")
        if i < nrow - 1:
            axS[i, 0].set_xticklabels([]); axS[i, 1].set_xticklabels([])
        else:
            axS[i, 0].set_xlabel(r"direction $\theta$ [deg]")
            axS[i, 1].set_xlabel(r"direction $\theta$ [deg]")
    for col, lab in ((0, cbl_f), (1, cbl_k)):
        cb = fig.colorbar(axS[0, col].collections[0], ax=axS[:, col].tolist(),
                          location="top", orientation="horizontal", shrink=0.9, aspect=30, pad=0.03)
        cb.set_label(lab)
    return fig

# directional spectrum E(X,theta) = D(X,theta) * S_omni(X) on the multiaperture backbone
mf_x, mf_y = M["frequency"].values, M["frequency_spectrum"].values
mk_x, mk_y = kw, M["wavenumber_spectrum"].values
EF = {c: (SFd[c][0], SFd[c][1], SFd[c][2] * np.interp(SFd[c][1], mf_x, mf_y)[:, None]) for c in cols}
EK = {c: (SKd[c][0], SKd[c][1], SKd[c][2] * np.interp(SKd[c][1], mk_x, mk_y)[:, None]) for c in cols}
vmaxEf = np.nanpercentile(np.concatenate([EF[c][2].ravel() for c in cols]), 99.5)
vmaxEk = np.nanpercentile(np.concatenate([EK[c][2].ravel() for c in cols]), 99.5)
_vir = plt.cm.viridis(np.linspace(0, 1, 256)); _vir[:, 3] = np.clip(np.linspace(0, 1, 256) / 0.06, 0, 1)
epss = mcolors.ListedColormap(_vir); epss.set_bad("0.9")
two_wide(EF, EK, r"$E(f,\theta)$  [m$^2$ Hz$^{-1}$ rad$^{-1}$]",
         r"$\Psi(k,\theta)$  [m$^3$ rad$^{-1}$]", epss,
         dict(norm=LogNorm(vmin=vmaxEf * 10 ** -2.5, vmax=vmaxEf)),
         dict(norm=LogNorm(vmin=vmaxEk * 10 ** -2.5, vmax=vmaxEk)), hline="red")

# %%
# Directional spreading (row-normalised to unit directional integral):

vmaxDf = np.nanpercentile(np.concatenate([SFd[c][2].ravel() for c in cols]), 99)
vmaxDk = np.nanpercentile(np.concatenate([SKd[c][2].ravel() for c in cols]), 99)
mag = plt.cm.magma.copy(); mag.set_bad("0.92")
two_wide(SFd, SKd, r"$D(f,\theta)$  [rad$^{-1}$]", r"$D(k,\theta)$  [rad$^{-1}$]", mag,
         dict(vmin=0, vmax=vmaxDf), dict(vmin=0, vmax=vmaxDk), hline="white")

#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

import xarray as xr
import numpy as np
from scipy import signal

from .helpers import get_sampling_frequency


NPERSEG = 256


# class ClassicSpectralAnalysis(object):
    # """This class implements the classic spectral directional analysis"""

    # def __init__(self, dataset: xr.Dataset, fs: float = 1.0):

        # # try to interpolate before computing fft
        # if (ds["z"].isnull().sum() / len(ds["time"])) > MAX_NAN_RATIO:
            # ds_interp = ds * np.nan
        # else:
            # ds_interp = ds.interpolate_na(
                # dim="time", method="quadratic", max_gap=MAX_GAP,
                # fill_value="extrapolate"
            # )
        # self.ds_interp = ds_interp
        # self.fs = get_sampling_frequency(self.ds)

    # def cross_spectral_matrix(self) -> xr.Dataset:

        # # extract variables from dataset
        # x, y, z = ds_interp["x"], ds_interp["y"], ds_interp["z"]

        # # constants
        # detrend = "constant"
        # nperseg = min(NPERSEG, len(ds.time)

        # # auto-spectra
        # frq, Sxx = signal.welch(x, fs=self.fs, nperseg=nperseg, detrend=detrend)
        # frq, Syy = signal.welch(y, fs=self.fs, nperseg=nperseg, detrend=detrend)
        # frq, Szz = signal.welch(z, fs=self.fs, nperseg=nperseg, detrend=detrend)

        # # cross-spectra
        # frq, Sxz = signal.csd(x, z, fs=self.fs, nperseg=nperseg, detrend=detrend)
        # frq, Syz = signal.csd(y, z, fs=self.fs, nperseg=nperseg, detrend=detrend)
        # frq, Sxy = signal.csd(x, y, fs=self.fs, nperseg=nperseg, detrend=detrend)

        # # find relative frequency
        # fp =  np.trapz(frq * Szz**4, x=frq) / np.trapz(Szz**4, x=frq)
        # ffp = frq / fp

        # return  xr.Dataset(
            # coords = {
                # "frqs": ("frqs", frq),
                # "ffp": ("frqs", ffp)
            # },
            # data_vars = {
                # "Sxx": ("frqs", Sxx),
                # "Syy": ("frqs", Syy),
                # "Szz": ("frqs", Szz),
                # "Sxz": ("frqs", Sxz),
                # "Syz": ("frqs", Syz),
                # "Sxy": ("frqs", Sxy),
            # }
        # )

    # def directional_moments(self, Phi: xr.Dataset) -> xr.Dataset:
        # Exx, Eyy, Ezz, Cxy, Qxz, Qyz = (
            # np.real(Phi["Sxx"]), np.real(Phi["Syy"]), np.real(Phi["Szz"]),
            # np.real(Phi["Sxy"]), np.imag(Phi["Sxz"]), np.imag(Phi["Syz"])
        # )

        # return  xr.Dataset(
            # {
                # "a1": Qxz / np.sqrt(Ezz * (Exx + Eyy)),
                # "b1": Qyz / np.sqrt(Ezz * (Exx + Eyy)),
                # "a2": (Exx - Eyy) / (Exx + Eyy),
                # "b2": 2 * Cxy / (Exx + Eyy)
            # }
        # )

    # def dftm_distribution(self, moments: xr.Dataset):
        # dirs =  xr.Variable(dims=("dirs"), data=np.arange(-180,180,5))

        # D = (
            # 1/2 +
            # moments["a1"] * np.cos(np.radians(dirs)) +
            # moments["b1"] * np.sin(np.radians(dirs)) +
            # moments["a2"] * np.cos(2*np.radians(dirs)) +
            # moments["b2"] * np.sin(2*np.radians(dirs))
        # )
        # D.coords["dirs"] = dirs
        # D = D.where(D > 0, 0)
        # return D / D.integrate("dirs")

    # def mem_distribution(self, Phi: xr.Dataset):
        # dirs =  xr.Variable(dims=("dirs"), data=np.arange(-180,180,5))

        # c1 = Phi["a1"] + 1j*Phi["b1"]
        # c2 = Phi["a2"] + 1j*Phi["b2"]

        # phi1 = (c1 - c2 * c1.conj()) / (1 - np.abs(c1)**2)
        # phi2 = c2 - c1.conj() * phi1

        # sigma_e = 1 - phi1 * c1.conj() - phi2 * c2.conj()

        # D_mem = (1/(2*np.pi)) * np.real(
            # sigma_e.expand_dims({"dirs": dirs}) /
            # np.abs(
                # 1 - phi1.expand_dims({"dirs": dirs}) * np.exp(-1j*dirs*np.pi/180)
                  # - phi2.expand_dims({"dirs": dirs}) * np.exp(-2j*dirs*np.pi/180)
            # )**2
        # )

        # return D_mem / D_mem.integrate("dirs")



# if True:
    # D_dft = dftm_distribution(circ) 
    # D_mem = mem_distribution(circ)
    # D_dft_smo = D_dft.rolling(frqs=32, center=True).median()
    # D_mem_smo = D_mem.rolling(frqs=32, center=True).median()

    # D_mem_smo = xr.apply_ufunc(
        # gaussian_filter, D_mem, input_core_dims=[["dirs", "frqs"]],
        # output_core_dims=[["dirs", "frqs"]], kwargs={"sigma": 1}
    # )

    # plot_directional_spectrum(
        # D_mem_smo.T.pipe(np.log10), levels=30, colorbar=True,
        # axes_kw={"rmax": 0.5, "is_period": True}, vmin=-3, vmax=-1
    # )

    




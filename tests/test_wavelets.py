#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

"""
Copyright © 2024 Daniel Pelaez-Zapata <http://github.com/dspelaez>
Distributed under terms of the GNU/GPL 3.0 license.

@author: Daniel Pelaez-Zapata
@github: http://github.com/dspelaez
@created: 2024-06-27
"""

import numpy as np
import xarray as xr
import scipy.signal as signal

from ewdm.wavelets import Morlet, cwt, xwt
from ewdm.sources import SpotterBuoysDataSource


def test_morlet_initialisation():
    morlet = Morlet()
    assert morlet.f0 == 6
    assert morlet.cdelta == 0.776
    assert morlet.gamma == 2.32
    assert morlet.deltaj0 == 0.60


def test_morlet_properties():
    morlet = Morlet(6.)
    assert np.isclose(morlet.flambda, 1.03304)
    assert np.isclose(morlet.coi, 0.70711)
    assert np.isclose(morlet.sup, 1.41421)


def test_cwt():
    # create input data array
    time = np.arange(0., 1024., 1)
    data = np.sin(2 * np.pi * time / 10)
    data_array = xr.DataArray(data, dims=["time"], coords={"time": time})
    omin, omax, nvoice = -8, -2, 16
    freqs = 2. ** np.linspace(omin, omax, nvoice*abs(omin-omax)+1)
    result = cwt(data_array, freqs)
    assert isinstance(result, xr.DataArray)
    assert "sampling_rate" in result.attrs
    assert np.isclose(result.sampling_rate, 1.0)


def test_xwt():
    # create input data arrays
    time = np.arange(0., 1024., 1)
    data_x = np.sin(2 * np.pi * time / 15)
    data_y = np.sin(2 * np.pi * time / 5)
    data_x_array = xr.DataArray(data_x, dims=["time"], coords={"time": time})
    data_y_array = xr.DataArray(data_y, dims=["time"], coords={"time": time})
    omin, omax, nvoice = -8, -2, 16
    freqs = 2. ** np.linspace(omin, omax, nvoice*abs(omin-omax)+1)
    result = xwt(data_x_array, data_y_array, freqs)
    assert isinstance(result, xr.DataArray)


def test_cwt_spotter_data():

    dataset = (
        SpotterBuoysDataSource("../data/displacements.csv")
        .read_dataset()
    )
    
    omin, omax, nvoice = -5, 0, 16
    freqs = 2. ** np.linspace(omin, omax, nvoice*abs(omin-omax)+1)
    result = cwt(
        dataset["surface_elevation"], freqs, fs=2.5
    )
    power = np.abs(result) ** 2

    fourier_freq, fourier_power = signal.welch(
        dataset["surface_elevation"], fs=2.5, nperseg=512
    )

    fourier_std = np.trapz(fourier_power, x=fourier_freq) ** 0.5
    wavelet_std = power.mean("time").integrate("frequency").item() ** 0.5
    signal_std = dataset["surface_elevation"].std("time").item()

    assert isinstance(power, xr.DataAarray)
    assert np.allclose(wavelet_std, fourier_std, rtol=0.5)
    assert np.allclose(wavelet_std, signal_std, rtol=0.5)


def test_cwt_std():
    # create input data arrays
    for input_std in [0.01, 0.1, 1, 10, 100, 1000]:
        fs = 2
        time = np.arange(4500) / fs
        data = np.random.rand(len(time))
        data_array = xr.DataArray(
            input_std * (data - data.mean()) / data.std(),
            dims=["time"], coords={"time": time})
        omin, omax, nvoice = -5, 0, 16
        freqs = 2. ** np.linspace(omin, omax, nvoice*abs(omin-omax)+1)
        result = cwt(data_array, freqs, fs=fs)
        power = np.abs(result) ** 2

        fourier_freq, fourier_power = signal.welch(
            data_array, fs=2.5, nperseg=512
        )

        scipy_cwt = signal.cwt(
            data_array, widths=1/freqs, wavelet=signal.morlet2
        )
        scipy_power = np.abs(scipy_cwt) ** 2

        scipy_std = np.trapz(np.mean(scipy_power, axis=1), x=freqs) ** 0.5
        fourier_std = np.trapz(fourier_power, x=fourier_freq) ** 0.5
        wavelet_std = power.mean("time").integrate("frequency").item() ** 0.5
        signal_std = data_array.std("time").item()
    
        print("Wavelet std: ", wavelet_std)
        print("Scipy cwt std: ", scipy_std)
        print("Fourier std: ", fourier_std)
        print("Signal std: ", signal_std)
        print(35*"=")

        assert np.allclose(wavelet_std, scipy_std, rtol=0.1)
        assert np.allclose(wavelet_std, fourier_std, rtol=0.1)
        assert np.allclose(wavelet_std, signal_std, rtol=0.1)

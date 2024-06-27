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

from ewdm.wavelets import Morlet, cwt, xwt
# from ewdm.sources import SpotterBuoysDataSource, CDIPDataSourceRealTime

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
    # Create input data array
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
    # Create input data arrays
    time = np.arange(0., 1024., 1)
    data_x = np.sin(2 * np.pi * time / 15)
    data_y = np.sin(2 * np.pi * time / 5)
    data_x_array = xr.DataArray(data_x, dims=["time"], coords={"time": time})
    data_y_array = xr.DataArray(data_y, dims=["time"], coords={"time": time})
    omin, omax, nvoice = -8, -2, 16
    freqs = 2. ** np.linspace(omin, omax, nvoice*abs(omin-omax)+1)
    result = xwt(data_x_array, data_y_array, freqs)
    assert isinstance(result, xr.DataArray)

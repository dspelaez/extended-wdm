#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

"""
@author: Daniel Santiago
@github: http://github.com/dspelaez
@created: 2024-06-27
"""

import numpy as np
import xarray as xr
import pytest
import os

from datetime import datetime

from ewdm.sources import SpotterBuoysDataSource, CDIPDataSourceRealTime

HERE = os.path.dirname(os.path.abspath(__file__))
SPOTTER_FILE = os.path.join(HERE, "../data/displacement.csv")

@pytest.fixture
def spotter_instance():
    return SpotterBuoysDataSource(SPOTTER_FILE)

@pytest.fixture
def cdip_instance():
    return CDIPDataSourceRealTime(station_id=189)

def test_spotter_wrong_filename():
    with pytest.raises(FileNotFoundError):
        spotter = SpotterBuoysDataSource("fake_displacement.csv")

def test_spotter_read_dataset(spotter_instance):
    dataset = spotter_instance.read_dataset()
    assert "eastward_displacement" in dataset.variables
    assert "northward_displacement" in dataset.variables
    assert "surface_elevation" in dataset.variables
    assert "sampling_rate" in dataset.attrs

def test_cdip_read_dataset(cdip_instance):
    dataset = cdip_instance.read_dataset(time_start="2024-06-22 08:00")
    assert "eastward_displacement" in dataset.variables
    assert "northward_displacement" in dataset.variables
    assert "surface_elevation" in dataset.variables
    assert "sampling_rate" in dataset.attrs

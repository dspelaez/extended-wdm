#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

"""
Using real time CDIP data
=========================

This example shows how to compute the directional wave spectrum from waverider wave buoys data obtained from CDIP dataset in real time.
"""


import numpy as np
import xarray as xr

from matplotlib import pyplot as plt

import ewdm
from ewdm.sources import CDIPDataSourceRealTime
from ewdm.plots import plot_directional_spectrum

cdip =  CDIPDataSourceRealTime(166)
dataset = cdip.read_dataset(time_start='2024-06-09T08:30')
spec = ewdm.Triplets(dataset)
output = spec.compute(
    omin=-5, omax=0, nvoice=16, dd=5, kappa=36, use="displacements"
)

fig, ax = plt.subplots()
plot_directional_spectrum(
    output.directional_spectrum, ax=ax, levels=None, colorbar=True
)

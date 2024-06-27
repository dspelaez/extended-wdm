#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

"""
Copyright © 2024 Daniel Santiago <http://github.com/dspelaez>
Distributed under terms of the GNU/GPL 3.0 license.

@author: Daniel Santiago
@github: http://github.com/dspelaez
@created: 2024-03-12
"""

import numpy as np
import xarray as xr


def get_sampling_frequency(time, precision=3):
    """Try to figure out the sampling frequency for a dataset"""

    if np.issubdtype(time.dtype, np.datetime64):
        dt = (
                time
                .diff("time")
                .mean("time")
                .astype("timedelta64[ns]")
                .astype('float')
                .item()
            )
        return np.round(1E9 / dt, precision)

    elif (
            np.issubdtype(time.dtype, np.floating) or 
            np.issubdtype(time.dtype, np.integer)
        ):
        # assume to be seconds
        return 1/time.diff("time").mean("time").item()


def wavenumber(f, d=100):
    """Compute wavenumber according to Hunt approximation"""

    if d < 0:
        raise ValueError("Depth must be positive")

    d0 = [0.666, 0.355, 0.161, 0.0632, 0.0218, 0.0065]
    g = 9.8
    w = 2.* np.pi * f
    y = (w**2)*d/g
    #
    poly = np.zeros_like(f)
    for n, dn in enumerate(d0):
        poly = poly + dn * y**(n+1)
    #
    return np.sqrt(y**2 + y/(1 + poly))/d

def phase_speed(T, d=100):
    """Return phase speed for the given period"""
    return (2*np.pi/T) / wavenumber(1/T, d=d)


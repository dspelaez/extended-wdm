#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

import numpy as np
import xarray as xr


def get_sampling_frequency(time: xr.DataArray, precision: int = 3) -> float:
    """Try to figure out the sampling frequency for a time coordinate.

    Args:
        time (xr.DataArray): Time data array.
        precision (int, optional): Number of decimal places to round to. Defaults to 3.

    Returns:
        float: Estimated sampling frequency.
    """
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


def wavenumber(f: np.ndarray, d: float = 100) -> np.ndarray:
    """Compute wavenumber according to Hunt approximation.

    Args:
        f (np.ndarray): Frequency array.
        d (float, optional): Depth. Defaults to 100.

    Returns:
        np.ndarray: Computed wavenumber array.

    Raises:
        ValueError: If depth is negative.
    """

    if d < 0:
        raise ValueError("Depth must be positive")

    d0 = [0.666, 0.355, 0.161, 0.0632, 0.0218, 0.0065]
    g = 9.8
    w = 2. * np.pi * f
    y = (w ** 2) * d / g
    #
    poly = np.zeros_like(f)
    for n, dn in enumerate(d0):
        poly = poly + dn * y ** (n + 1)
    #
    return np.sqrt(y ** 2 + y / (1 + poly)) / d

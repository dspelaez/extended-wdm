#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8


import numpy as np
import xarray as xr
import logging

from .parameters import VARIABLE_NAMES

# von mises kernel desnsity estimation
def _vonmises_kde(arr, bins, kappa):
    """Return Kernel-Density Estimation using von Mises distribution"""

    # define input parameters
    x = np.radians(bins[:,None] - arr[None,:])

    # integrate vonmises kernels
    kde = (
        np.exp(kappa * np.cos(x)).sum(axis=1) / (2 * np.pi * np.i0(kappa))
    )
    kde /= np.trapz(kde, x=bins)
    return kde


# function to get histogram along freq axis
def _get_density(arr, bins, kappa):
    if np.isnan(arr).all():
        return np.zeros_like(bins, dtype="float") * np.nan
    else:
        if kappa is not None:
            return _vonmises_kde(arr, bins=bins, kappa=kappa)
        else:
            bins_edges = np.r_[bins, bins[-1]+np.diff(bins)[0]]
            return np.histogram(arr, bins=bins_edges, density=True)[0]


# function to actually estimate the spectrum
def estimate_directional_distribution(power, theta, dd, kappa):
    """Construct directional distribution function from local wave directions"""

    # array of directiontions where dd is the directional resolution
    bins = np.arange(-180, 180, dd)

    # directional distribution function
    D = np.apply_along_axis(
        _get_density, arr=theta, bins=bins, kappa=kappa, axis=1
    )

    # determine average wavelet power
    S = power.mean("time").data

    # array containing directional wave spectra
    E = S[:,None] * D

    # return dataset
    output_ds = xr.Dataset(
        data_vars = {
            "directional_spectrum": (["frequency", "direction"], E.data),
            "directional_distribution": (["frequency", "direction"], D.data),
            "frequency_spectrum": (["frequency"], S)
        },
        coords = {
            "frequency": power["frequency"].data,
            "direction": bins,
        }
    )
    output_ds["frequency"].attrs = VARIABLE_NAMES["frequency"]
    output_ds["direction"].attrs = VARIABLE_NAMES["direction"]
    for var in output_ds:
        output_ds[var].attrs = VARIABLE_NAMES[var]

    return output_ds


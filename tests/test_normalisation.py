#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

"""
Copyright © 2024 Daniel Pelaez-Zapata <http://github.com/dspelaez>
Distributed under terms of the GNU/GPL 3.0 license.

Regression tests for the wavelet-power normalisation in `Triplets` and
`Arrays`. When `normalise=True`, the frequency spectrum must integrate to
the surface-elevation variance, int S(f) df = sigma^2. The correct rescale
factor is variance / wavelet-energy (sigma^2 / E), not sigma / sqrt(E),
which only coincides when the raw wavelet integral E already equals sigma^2.
"""

import numpy as np
import pytest

import ewdm

# numpy 2.0 renamed np.trapz to np.trapezoid; support both.
_trapezoid = getattr(np, "trapezoid", getattr(np, "trapz", None))

GRAV = 9.8

# a few octaves that comfortably bracket the synthetic wave frequencies
COMMON = dict(omin=-4, omax=1, nvoice=16)


def _synthetic_sea_state(fs=4.0, duration=1200.0, seed=0):
    """A multi-component synthetic sea state on a Triplets-friendly grid.

    Returns time, surface elevation and the eastward/northward
    displacements of a set of plane waves travelling in different
    directions, so the variance is well defined and non-trivial.
    """
    rng = np.random.default_rng(seed)
    t = np.arange(0, duration, 1 / fs)

    components = [
        # (frequency Hz, amplitude m, direction deg)
        (0.12, 0.60, 30.0),
        (0.20, 0.35, 75.0),
        (0.33, 0.15, 200.0),
    ]

    eta = np.zeros_like(t)
    east = np.zeros_like(t)
    north = np.zeros_like(t)
    for f0, amp, drc in components:
        omega = 2 * np.pi * f0
        phase = rng.uniform(0, 2 * np.pi)
        arg = omega * t + phase
        theta = np.radians(drc)
        eta += amp * np.cos(arg)
        # deep-water: horizontal displacement leads elevation by 90 deg,
        # projected onto the direction of propagation
        east += amp * np.cos(theta) * np.sin(arg)
        north += amp * np.sin(theta) * np.sin(arg)

    return t, eta, east, north


def _plane_wave_array(f0=0.2, amplitude=0.5, direction_deg=40.0,
                      fs=4.0, duration=1200.0):
    """A single plane wave sampled by a small pentagon-plus-centre array."""
    t = np.arange(0, duration, 1 / fs)
    omega = 2 * np.pi * f0
    k0 = omega**2 / GRAV
    theta = np.radians(direction_deg)
    kx0, ky0 = k0 * np.cos(theta), k0 * np.sin(theta)

    radius = 0.5 * (np.pi / k0) / 2
    ang = np.radians(np.array([0, 72, 144, 216, 288]))
    px = np.r_[0.0, radius * np.cos(ang)]
    py = np.r_[0.0, radius * np.sin(ang)]

    eta = np.stack(
        [amplitude * np.cos(kx0 * x + ky0 * y - omega * t)
         for x, y in zip(px, py)],
        axis=-1
    )
    return t, eta, px, py


def _integrated_variance(spectrum):
    """Integrate a frequency_spectrum DataArray over frequency."""
    return _trapezoid(spectrum.data, x=spectrum["frequency"].data)


def test_triplets_normalisation_conserves_variance():
    """int S(f) df must equal the surface-elevation variance."""
    t, eta, east, north = _synthetic_sea_state()
    sigma2 = eta.var()

    trip = ewdm.Triplets.from_numpy(
        time=t, surface_elevation=eta,
        eastward_displacement=east, northward_displacement=north,
        fs=4.0, normalise=True, interpolate=False
    )
    out = trip.compute(use="displacements", **COMMON)

    integrated = _integrated_variance(out["frequency_spectrum"])
    assert np.isclose(integrated, sigma2, rtol=1e-6)


def test_arrays_normalisation_conserves_variance():
    """int S(f) df must equal the surface-elevation variance."""
    t, eta, px, py = _plane_wave_array()
    sigma2 = eta.var(axis=0).mean()  # per-element variance, averaged

    arr = ewdm.Arrays.from_numpy(
        time=t, surface_elevation=eta,
        position_x=px, position_y=py, fs=4.0, normalise=True, interpolate=False
    )
    out = arr.compute(coordinate="frequency", **COMMON)

    integrated = _integrated_variance(out["frequency_spectrum"])
    assert np.isclose(integrated, sigma2, rtol=1e-6)


def test_normalise_flag_changes_amplitude():
    """Turning normalisation off must change the spectral amplitude.

    Guards against the rescale factor silently collapsing to 1, which is
    what masked the original sigma / sqrt(E) bug.
    """
    t, eta, east, north = _synthetic_sea_state()

    common = dict(use="displacements", **COMMON)
    on = ewdm.Triplets.from_numpy(
        time=t, surface_elevation=eta,
        eastward_displacement=east, northward_displacement=north,
        fs=4.0, normalise=True, interpolate=False
    ).compute(**common)
    off = ewdm.Triplets.from_numpy(
        time=t, surface_elevation=eta,
        eastward_displacement=east, northward_displacement=north,
        fs=4.0, normalise=False, interpolate=False
    ).compute(**common)

    var_on = _integrated_variance(on["frequency_spectrum"])
    var_off = _integrated_variance(off["frequency_spectrum"])
    assert not np.isclose(var_on, var_off, rtol=1e-3)
    # only the normalised spectrum is tied to the physical variance
    assert np.isclose(var_on, eta.var(), rtol=1e-6)

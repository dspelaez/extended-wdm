#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

"""
Copyright © 2024 Daniel Pelaez-Zapata <http://github.com/dspelaez>
Distributed under terms of the GNU/GPL 3.0 license.

Tests for the multi-aperture estimator `MultiApertureArrays`. Both input
modes are exercised: the dense elevation field, from which virtual staffs are
seeded into nested apertures, and the sparse sensor array, whose sub-apertures
are grouped by baseline. The wavenumber and inverse-phase-speed (nu) spectra
follow Björkqvist et al. (2019).
"""

import numpy as np
import pytest

import ewdm
from ewdm import MultiApertureArrays
from ewdm.multiaperture import _cw_from_N_to_math_angle

# numpy 2.0 renamed np.trapz to np.trapezoid; support both in tests.
_trapezoid = getattr(np, "trapezoid", getattr(np, "trapz", None))


GRAV = 9.8


def _wavenumber(f0, depth):
    """Finite-depth linear wavenumber, solved iteratively."""
    omega = 2 * np.pi * f0
    k = omega**2 / GRAV
    for _ in range(100):
        k = omega**2 / (GRAV * np.tanh(k * depth))
    return k


def _monochromatic_field(f0=0.25, direction_deg=40.0, amplitude=0.5,
                         dx=1.0, fs=4.0, ny=48, nx=48, duration=128.0,
                         depth=100.0):
    """Build a dense plane-wave elevation field.

    Returns the field together with the true wavenumber and inverse phase speed
    so the tests can check the peak locations. The window is a few wavelengths
    wide so the coarse apertures alias the wave, which is the regime the
    anti-alias gate is meant to handle.
    """
    omega = 2 * np.pi * f0
    k0 = _wavenumber(f0, depth)
    theta = np.radians(direction_deg)
    kx0, ky0 = k0 * np.cos(theta), k0 * np.sin(theta)

    T = int(duration * fs)
    t = np.arange(T) / fs
    # build the wave in the same East/North frame as `seed_aperture`
    # (East = (cxp - j) dx, North = (cyp - i) dx), so the planted propagation
    # direction is exactly `direction_deg` in math angle
    cxp, cyp = (nx - 1) / 2, (ny - 1) / 2
    east = (cxp - np.arange(nx))[None, :, None] * dx
    north = (cyp - np.arange(ny))[:, None, None] * dx
    eta = amplitude * np.cos(kx0 * east + ky0 * north - omega * t[None, None, :])
    return eta, k0, k0 / omega


def _ring_array(rings=(3.0, 6.0, 10.0, 15.0), n_azimuth=6):
    """Concentric-ring sparse array with a centre element."""
    px, py = [0.0], [0.0]
    for r in rings:
        for a in np.linspace(0, 2 * np.pi, n_azimuth, endpoint=False):
            px.append(r * np.cos(a))
            py.append(r * np.sin(a))
    return np.array(px), np.array(py)


def _sparse_monochromatic(f0=0.12, direction_deg=30.0, amplitude=0.4,
                          fs=4.0, duration=256.0, depth=100.0):
    """Plane wave sampled at the concentric-ring array elements."""
    px, py = _ring_array()
    omega = 2 * np.pi * f0
    k0 = _wavenumber(f0, depth)
    theta = np.radians(direction_deg)
    kx0, ky0 = k0 * np.cos(theta), k0 * np.sin(theta)
    T = int(duration * fs)
    t = np.arange(T) / fs
    elev = np.stack(
        [amplitude * np.cos(kx0 * x + ky0 * y - omega * t)
         for x, y in zip(px, py)],
        axis=-1
    )
    return t, elev, px, py, k0, k0 / omega


# common compute settings: an explicit frequency grid below the Nyquist and the
# anti-alias gate on, as the dense field aliases the short test wave
FREQS = np.logspace(np.log10(0.06), np.log10(1.9), 52)
DENSE = dict(seed=1, n_staff=12, freqs=FREQS, antialias_gate=True)


def test_field_schema_and_units():
    """The dense path should return the documented variables and units."""
    eta, _, _ = _monochromatic_field()
    spec = MultiApertureArrays.from_field(eta, dx=1.0, depth=100.0, fs=4.0)
    out = spec.compute(**DENSE)

    for var in ("frequency_spectrum", "wavenumber_spectrum", "nu_spectrum",
                "directional_spectrum_f", "directional_spectrum_k",
                "directional_spectrum_nu", "mean_direction",
                "directional_spread", "sign_reference", "var_eta"):
        assert var in out

    assert out["directional_spectrum_f"].dims == ("frequency", "direction")
    assert out["directional_spectrum_k"].dims == ("wavenumber", "direction")
    assert out["directional_spectrum_nu"].dims == ("nu", "direction")
    assert out["wavenumber"].attrs["units"] == "rad/m"
    assert out["nu"].attrs["units"] == "s/m"
    assert out["directional_spectrum_k"].attrs["units"] == "m^4/rad"
    assert out["directional_spectrum_nu"].attrs["units"] == "m^4/(s^2 rad)"
    assert out.attrs["method"] == "multiaperture EWDM"
    assert out.attrs["mode"] == "field"


def test_field_variance_calibrated():
    """The integral of S(f) over frequency should equal the field variance."""
    eta, _, _ = _monochromatic_field()
    out = MultiApertureArrays.from_field(eta, dx=1.0, depth=100.0,
                                         fs=4.0).compute(**DENSE)
    var_eta = float(out["var_eta"])
    m0 = _trapezoid(out["frequency_spectrum"].values, out["frequency"].values)
    assert np.isclose(m0, var_eta, rtol=1e-6)
    # the plane-wave variance is amplitude**2 / 2
    assert np.isclose(var_eta, 0.5**2 / 2, rtol=0.05)


def test_field_wavenumber_peak():
    """The wavenumber spectrum should peak near the true wavenumber."""
    eta, k0, _ = _monochromatic_field()
    out = MultiApertureArrays.from_field(eta, dx=1.0, depth=100.0,
                                         fs=4.0).compute(**DENSE)
    k = out["wavenumber"].values
    Fk = out["wavenumber_spectrum"].values
    k_peak = k[np.nanargmax(Fk)]
    assert abs(k_peak - k0) / k0 < 0.15


def test_field_nu_peak():
    """The nu spectrum should peak near k0 / omega0 with the gates disabled."""
    eta, _, nu0 = _monochromatic_field()
    out = MultiApertureArrays.from_field(eta, dx=1.0, depth=100.0, fs=4.0).compute(
        nu_f_lim=None, nu_k_lim=None, **DENSE)
    nu = out["nu"].values
    Q = out["nu_spectrum"].values
    nu_peak = nu[np.nanargmax(Q)]
    assert abs(nu_peak - nu0) / nu0 < 0.15


def test_direction_conventions_consistent():
    """The math and compass conventions must be an exact transform of one
    another, and both must recover the planted direction."""
    eta, _, _ = _monochromatic_field(direction_deg=40.0)
    out_cw = MultiApertureArrays.from_field(
        eta, dx=1.0, depth=100.0, fs=4.0,
        direction_convention="cw_from_north").compute(**DENSE)
    out_math = MultiApertureArrays.from_field(
        eta, dx=1.0, depth=100.0, fs=4.0,
        direction_convention="math").compute(**DENSE)

    fpk = np.nanargmax(out_cw["frequency_spectrum"].values)
    dir_cw = out_cw["mean_direction"].values[fpk]
    dir_math = out_math["mean_direction"].values[fpk]

    # math angle 40 deg corresponds to 50 deg clockwise from North
    assert abs(((dir_cw - 50.0 + 180) % 360) - 180) < 12.0
    # the two conventions are the exact analytic transform of each other
    assert np.isclose(_cw_from_N_to_math_angle(dir_cw), dir_math, atol=1e-6)


def test_antialias_gate_removes_lowk_pileup():
    """Without the gate the aliased coarse apertures pile energy at low k; the
    gate moves the peak back onto the true wavenumber."""
    eta, k0, _ = _monochromatic_field()
    spec = MultiApertureArrays.from_field(eta, dx=1.0, depth=100.0, fs=4.0)
    # disable the reliability gate (default on) to isolate the anti-alias gate
    open_ = spec.compute(seed=1, n_staff=12, freqs=FREQS,
                         antialias_gate=False, reliability_gate=None)
    gated = spec.compute(seed=1, n_staff=12, freqs=FREQS,
                         antialias_gate=True, reliability_gate=None)

    k = gated["wavenumber"].values
    k_open = k[np.nanargmax(open_["wavenumber_spectrum"].values)]
    k_gated = k[np.nanargmax(gated["wavenumber_spectrum"].values)]
    assert k_open < 0.5 * k0          # spurious low-k peak without the gate
    assert abs(k_gated - k0) / k0 < 0.15


def test_reliability_gate_removes_lowk_pileup():
    """The dispersion-independent reliability gate removes the aliased low-k
    pileup just like the anti-alias gate, recovering the true peak wavenumber."""
    eta, k0, _ = _monochromatic_field()
    spec = MultiApertureArrays.from_field(eta, dx=1.0, depth=100.0, fs=4.0)
    open_ = spec.compute(seed=1, n_staff=12, freqs=FREQS, reliability_gate=None)
    gated = spec.compute(seed=1, n_staff=12, freqs=FREQS, reliability_gate=0.5)

    k = gated["wavenumber"].values
    k_open = k[np.nanargmax(open_["wavenumber_spectrum"].values)]
    k_gated = k[np.nanargmax(gated["wavenumber_spectrum"].values)]
    assert k_open < 0.5 * k0          # spurious low-k peak without any gate
    assert abs(k_gated - k0) / k0 < 0.15


def test_depiston_is_solve_only():
    """Built-in de-piston changes only the wavevector solve, so the frequency
    spectrum (and hence the variance) is left untouched."""
    eta, _, _ = _monochromatic_field()
    spec = MultiApertureArrays.from_field(eta, dx=1.0, depth=100.0, fs=4.0)
    base = spec.compute(**DENSE)
    dep = spec.compute(depiston=True, **DENSE)
    assert np.allclose(base["frequency_spectrum"].values,
                       dep["frequency_spectrum"].values)
    assert np.isfinite(dep["wavenumber_spectrum"].values).any()


def test_solver_misfit_low_for_clean_wave():
    """`solve_wavevectors` returns a misfit that is small for a cleanly resolved
    plane wave (consistent phase differences across pairs)."""
    from ewdm.multiaperture import solve_wavevectors, cwt_stack, seed_aperture
    eta, k0, _ = _monochromatic_field()
    ny, nx, _ = eta.shape
    ii, jj, px, py, _ = seed_aperture(ny, nx, 1.0, 10, 12, 0)
    es = np.stack([eta[i, j, :] for i, j in zip(ii, jj)])
    W = cwt_stack(es - es.mean(1, keepdims=True), FREQS, 4.0, 12.0)
    kx, ky, resid, misfit = solve_wavevectors(W, px, py)
    # at the wave frequency the misfit should be small (well-resolved)
    fi = np.argmin(np.abs(FREQS - 0.25))
    assert np.nanmedian(misfit[fi]) < 0.3


def test_stitch_no_interior_notch():
    """The stitched F(k) must have no NaN or zero notch inside the band actually
    covered by the apertures (where the stitch weight is non-zero)."""
    eta, _, _ = _monochromatic_field()
    out = MultiApertureArrays.from_field(eta, dx=1.0, depth=100.0, fs=4.0).compute(
        return_apertures=True, **DENSE)
    Fk = out["wavenumber_spectrum"].values
    ck = out["stitch_weight_k"].values
    covered = ck > 1e-6
    assert np.all(np.isfinite(Fk[covered]))
    assert np.all(Fk[covered] > 0)


def test_return_apertures_diagnostics():
    """The per-aperture diagnostics should appear on an aperture dimension."""
    eta, _, _ = _monochromatic_field()
    out = MultiApertureArrays.from_field(eta, dx=1.0, depth=100.0, fs=4.0).compute(
        return_apertures=True, **DENSE)
    assert "aperture" in out.dims
    assert out["aperture_bmax"].dims == ("aperture",)
    assert out["aperture_Fk"].dims == ("aperture", "wavenumber")
    # the baselines should shrink from the coarsest to the tightest aperture
    bmax = out["aperture_bmax"].values
    assert bmax[0] > bmax[-1]


def test_sparse_wavenumber_peak():
    """The sparse path should recover the peak wavenumber of a plane wave."""
    t, elev, px, py, k0, _ = _sparse_monochromatic()
    spec = MultiApertureArrays.from_numpy(
        time=t, surface_elevation=elev, position_x=px, position_y=py,
        depth=100.0, fs=4.0)
    k_grid = 2.0**np.linspace(np.log2(0.02), np.log2(1.0), 64)
    nu_grid = 2.0**np.linspace(np.log2(0.02), np.log2(2.0), 64)
    out = spec.compute(freqs=FREQS, k_grid=k_grid, nu_grid=nu_grid,
                       nu_f_lim=None, nu_k_lim=None)

    assert out.attrs["mode"] == "sparse"
    k = out["wavenumber"].values
    k_peak = k[np.nanargmax(out["wavenumber_spectrum"].values)]
    assert abs(k_peak - k0) / k0 < 0.2


def test_sparse_collinear_raises():
    """A collinear array cannot form a two-dimensional aperture."""
    t = np.arange(0, 256, 0.25)
    px = np.array([0.0, 1.0, 2.0, 3.0])
    py = np.zeros_like(px)
    elev = np.stack([np.cos(0.1 * x + t) for x in px], axis=-1)
    spec = MultiApertureArrays.from_numpy(
        time=t, surface_elevation=elev, position_x=px, position_y=py,
        depth=100.0, fs=4.0)
    with pytest.raises(ValueError):
        spec.compute(freqs=FREQS)


def test_sparse_apertures_nested():
    """The sparse sub-apertures should be nested, coarse to tight, each with at
    least three non-collinear elements."""
    px, py = _ring_array()
    t = np.arange(0, 64, 0.25)
    elev = np.zeros((len(t), len(px)))
    spec = MultiApertureArrays.from_numpy(
        time=t, surface_elevation=elev, position_x=px, position_y=py,
        depth=100.0, fs=4.0)
    aps = spec._sparse_apertures()
    assert len(aps) >= 2
    bmax = []
    for _, idx in aps:
        assert len(idx) >= 3
        sx, sy = px[idx], py[idx]
        bmax.append(np.hypot(sx[:, None] - sx[None, :],
                             sy[:, None] - sy[None, :]).max())
    # coarse -> tight: the longest baseline shrinks from one aperture to the next
    assert all(bmax[i] >= bmax[i + 1] for i in range(len(bmax) - 1))


def test_sparse_lh_sign_with_slopes():
    """Co-located slopes let the sparse path resolve the 180-degree ambiguity
    through the Longuet-Higgins first moment."""
    import xarray as xr
    px, py = _ring_array()
    f0, direction_deg, depth, fs = 0.12, 30.0, 100.0, 4.0
    omega = 2 * np.pi * f0
    k0 = _wavenumber(f0, depth)
    theta = np.radians(direction_deg)
    kx0, ky0 = k0 * np.cos(theta), k0 * np.sin(theta)
    t = np.arange(int(256 * fs)) / fs
    ph = kx0 * px[:, None] + ky0 * py[:, None] - omega * t[None, :]
    elev = (0.4 * np.cos(ph)).T                       # (time, element)
    se = (-0.4 * kx0 * np.sin(ph)).T                  # eastward slope
    sn = (-0.4 * ky0 * np.sin(ph)).T                  # northward slope
    ds = xr.Dataset(
        data_vars={
            "surface_elevation": (["time", "element"], elev),
            "eastward_slope": (["time", "element"], se),
            "northward_slope": (["time", "element"], sn),
            "position_x": ("element", px),
            "position_y": ("element", py),
        },
        coords={"time": t, "element": np.arange(len(px))},
    )
    spec = MultiApertureArrays.from_arrays(ds, depth=depth, fs=fs,
                                           direction_convention="math")
    k_grid = 2.0**np.linspace(np.log2(0.02), np.log2(1.0), 64)
    nu_grid = 2.0**np.linspace(np.log2(0.02), np.log2(2.0), 64)
    out = spec.compute(freqs=FREQS, k_grid=k_grid, nu_grid=nu_grid,
                       nu_f_lim=None, nu_k_lim=None)
    # the Longuet-Higgins reference is populated (sign path ran)
    assert np.isfinite(out["lh_direction"].values).any()
    # the peak-frequency direction lands in the planted hemisphere
    fpk = np.nanargmax(out["frequency_spectrum"].values)
    md = out["mean_direction"].values[fpk]
    assert np.cos(np.radians(md - direction_deg)) > 0


def test_sparse_reliability_gate():
    """The reliability gate (default on) works through the sparse path and still
    recovers the peak wavenumber."""
    t, elev, px, py, k0, _ = _sparse_monochromatic()
    spec = MultiApertureArrays.from_numpy(
        time=t, surface_elevation=elev, position_x=px, position_y=py,
        depth=100.0, fs=4.0)
    k_grid = 2.0**np.linspace(np.log2(0.02), np.log2(1.0), 64)
    nu_grid = 2.0**np.linspace(np.log2(0.02), np.log2(2.0), 64)
    out = spec.compute(freqs=FREQS, k_grid=k_grid, nu_grid=nu_grid,
                       nu_f_lim=None, nu_k_lim=None)   # default reliability_gate
    k = out["wavenumber"].values
    assert abs(k[np.nanargmax(out["wavenumber_spectrum"].values)] - k0) / k0 < 0.25


def test_sign_anchor_flips_hemisphere():
    """An external sign anchor pointing the other way must flip the reported
    hemisphere by roughly 180 degrees."""
    eta, _, _ = _monochromatic_field(direction_deg=40.0)
    spec = MultiApertureArrays.from_field(eta, dx=1.0, depth=100.0, fs=4.0)
    base = spec.compute(**DENSE)

    fpk = np.nanargmax(base["frequency_spectrum"].values)
    d0 = base["mean_direction"].values[fpk]

    # anchor pointing 180 degrees away over the whole band, fully reliable
    anchor = (FREQS, np.full_like(FREQS, (d0 + 180.0) % 360.0),
              np.ones_like(FREQS))
    flipped = spec.compute(sign_anchor=anchor, **DENSE)
    d1 = flipped["mean_direction"].values[fpk]

    sep = abs(((d1 - d0 + 180) % 360) - 180)
    assert sep > 120.0


def test_backward_compatibility_arrays():
    """Adding the multi-aperture class must not disturb the existing Arrays
    radial path."""
    t = np.arange(0, 256, 0.25)
    f0, omega = 0.2, 2 * np.pi * 0.2
    k0 = omega**2 / GRAV
    radius = 0.5 * (np.pi / k0) / 2
    ang = np.radians(np.array([0, 72, 144, 216, 288]))
    px = np.r_[0.0, radius * np.cos(ang)]
    py = np.r_[0.0, radius * np.sin(ang)]
    eta = np.stack([0.5 * np.cos(k0 * np.cos(0.7) * x + k0 * np.sin(0.7) * y
                                 - omega * t) for x, y in zip(px, py)], axis=-1)
    arr = ewdm.Arrays.from_numpy(time=t, surface_elevation=eta,
                                 position_x=px, position_y=py, fs=4.0)
    out = arr.compute(coordinate="wavenumber", omin=-4, omax=1, nvoice=8)
    assert "wavenumber_spectrum" in out


def test_auto_apertures_default_ladder():
    """The default (no `apertures`) builds a data-aware ladder of coarse windows
    down to tight crosses, spanning the field of view to near the pixel scale."""
    from ewdm.multiaperture import auto_apertures
    ny, nx, dx = 64, 64, 0.5
    aps = auto_apertures(ny, nx, dx)
    bl = [(2 * e[1] * dx if isinstance(e, tuple) else e * dx) for _, e in aps]
    assert all(a > b for a, b in zip(bl, bl[1:]))            # coarse -> fine
    assert not isinstance(aps[0][1], tuple)                  # coarsest is a window
    assert bl[0] <= min(ny, nx) * dx                         # fits the field
    assert bl[-1] <= 6 * dx                                  # tight high-k rung
    assert any(isinstance(e, tuple) and e[0] == "plus" for _, e in aps)

    # compute() with no apertures runs and the ladder reaches near the pixel scale
    eta, k0, _ = _monochromatic_field(dx=dx, ny=ny, nx=nx)
    out = MultiApertureArrays.from_field(eta, dx, depth=100.0, fs=4.0).compute(
        seed=0, return_apertures=True)
    assert out["aperture_khi"].values.max() > np.pi / (6 * dx)
    k = out["wavenumber"].values
    F = out["wavenumber_spectrum"].values
    assert abs(k[np.nanargmax(F)] - k0) / k0 < 0.35         # peak still recovered

#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

"""
Composite multi-aperture directional wave spectra.

Each wavenumber octave is resolved with a matched aperture: coarse, widely
spaced staffs for the long waves (low wavenumber) and tight clusters for the
short waves (high wavenumber). Every aperture is only trusted within its own
anti-alias band, k < pi / baseline_max, and the per-aperture spectra are then
stitched together. Two kinds of input are supported:

* a dense elevation field eta(y, x, t) together with the grid spacing dx, from
  which virtual staffs are seeded into nested apertures. This is built with the
  `from_field` constructor.

* a sparse sensor array following the same convention as the `Arrays` class,
  from which nested sub-apertures are grouped by baseline. This is built with
  the `from_arrays` (or `from_numpy`) constructor.

The estimator only deals with the array-measured (wavelet directional method)
machinery. The camera slope-to-elevation reconstruction and any external sign
anchor (for example a 3D FFT estimate) are left to the caller and passed in as
hooks, namely the `solve_eta` and `sign_anchor` arguments of `compute`.

Wave direction is reported in degrees clockwise from North (the oceanographic
convention) by default. Set `direction_convention="math"` to obtain instead the
eastward-northward math angle used elsewhere in ewdm.
"""

import logging

import numpy as np
import xarray as xr

from .wavelets import cwt, Morlet
from .density import (estimate_directional_distribution,
                      estimate_radial_distribution)
from .helpers import get_sampling_frequency
from .parameters import VARIABLE_NAMES, GRAV

logger = logging.getLogger(__name__)

# numpy 2.0 renamed `np.trapz` to `np.trapezoid`; use the new name when
# available, fall back on older numpy. mirrors the shim in `density.py`.
_trapezoid = getattr(np, "trapezoid", getattr(np, "trapz", None))

RADTODEG = 180. / np.pi


# dispersion relation and circular statistics {{{
def _nu_of_k(k, depth):
    """Inverse phase speed nu = k / omega along the linear dispersion relation."""
    return np.sqrt(k / (GRAV * np.tanh(np.clip(k * depth, 1e-9, 50.0))))


def k_dispersion(f, depth):
    """Linear gravity-wave wavenumber k(f) in rad/m at finite depth.

    The dispersion relation is solved iteratively because the wavenumber appears
    on both sides through the hyperbolic tangent.
    """
    w = 2 * np.pi * np.asarray(f, float)
    k = w**2 / GRAV
    for _ in range(100):
        k_new = w**2 / (GRAV * np.tanh(np.clip(k * depth, 1e-9, 50.0)))
        converged = np.allclose(k_new, k, rtol=1e-12, atol=0.0)
        k = k_new
        if converged:
            break
    return k


def circ_stats(F, theta_deg, axis=-1):
    """Energy-weighted circular mean direction and spread in degrees.

    The mean direction and angular spread are computed from the first
    trigonometric moment of the directional spectrum F over the direction axis.
    """
    a = np.radians(theta_deg)
    C = np.sum(F * np.cos(a), axis=axis)
    S = np.sum(F * np.sin(a), axis=axis)
    tot = np.sum(F, axis=axis)
    tot = np.where(tot == 0, np.nan, tot)
    R = np.hypot(C, S) / tot
    mean = np.degrees(np.arctan2(S, C)) % 360.0
    spread = np.degrees(np.sqrt(np.clip(2.0 * (1.0 - R), 0, None)))
    return mean, spread


def _math_angle_to_cw_from_N(theta_en_rad, flip=False):
    """Convert an eastward-northward math angle (radians, counter-clockwise from
    East) into a compass direction in degrees (clockwise from North), wrapped to
    the [-180, 180) interval."""
    deg = (90.0 - np.degrees(theta_en_rad)) % 360.0
    if flip:
        deg = (deg + 180.0) % 360.0
    return ((deg + 180.0) % 360.0) - 180.0


def _cw_from_N_to_math_angle(deg_cwN):
    """Inverse of `_math_angle_to_cw_from_N`: convert a compass direction in
    degrees (clockwise from North) into the eastward-northward math angle in
    degrees, wrapped to the [-180, 180) interval."""
    a = (90.0 - np.asarray(deg_cwN, float)) % 360.0
    return ((a + 180.0) % 360.0) - 180.0


def _wrap180(x):
    """Wrap an angle in degrees to the [-180, 180) interval."""
    return ((x + 180.0) % 360.0) - 180.0
# }}}


# array geometry, aperture seeding and stitching {{{
def erode_valid(valid):
    """Erode a boolean footprint by one pixel.

    Only pixels whose four neighbours are all inside the footprint are kept, so
    that central-difference gradients stay within it.
    """
    ev = valid.copy()
    ev[1:-1, 1:-1] &= (valid[:-2, 1:-1] & valid[2:, 1:-1]
                       & valid[1:-1, :-2] & valid[1:-1, 2:])
    ev[[0, -1], :] = False
    ev[:, [0, -1]] = False
    return ev


def seed_aperture(ny, nx, dx, extent_px, n_staff, seed, valid=None):
    """Draw random virtual staffs inside a centred window.

    Args:
        ny, nx (int): Number of grid rows and columns.
        dx (float): Grid spacing in metres.
        extent_px (int or tuple): Window side in pixels. A scalar gives a square
            window, a (rows, cols) tuple a rectangular one.
        n_staff (int): Number of virtual staffs to draw.
        seed (int): Seed for the random number generator.
        valid (np.ndarray, optional): Boolean footprint. When given, staffs are
            drawn only from the True pixels, the window is recentred on the valid
            centroid and the staffs are sampled without replacement (with
            replacement only if fewer than `n_staff` pixels are available).

    Returns:
        tuple: The row and column indices of the staffs, their easting and
        northing positions in metres and the longest baseline. The column axis
        is West-positive, so East = (cxp - j) dx and North = (cyp - i) dx.
    """
    rng = np.random.default_rng(seed)
    cxp, cyp = (nx - 1) / 2, (ny - 1) / 2
    # deterministic plus (cross) aperture: centre + four arms at +/- arm_px, a
    # tight, well-conditioned high-wavenumber rung. arm_px may be fractional
    # (the field is bilinearly sampled when the staff series are extracted).
    if isinstance(extent_px, (tuple, list)) and len(extent_px) == 2 and extent_px[0] == "plus":
        a = float(extent_px[1])
        ii = np.array([cyp, cyp, cyp, cyp + a, cyp - a])
        jj = np.array([cxp, cxp + a, cxp - a, cxp, cxp])
        px = (cxp - jj) * dx; py = (cyp - ii) * dx
        b = np.hypot(px[:, None] - px[None, :], py[:, None] - py[None, :])
        return ii, jj, px, py, float(b.max())
    ey, ex = extent_px if isinstance(extent_px, (tuple, list)) else (extent_px, extent_px)
    ey = int(min(ey, ny - 1))
    ex = int(min(ex, nx - 1))

    # draw staffs either from the whole centred window or, when a footprint is
    # given, only from the valid pixels around the footprint centroid
    if valid is None:
        i0 = (ny - ey) // 2
        j0 = (nx - ex) // 2
        ii = rng.integers(i0, i0 + ey + 1, n_staff)
        jj = rng.integers(j0, j0 + ex + 1, n_staff)
    else:
        vi, vj = np.where(valid)
        ci, cj = int(round(vi.mean())), int(round(vj.mean()))
        i0, j0 = max(ci - ey // 2, 0), max(cj - ex // 2, 0)
        wi, wj = np.where(valid[i0:i0 + ey + 1, j0:j0 + ex + 1])
        # if the window missed the footprint, fall back to all valid pixels
        if wi.size == 0:
            wi, wj, i0, j0 = vi, vj, 0, 0
        pick = rng.choice(wi.size, size=n_staff, replace=wi.size < n_staff)
        ii, jj = wi[pick] + i0, wj[pick] + j0

    # convert pixel indices to physical positions and find the longest baseline
    px = (cxp - jj) * dx
    py = (cyp - ii) * dx
    b = np.hypot(px[:, None] - px[None, :], py[:, None] - py[None, :])
    return ii, jj, px, py, float(b.max())


def _bilinear_stack(arr, ii, jj):
    """Sample arr (ny, nx, T) at (possibly fractional) staff rows ii and columns
    jj by bilinear interpolation. Reduces to plain indexing for integer ii, jj."""
    ny, nx = arr.shape[:2]
    ii = np.asarray(ii, dtype=float); jj = np.asarray(jj, dtype=float)
    i0 = np.floor(ii).astype(int); j0 = np.floor(jj).astype(int)
    fi = ii - i0; fj = jj - j0
    i0 = np.clip(i0, 0, ny - 1); j0 = np.clip(j0, 0, nx - 1)
    i1 = np.clip(i0 + 1, 0, ny - 1); j1 = np.clip(j0 + 1, 0, nx - 1)
    return np.stack([
        (1 - fi[s]) * (1 - fj[s]) * arr[i0[s], j0[s]]
        + (1 - fi[s]) * fj[s] * arr[i0[s], j1[s]]
        + fi[s] * (1 - fj[s]) * arr[i1[s], j0[s]]
        + fi[s] * fj[s] * arr[i1[s], j1[s]]
        for s in range(len(i0))])


def aperture_band(b_max, lo_frac=0.30, hi_frac=1.0):
    """Trusted wavenumber band of an aperture.

    The band runs from the smallest resolvable phase over the longest baseline,
    lo_frac / b_max, up to the anti-alias ceiling, hi_frac * pi / b_max.
    """
    return lo_frac / b_max, hi_frac * np.pi / b_max


def _log_edge_taper(grid, lo, hi, width):
    """Cosine edge weight on a logarithmic grid.

    The weight is one inside the band [lo, hi] and tapers smoothly to zero over
    `width` octaves at each edge, and is zero outside. This down-weights each
    aperture's anti-alias tail where neighbouring bands overlap and cancels in
    single-aperture regions.
    """
    lg = np.log2(grid)
    llo, lhi = np.log2(lo), np.log2(hi)
    ramp = np.clip(np.minimum((lg - llo) / width, (lhi - lg) / width), 0.0, 1.0)
    t = 0.5 * (1.0 - np.cos(np.pi * ramp))
    return np.where((lg > llo) & (lg < lhi), t, 0.0)


def default_apertures():
    """Default coarse-to-tight ladder of apertures.

    Returns a list of (name, extent_px) tuples spaced by roughly a factor 1.3 so
    that every wavenumber octave is covered by at least two apertures, leaving no
    single-aperture notch at the stitches.
    """
    return [('A0', 28), ('A1', 22), ('A2', 17), ('A3', 13),
            ('A4', 10), ('A5', 8), ('A6', 6), ('A7', 4)]


def auto_apertures(ny, nx, dx, fov_frac=0.55, b_min_px=3.0, ratio=1.4,
                   min_window_px=12):
    """Data-aware coarse-to-tight aperture ladder sized to the field.

    Builds a half-octave ladder of baselines from a fraction of the field of
    view down to a few pixels, so the stitched spectrum spans the swell through
    the wind sea without the caller hand-tuning a pixel ladder. Coarse rungs are
    random windows (well sampled); the tight high-wavenumber rungs are
    deterministic plus (cross) apertures, which stay well conditioned where a
    random window would have too few pixels.

    Args:
        ny, nx (int): Field shape in pixels.
        dx (float): Grid spacing in metres.
        fov_frac (float): Coarsest baseline as a fraction of the (smaller) field
            extent; sets the longest resolved wave (default 0.55).
        b_min_px (float): Tightest baseline in pixels; sets the anti-alias
            ceiling k_max = pi / (b_min_px * dx) (default 3).
        ratio (float): Baseline ratio between successive rungs (default 1.4, a
            half octave, so neighbouring trusted bands overlap).
        min_window_px (int): Rungs with an extent below this many pixels are
            emitted as plus crosses instead of random windows (default 12).

    Returns:
        list: ``(name, extent)`` tuples, where ``extent`` is an integer window
        side in pixels or ``("plus", arm_px)`` for a cross aperture.
    """
    L = min(ny, nx) * dx
    b0 = fov_frac * L
    b_min = b_min_px * dx
    n = max(2, int(np.ceil(np.log(b0 / b_min) / np.log(ratio))) + 1)
    aps = []
    for i, b in enumerate(np.geomspace(b0, b_min, n)):
        ext_px = b / dx
        if ext_px >= min_window_px:
            aps.append((f"A{i}", int(round(ext_px))))
        else:
            aps.append((f"A{i}", ("plus", round(b / 2.0 / dx, 3))))
    return aps


def default_grids(dx, fmin=0.06, fmax=3.5, nfreq=52,
                  kmin=0.08, numin=0.01, numax=2.0, nrad=64):
    """Default log-spaced frequency, wavenumber and inverse phase speed grids.

    The wavenumber grid runs up to the pixel Nyquist pi / dx.
    """
    freqs = np.logspace(np.log10(fmin), np.log10(fmax), nfreq)
    kmax = np.pi / dx
    k_grid = 2.0**np.linspace(np.log2(kmin), np.log2(kmax), nrad)
    nu_grid = 2.0**np.linspace(np.log2(numin), np.log2(numax), nrad)
    return freqs, k_grid, nu_grid
# }}}


# wavelet wavevectors and sign resolution {{{
def cwt_stack(staff_series, freqs, fs, omega0=6.0):
    """Continuous wavelet transform of each staff time series.

    Args:
        staff_series (np.ndarray): Stack of staff time series.
        freqs (np.ndarray): Frequency array.
        fs (float): Sampling frequency.
        omega0 (float): Morlet central frequency. Larger values give a finer
            frequency resolution (df / f is of order 1 / omega0) at the expense
            of a coarser time resolution and therefore noisier statistics.

    Returns:
        np.ndarray: Complex coefficients with shape (n_staff, n_freq, n_time).
    """
    mother = Morlet(omega0)
    return np.stack([cwt(s.astype(float), freqs=freqs, fs=fs, mother=mother).values
                     for s in staff_series])


def solve_wavevectors(W, px, py):
    """Least-squares wavevector from the cross-staff phase differences.

    For each frequency and time, the wavevector is solved from the phase
    differences of every staff pair. Returns a root-mean-square residual rather
    than numpy's least-squares residual; kept separate from
    `Arrays.compute_wavenumbers` so this estimator stays self-contained.

    Args:
        W (np.ndarray): Wavelet coefficients with shape (n_staff, n_freq, n_time).
        px, py (np.ndarray): Staff positions in metres.

    Returns:
        tuple: The two wavevector components, the root-mean-square residual and a
        normalised misfit, each with shape (n_freq, n_time). The misfit is the
        residual divided by the phase-difference magnitude; it is near zero for a
        cleanly resolved plane wave and approaches one when the phase differences
        are inconsistent across pairs (a short wave aliased by a coarse aperture,
        or a noise-dominated sample). It drives the optional reliability gate in
        :meth:`MultiApertureArrays.compute`.

    Notes:
        The misfit catches the dominant spurious shift of energy towards zero
        wavenumber from short waves aliased by an over-large aperture: those
        wrap the phase inconsistently across pairs of different baseline, giving
        a large residual. It makes no dispersion assumption. It does not flag the
        weaker sub-footprint bias, where a wave longer than the aperture gives
        small consistent phase differences and a low misfit; that mode is
        addressed by the de-pistoning hook (`solve_eta`).
    """
    nstaff, nf, T = W.shape
    pairs = [(m, n) for m in range(nstaff) for n in range(m + 1, nstaff)]
    A = np.array([[px[n] - px[m], py[n] - py[m]] for m, n in pairs])

    # initialise the wavevector components, residual and normalised misfit
    kx = np.empty((nf, T))
    ky = np.empty((nf, T))
    resid = np.empty((nf, T))
    misfit = np.empty((nf, T))

    # loop for each frequency and solve the least-squares system over all pairs
    for fi in range(nf):
        dphi = np.empty((len(pairs), T))
        for q, (m, n) in enumerate(pairs):
            dphi[q] = np.angle(W[n, fi] * np.conj(W[m, fi]))
        sol, *_ = np.linalg.lstsq(A, dphi, rcond=None)
        kx[fi], ky[fi] = sol[0], sol[1]
        resid[fi] = np.sqrt(np.mean((A @ sol - dphi)**2, axis=0))
        misfit[fi] = resid[fi] / (np.sqrt(np.mean(dphi**2, axis=0)) + 1e-9)
    return kx, ky, resid, misfit


def lh_direction(We, Wsx, Wsy, flip=False):
    """Longuet-Higgins first-moment mean direction.

    The unambiguous (0 to 360 degrees) mean direction is computed from the
    co-located heave and slope wavelet coefficients, where the slopes are the
    gradient of the elevation. This resolves the 180-degree sign ambiguity of the
    array solution. The result is returned in degrees clockwise from North.
    """
    cross_x = np.mean(We * np.conj(Wsx), axis=(0, 2))
    cross_y = np.mean(We * np.conj(Wsy), axis=(0, 2))
    Cee = np.mean(np.abs(We)**2, axis=(0, 2))
    Cxx = np.mean(np.abs(Wsx)**2, axis=(0, 2))
    Cyy = np.mean(np.abs(Wsy)**2, axis=(0, 2))
    den = np.sqrt(np.clip(Cee * (Cxx + Cyy), 1e-30, None))
    return _math_angle_to_cw_from_N(np.arctan2(np.imag(cross_y) / den,
                                               np.imag(cross_x) / den), flip=flip)
# }}}


# multi-aperture estimator {{{
class MultiApertureArrays:
    """Perform multi-aperture EWDM for dense fields or sparse spatial arrays.

    The class is constructed with one of the factory methods, `from_field` for a
    dense elevation grid or `from_arrays` (`from_numpy`) for a sparse sensor
    array, and the spectra are produced by the `compute` method. See the module
    docstring for the input conventions.

    Arguments:
        field (np.ndarray, optional): Dense elevation field with shape
            (ny, nx, time). Rows are North and columns are East, with a
            West-positive column axis as used by `seed_aperture`.
        dx (float, optional): Grid spacing in metres. Required for a dense field.
        valid (np.ndarray, optional): Boolean footprint of the dense field.
        dataset (xr.Dataset, optional): Sparse sensor array following the
            `Arrays` convention.
        depth (float): Mean water depth in metres.
        fs (float, optional): Sampling frequency in Hz. For a sparse array it is
            inferred from the time coordinate when not given.
        direction_convention (str): Either `cw_from_north` (default) or `math`.
    """

    def __init__(self, *, field=None, dx=None, valid=None,
                 dataset=None, depth=None, fs=None,
                 direction_convention='cw_from_north'):

        # exactly one of a dense field or a sparse dataset must be supplied
        if (field is None) == (dataset is None):
            raise ValueError(
                "Provide exactly one of a dense `field` or a sparse `dataset`. "
                "Use the `from_field` or `from_arrays` factory methods."
            )
        if depth is None:
            raise ValueError("`depth`, the mean water depth in metres, is required.")
        if direction_convention not in ("cw_from_north", "math"):
            raise ValueError(
                "`direction_convention` must be `cw_from_north` or `math`."
            )

        self.field = field
        self.dx = dx
        self.valid = valid
        self.dataset = dataset
        self.depth = float(depth)
        self.fs = fs
        self.direction_convention = direction_convention
        self.mode = "field" if field is not None else "sparse"

        # validate and derive the geometry depending on the input mode
        if self.mode == "field":
            if dx is None:
                raise ValueError(
                    "`dx`, the grid spacing in metres, is required for a dense field."
                )
            self.field = np.asarray(field, dtype=float)
            if self.field.ndim != 3:
                raise ValueError("`field` must have shape (ny, nx, time).")
            ny, nx, _ = self.field.shape
            # frame-fundamental wavenumber, used later by the `auto` nu gate
            self._kfov = 2 * np.pi / (nx * self.dx)
        else:
            self._prepare_sparse(dataset)


    @classmethod
    def from_field(cls, eta, dx, depth, fs, *, valid=None,
                   direction_convention='cw_from_north'):
        """Create an instance from a dense elevation grid.

        Arguments:
            eta: Elevation field with shape (ny, nx, time). Rows are North and
                columns are East, with a West-positive column axis.
            dx: Grid spacing in metres.
            depth: Mean water depth in metres.
            fs: Sampling frequency in Hz.
            valid: Optional (ny, nx) boolean footprint. The elevation may be NaN
                outside it; the variance is taken over the footprint and the
                staffs are drawn from its one-pixel-eroded pixels.
            direction_convention: Either `cw_from_north` (default) or `math`.
        """
        return cls(field=eta, dx=dx, valid=valid, depth=depth, fs=fs,
                   direction_convention=direction_convention)


    @classmethod
    def from_arrays(cls, dataset, depth, fs=None, *,
                    direction_convention='cw_from_north'):
        """Create an instance from a sparse sensor array.

        Arguments:
            dataset: Dataset following the `Arrays` convention, with
                `surface_elevation` (time, element), `position_x` (element) and
                `position_y` (element). Optional `eastward_slope` and
                `northward_slope` (time, element) enable the Longuet-Higgins sign
                fallback.
            depth: Mean water depth in metres.
            fs: Sampling frequency in Hz. If None, it is inferred from the time
                coordinate or the `sampling_rate` attribute.
            direction_convention: Either `cw_from_north` (default) or `math`.
        """
        return cls(dataset=dataset, depth=depth, fs=fs,
                   direction_convention=direction_convention)


    @classmethod
    def from_numpy(cls, time, surface_elevation, position_x, position_y,
                   depth, fs=None, **kwargs):
        """Create an instance of a sparse array from numpy arrays.

        This mirrors `Arrays.from_numpy`. The `surface_elevation` array has shape
        (time, n_elements) and `position_x` and `position_y` have length
        n_elements.
        """
        position_x = np.asarray(position_x)
        position_y = np.asarray(position_y)
        if len(position_x) != len(position_y):
            raise ValueError("Number of elements in `x` and `y` are not consistent.")

        # element numbering is always assumed from 0 to number of elements
        elements = np.arange(len(position_x))

        # create dataset from numpy arrays
        dataset = xr.Dataset(
            data_vars={
                "surface_elevation": (["time", "element"], surface_elevation),
                "position_x": ("element", position_x),
                "position_y": ("element", position_y),
            },
            coords={"time": time, "element": elements},
        )
        return cls.from_arrays(dataset, depth=depth, fs=fs, **kwargs)


    def _prepare_sparse(self, dataset):
        """Validate a sparse dataset and derive its geometry."""

        # the dataset must contain elevation and the element positions
        for var in ("surface_elevation", "position_x", "position_y"):
            if var not in dataset:
                raise ValueError(f"`dataset` is missing the `{var}` variable.")

        # determine the sampling frequency from the argument, the global
        # attribute or the time coordinate, in that order
        if self.fs is None:
            self.fs = dataset.attrs.get("sampling_rate")
        if self.fs is None:
            self.fs = get_sampling_frequency(dataset["time"])
        self.fs = float(self.fs)

        # element positions and the pairwise baselines
        self._px = dataset["position_x"].values.astype(float)
        self._py = dataset["position_y"].values.astype(float)
        b = np.hypot(self._px[:, None] - self._px[None, :],
                     self._py[:, None] - self._py[None, :])

        # the element span and the smallest non-zero baseline set the default
        # grid bounds, since a sparse array has no pixel grid of its own
        span = float(b.max())
        bmin = float(b[b > 0].min()) if np.any(b > 0) else span
        self.dx = bmin / 2.0
        self._kfov = 2 * np.pi / span


    def _sparse_apertures(self):
        """Group elements into nested sub-apertures by baseline.

        The sub-apertures are element subsets within a shrinking radius about the
        array centroid, from coarse to tight, spaced by roughly a factor 1.3.
        Each subset must have at least three non-collinear elements so that the
        two-dimensional least-squares solve is well posed.
        """
        px, py = self._px, self._py
        cx, cy = px.mean(), py.mean()
        dist = np.hypot(px - cx, py - cy)
        rmax = float(dist.max())
        bmin = dist[dist > 0].min() if np.any(dist > 0) else 0.0

        # shrink the radius and keep every subset with enough non-collinear
        # elements, dropping duplicates of the previous subset
        out, r, ai = [], rmax, 0
        while r > 0:
            idx = np.where(dist <= r * 1.0000001)[0]
            if idx.size >= 3:
                pts = np.c_[px[idx] - px[idx].mean(), py[idx] - py[idx].mean()]
                if np.linalg.matrix_rank(pts, tol=1e-9) >= 2:
                    if not out or set(idx) != set(out[-1][1]):
                        out.append((f"A{ai}", idx))
                        ai += 1
            r /= 1.3
            if r < bmin:
                break

        if not out:
            raise ValueError(
                "Sparse array has no aperture with at least three "
                "non-collinear elements."
            )
        return out


    def _build_apertures_field(self, freqs, apertures, n_staff, seed, omega0,
                               solve_eta, lo_frac, lo_frac_broad, hi_frac,
                               flip_dir):
        """Seed virtual staffs and solve a wavevector for each aperture.

        This is the dense-field path. It returns the list of per-aperture
        dictionaries, the Longuet-Higgins reference direction, the footprint
        variance and the time array.
        """
        eta = self.field
        ny, nx, T = eta.shape
        valid = self.valid

        # variance over the whole frame, or over the footprint when one is given.
        # in the latter case the staffs are drawn from a one-pixel-eroded mask so
        # that the gradient central differences stay inside the footprint
        if valid is None:
            var_eta = float(eta.var(axis=2).mean())
            seed_valid = None
        else:
            valid = np.asarray(valid, dtype=bool)
            var_eta = float(np.nanmean(np.where(valid, eta.var(axis=2), np.nan)))
            seed_valid = erode_valid(valid)
        time = np.arange(T) / self.fs

        # loop over the apertures, seeding the staffs and solving the wavevector
        ap = []
        a0_ij = None
        for ai, (name, ext) in enumerate(apertures):
            ii, jj, px, py, bmax = seed_aperture(ny, nx, self.dx, ext, n_staff,
                                                 seed + ai, valid=seed_valid)
            es = _bilinear_stack(eta, ii, jj)
            es = es - es.mean(1, keepdims=True)
            We = cwt_stack(es, freqs, self.fs, omega0)

            # solve the wavevector on the optional de-pistoned field; the power
            # is always taken from the wavelet coefficients of the elevation
            if solve_eta is not None:
                ek = _bilinear_stack(solve_eta, ii, jj)
                Wk = cwt_stack(ek - ek.mean(1, keepdims=True), freqs, self.fs,
                               omega0)
            else:
                Wk = We
            kx, ky, resid, misfit = solve_wavevectors(Wk, px, py)

            # the coarsest aperture uses the broad lower fraction to reach lower k
            klo, khi = aperture_band(
                bmax, lo_frac=(lo_frac_broad if ai == 0 else lo_frac),
                hi_frac=hi_frac)
            ap.append(dict(name=name, bmax=bmax, klo=klo, khi=khi,
                           kmag=np.hypot(kx, ky), resid=resid, misfit=misfit, We=We,
                           dir=_math_angle_to_cw_from_N(np.arctan2(ky, kx),
                                                        flip=flip_dir),
                           P=(np.abs(We)**2).mean(0)))
            if ai == 0:
                a0_ij = (ii, jj)

        # Longuet-Higgins first moment from the coarsest aperture's staffs, using
        # the spatial gradient as the local slope. `seed_aperture` uses an
        # East/North frame (East = (cxp - j) dx, North = (cyp - i) dx) opposite to
        # the pixel row/column indices, so the slopes must match:
        # d/dEast = -d/d(col), d/dNorth = -d/d(row). Otherwise the heave-slope
        # cross spectrum is sign-flipped and the resolved hemisphere is 180
        # degrees from the measured wavevector.
        gy, gx = np.gradient(eta, self.dx, axis=(0, 1))
        gx, gy = -gx, -gy
        ii, jj = a0_ij
        sxs = np.stack([gx[i, j, :] for i, j in zip(ii, jj)])
        sys = np.stack([gy[i, j, :] for i, j in zip(ii, jj)])
        Wsx = cwt_stack(sxs - sxs.mean(1, keepdims=True), freqs, self.fs, omega0)
        Wsy = cwt_stack(sys - sys.mean(1, keepdims=True), freqs, self.fs, omega0)
        lh = lh_direction(ap[0]['We'], Wsx, Wsy, flip=flip_dir)
        return ap, lh, var_eta, time


    def _build_apertures_sparse(self, freqs, omega0, solve_eta, lo_frac,
                                lo_frac_broad, hi_frac, flip_dir):
        """Group elements into nested sub-apertures and solve a wavevector each.

        This is the sparse-array path. It returns the same per-aperture
        dictionaries, Longuet-Higgins reference, footprint variance and time
        array as the dense-field path.
        """
        ds = self.dataset
        elev = ds["surface_elevation"].transpose("element", "time").values
        elev = elev.astype(float)
        nelem, T = elev.shape
        time = np.arange(T) / self.fs

        # mean variance over the elements
        var_eta = float(np.nanmean(np.nanvar(elev, axis=1)))

        # optional de-pistoned solve field, in the same element ordering
        if solve_eta is not None:
            solve = np.asarray(solve_eta, dtype=float)
            if solve.shape != elev.shape:
                raise ValueError(
                    "`solve_eta` must have shape (n_element, time)."
                )
        else:
            solve = None

        # loop over the baseline-grouped sub-apertures and solve the wavevector
        groups = self._sparse_apertures()
        ap, a0_idx = [], None
        for ai, (name, idx) in enumerate(groups):
            px, py = self._px[idx], self._py[idx]
            bmax = float(np.hypot(px[:, None] - px[None, :],
                                  py[:, None] - py[None, :]).max())
            es = elev[idx] - elev[idx].mean(1, keepdims=True)
            We = cwt_stack(es, freqs, self.fs, omega0)
            if solve is not None:
                ek = solve[idx] - solve[idx].mean(1, keepdims=True)
                Wk = cwt_stack(ek, freqs, self.fs, omega0)
            else:
                Wk = We
            kx, ky, resid, misfit = solve_wavevectors(Wk, px, py)
            klo, khi = aperture_band(
                bmax, lo_frac=(lo_frac_broad if ai == 0 else lo_frac),
                hi_frac=hi_frac)
            ap.append(dict(name=name, bmax=bmax, klo=klo, khi=khi,
                           kmag=np.hypot(kx, ky), resid=resid, misfit=misfit, We=We,
                           dir=_math_angle_to_cw_from_N(np.arctan2(ky, kx),
                                                        flip=flip_dir),
                           P=(np.abs(We)**2).mean(0)))
            if ai == 0:
                a0_idx = idx

        # the Longuet-Higgins sign fallback is only available when co-located
        # slopes are present in the dataset; otherwise the caller must resolve
        # the sign through the `sign_anchor` hook
        if "eastward_slope" in ds and "northward_slope" in ds:
            sx = ds["eastward_slope"].transpose("element", "time").values[a0_idx]
            sy = ds["northward_slope"].transpose("element", "time").values[a0_idx]
            Wsx = cwt_stack(sx - sx.mean(1, keepdims=True), freqs, self.fs, omega0)
            Wsy = cwt_stack(sy - sy.mean(1, keepdims=True), freqs, self.fs, omega0)
            lh = lh_direction(ap[0]['We'], Wsx, Wsy, flip=flip_dir)
        else:
            logger.warning(
                "Sparse array has no `eastward_slope` or `northward_slope`, so "
                "the Longuet-Higgins sign fallback is unavailable. Pass a "
                "`sign_anchor` to resolve the 180-degree ambiguity."
            )
            lh = None
        return ap, lh, var_eta, time


    def _run_multiaperture(self, ap, lh, var_eta, time, freqs, k_grid, nu_grid,
                           dd, kappa, rel_bandwidth, radial_bandwidth_mode,
                           power_weighted, nu_lo_broad, nu_f_lim, nu_k_lim,
                           stitch_taper, antialias_gate, antialias_mult,
                           sign_anchor, sign_anchor_rmin, sign_coh_min,
                           reliability_gate):
        """Resolve the sign, gate and stitch the per-aperture solutions.

        This is the common estimator shared by both input modes. It works only on
        the assembled per-aperture dictionaries and returns a dictionary of the
        composite spectra and diagnostics.
        """
        depth = self.depth
        coords = {'frequency': freqs, 'time': time}

        # sign reference theta(f): per-frequency axis from the energy-weighted
        # measured wavevector (smoothed double-angle mean), unwrapped across
        # frequency into a signed direction, with the global sign set from the
        # anchor or, failing that, the energy-weighted Longuet-Higgins hemisphere
        P0 = ap[0]['P']
        dir0 = ap[0]['dir']
        wf = P0.mean(1)
        ca = (P0 * np.cos(np.radians(2 * dir0))).sum(1)
        sa = (P0 * np.sin(np.radians(2 * dir0))).sum(1)
        R2 = np.hypot(ca, sa) / np.maximum(P0.sum(1), 1e-30)   # axis coherence
        sm = np.array([0.25, 0.5, 0.25])
        axis = np.degrees(np.arctan2(np.convolve(sa, sm, mode='same'),
                                     np.convolve(ca, sm, mode='same'))) / 2.0

        # pick whichever of the two axis branches is closest to the running value
        def _branch(a, prev):
            c1, c2 = _wrap180(a), _wrap180(a + 180.0)
            return c1 if abs(_wrap180(c1 - prev)) <= abs(_wrap180(c2 - prev)) else c2

        # seed the unwrap at the spectral peak and propagate outward. above the
        # peak the axis is always unwrapped; below the peak only the coherent
        # frequencies update the running direction, the incoherent ones hold it
        pk = int(np.argmax(wf))
        theta = np.empty(len(freqs))
        cur = axis[pk]
        theta[pk] = cur
        for i in range(pk + 1, len(freqs)):
            cur = _branch(axis[i], cur)
            theta[i] = cur
        cur = theta[pk]
        for i in range(pk - 1, -1, -1):
            if R2[i] >= sign_coh_min:
                cur = _branch(axis[i], cur)
            theta[i] = cur

        # energy-weighted circular mean of a set of directions
        def _ewmean(deg, w):
            return np.degrees(np.arctan2((w * np.sin(np.radians(deg))).sum(),
                                         (w * np.cos(np.radians(deg))).sum()))

        # set the global sign from the external anchor when it is reliable, then
        # from the Longuet-Higgins direction, and otherwise leave the axis as is
        if sign_anchor is not None and (
                np.isfinite(np.asarray(sign_anchor[1], float))
                & (np.asarray(sign_anchor[2], float) >= sign_anchor_rmin)).any():
            af, ad, aR = (np.asarray(x, float) for x in sign_anchor)
            ok = np.isfinite(ad) & (aR >= sign_anchor_rmin)
            jc = np.argmin(np.abs(af[ok, None] - freqs[None, :]), axis=1)
            w = aR[ok] * wf[jc]
            flip = np.cos(np.radians(_ewmean(theta[jc], w) - _ewmean(ad[ok], w))) < 0
        elif lh is not None:
            flip = np.cos(np.radians(_ewmean(theta, wf) - _ewmean(lh, wf))) < 0
        else:
            flip = False
        ref = _wrap180(theta + 180.0) if flip else theta

        # fold every aperture's direction into the signed hemisphere
        def _fold(dir_cwN):
            mis = np.cos(np.radians(dir_cwN - ref[:, None])) < 0
            return np.where(mis, _wrap180(dir_cwN + 180.0), dir_cwN)
        for d in ap:
            d['dirf'] = _fold(d['dir'])

        # one variance calibration from the full-frame power so that the integral
        # of S(f) over frequency equals the elevation variance
        cal = var_eta / _trapezoid(ap[0]['P'].mean(1), freqs)

        # anti-alias gate threshold (see the `compute` docstring)
        kdisp = k_dispersion(freqs, depth)
        gate_k = antialias_mult * float(k_dispersion(freqs[pk], depth))

        # the Q(nu) deposit is restricted to a trusted scale window from the array
        # baselines, with a high cut at the frame-fundamental wavenumber, and the
        # matching dispersion band in frequency. Only Q(nu) is gated this way
        bmx = [d['bmax'] for d in ap]
        if nu_k_lim == 'auto':
            nu_k_lim = (2 * np.pi / (50.0 * max(bmx)), self._kfov)
        if nu_f_lim == 'auto':
            _fk = lambda kk: float(np.sqrt(
                GRAV * kk * np.tanh(np.clip(kk * depth, 1e-9, 50)))
                / (2 * np.pi))
            nu_f_lim = ((_fk(nu_k_lim[0]), _fk(nu_k_lim[1]))
                        if nu_k_lim is not None else None)

        # accumulate the stitched wavenumber and inverse phase speed spectra
        bins_dir = np.arange(-180.0, 180.0, dd)
        Fk = np.zeros(len(k_grid)); ck = np.zeros(len(k_grid))
        Qn = np.zeros(len(nu_grid)); cn = np.zeros(len(nu_grid))
        Fkd = np.zeros((len(k_grid), len(bins_dir)))
        Qnd = np.zeros((len(nu_grid), len(bins_dir)))
        Ffd = None        # stitched frequency-direction deposit (filled in loop)
        ap_ok_omni = []   # per-aperture omnidirectional F(k), for diagnostics
        for d in ap:
            # anti-alias frequency mask: send the frequencies that would alias off
            # the k and nu grids so they never deposit, while keeping the power
            # intact so the power-weighting is undisturbed
            fok = (kdisp <= d['khi']) if (antialias_gate and d['khi'] < gate_k) \
                else np.ones(len(freqs), bool)
            # dispersion-independent reliability mask: drop (freq, time) samples
            # whose least-squares phase misfit is high (a short wave aliased by
            # this aperture, or noise) so their spurious low wavenumber never
            # deposits. The power is kept intact, so the power-weighting is
            # undisturbed. On by default (reliability_gate=0.6); None disables it.
            relok = (d['misfit'] <= reliability_gate) \
                if reliability_gate is not None else np.ones_like(d['kmag'], bool)
            power = xr.DataArray(d['P'] * cal, dims=['frequency', 'time'],
                                 coords=coords)
            thd = xr.DataArray(d['dirf'], dims=['frequency', 'time'],
                               coords=coords)

            # deposit onto the wavenumber grid
            k_vals = np.where(fok[:, None] & relok, d['kmag'], 1e12)
            kk = xr.DataArray(k_vals, dims=['frequency', 'time'], coords=coords)
            ok = estimate_radial_distribution(power, thd, kk, 'wavenumber',
                                              k_grid, dd, kappa,
                                              bandwidth=rel_bandwidth,
                                              bandwidth_mode=radial_bandwidth_mode,
                                              power_weighted=power_weighted)

            # deposit onto the inverse phase speed grid, gated to the trusted
            # (f, k) window. out-of-window samples are sent off the nu grid, again
            # keeping the power intact so F(k) and S(f) are unaffected
            nu_vals = d['kmag'] / (2 * np.pi * freqs[:, None])
            if nu_f_lim is not None or nu_k_lim is not None:
                keep = np.ones(nu_vals.shape, bool)
                if nu_f_lim is not None:
                    keep &= ((freqs >= nu_f_lim[0]) & (freqs <= nu_f_lim[1]))[:, None]
                if nu_k_lim is not None:
                    keep &= (d['kmag'] >= nu_k_lim[0]) & (d['kmag'] <= nu_k_lim[1])
                nu_vals = np.where(keep, nu_vals, 1e6)
            nu_vals = np.where(fok[:, None] & relok, nu_vals, 1e12)
            nu = xr.DataArray(nu_vals, dims=['frequency', 'time'], coords=coords)
            on = estimate_radial_distribution(power, thd, nu, 'nu', nu_grid,
                                              dd, kappa, bandwidth=rel_bandwidth,
                                              bandwidth_mode=radial_bandwidth_mode,
                                              power_weighted=power_weighted)

            ok_omni = ok['wavenumber_spectrum'].values
            ok_dir = ok['directional_spectrum'].values
            on_omni = on['nu_spectrum'].values
            on_dir = on['directional_spectrum'].values

            # the nu band follows the aperture's k-band along the dispersion
            # relation; the broadest aperture is extended down to `nu_lo_broad` to
            # show the measured inverse phase speed tail
            ink = (k_grid >= d['klo']) & (k_grid <= d['khi'])
            nlo, nhi = _nu_of_k(d['khi'], depth), _nu_of_k(d['klo'], depth)
            lo = min(nlo, nhi)
            if d is ap[0]:
                lo = min(lo, nu_lo_broad)

            # stitch weights: a cosine edge taper in log-k and log-nu when one is
            # requested, otherwise the hard in-band indicator. The weight cancels
            # in the single-aperture regions
            wk = _log_edge_taper(k_grid, d['klo'], d['khi'], stitch_taper) \
                if stitch_taper else ink.astype(float)
            wn = _log_edge_taper(nu_grid, lo, max(nlo, nhi), stitch_taper) \
                if stitch_taper else \
                ((nu_grid >= lo) & (nu_grid <= max(nlo, nhi))).astype(float)

            ap_ok_omni.append((ok_omni.copy(), ink.copy()))
            Fk += ok_omni * wk; ck += wk
            Fkd += ok_dir * wk[:, None]
            Qn += on_omni * wn; cn += wn
            Qnd += on_dir * wn[:, None]

            # stitched frequency-direction deposit (dispersion-free): this
            # aperture contributes E(f, theta) from its reliable, non-aliased
            # samples only (the `relok` phase-misfit gate). summed over apertures,
            # each frequency takes its direction from whichever apertures resolve
            # it, so the coarsest aperture's aliased high-frequency direction does
            # not flip the estimate.
            pgd = xr.DataArray(d['P'] * cal * relok, dims=['frequency', 'time'],
                               coords=coords)
            ofd = estimate_directional_distribution(pgd, thd, dd, kappa)
            Ffd = ofd['directional_spectrum'].values if Ffd is None \
                else Ffd + ofd['directional_spectrum'].values

        # normalise the stitched spectra by the accumulated weights
        Fk = np.where(ck > 0, Fk / np.maximum(ck, 1e-30), np.nan)
        Fkd /= np.maximum(ck, 1e-30)[:, None]
        Qn = np.where(cn > 0, Qn / np.maximum(cn, 1e-30), np.nan)
        Qnd /= np.maximum(cn, 1e-30)[:, None]

        # frequency spectrum S(f) from the coarsest aperture: its power is valid
        # at every frequency (only its direction aliases above the baseline
        # limit). F(f, theta) is the reliability-stitched deposit, renormalised
        # so it integrates over direction to S(f).
        power0 = xr.DataArray(ap[0]['P'] * cal, dims=['frequency', 'time'],
                              coords=coords)
        theta0 = xr.DataArray(ap[0]['dirf'], dims=['frequency', 'time'],
                              coords=coords)
        of = estimate_directional_distribution(power0, theta0, dd, kappa)
        Sf = of['frequency_spectrum'].values
        thg = of['direction'].values
        row = Ffd.sum(1) * np.radians(dd)
        Fft = np.divide(Ffd, np.maximum(row, 1e-30)[:, None]) * Sf[:, None]
        thbar, sigma = circ_stats(Fft, thg)

        return dict(freqs=freqs, theta=thg, k=k_grid, nu=nu_grid,
                    var_eta=var_eta, Sf=Sf, Fft=Fft, Fk=Fk, Qn=Qn, Fkd=Fkd,
                    Qnd=Qnd, thbar=thbar, sigma=sigma,
                    lh=(lh if lh is not None else np.full(len(freqs), np.nan)),
                    sign_ref=ref, ap_names=[d['name'] for d in ap],
                    ap_bands=[(d['klo'], d['khi']) for d in ap],
                    ap_bmax=[d['bmax'] for d in ap],
                    ap_ok_omni=ap_ok_omni, ck=ck)


    def compute(self, *, freqs=None, k_grid=None, nu_grid=None,
                apertures=None, n_staff=16, seed=20,
                dd=4.0, kappa=36.0, omega0=12.0,
                lo_frac=1.0, lo_frac_broad=0.05, hi_frac=1.0, flip_dir=False,
                rel_bandwidth=0.08, radial_bandwidth_mode='relative',
                power_weighted=True, nu_lo_broad=0.04,
                nu_f_lim='auto', nu_k_lim='auto',
                stitch_taper=0.7, solve_eta=None, depiston=False,
                antialias_gate=False, antialias_mult=3.0,
                reliability_gate=0.6,
                sign_anchor=None, sign_anchor_rmin=0.15, sign_coh_min=0.6,
                return_apertures=False) -> xr.Dataset:
        """Perform the multi-aperture computation.

        The result is a dataset with the omnidirectional spectra S(f), F(k) and
        Q(nu), the polar directional spectra F(f, theta), Psi(k, theta) and
        Q(nu, theta), and the mean direction, directional spread, Longuet-Higgins
        reference and unwrapped sign reference. The direction is in degrees
        clockwise from North unless the instance was built with
        `direction_convention="math"`.

        Args:
            freqs, k_grid, nu_grid (np.ndarray, optional): Log-spaced radial
                grids. When None they default to `default_grids`; the sparse path
                derives the pixel-Nyquist ceiling from the smallest baseline, so
                explicit grids are recommended for a sparse array.
            apertures (list, optional): List of (name, extent) tuples for the
                dense path, where the extent is a window side in pixels (a scalar
                square or a (rows, cols) tuple) or ("plus", arm_px) for a cross
                aperture. Defaults to `auto_apertures`, a data-aware ladder sized
                to the field of view (coarse windows down to tight crosses), so
                the stitched spectrum spans the full wavenumber range without
                hand-tuning. Ignored for a sparse array, whose sub-apertures are
                grouped by baseline.
            n_staff (int): Number of virtual staffs per aperture (dense path).
            seed (int): Seed for the random staff placement (dense path).
            dd (float): Directional resolution in degrees (default 4 degrees).
            kappa (float): Smoothness parameter of the von Mises kernel.
            omega0 (float): Morlet central frequency.
            lo_frac, lo_frac_broad, hi_frac (float): Trusted-band fractions. The
                coarsest aperture uses `lo_frac_broad` to reach lower wavenumbers.
            flip_dir (bool): Globally flip the measured direction convention.
            rel_bandwidth (float): Radial kernel width forwarded to
                `estimate_radial_distribution`.
            radial_bandwidth_mode (str): Radial kernel mode, one of `relative`
                (a log-space fractional width), `histogram` or `absolute`.
            power_weighted (bool): Bin the variance rather than the sample
                occurrence into F(k) and Q(nu).
            nu_lo_broad (float): Extend the coarsest aperture's Q(nu) band down to
                this inverse phase speed to show the measured tail.
            nu_f_lim, nu_k_lim (str, tuple or None): Gates applied only to Q(nu).
                `auto` derives them from the array baselines, with the high cut at
                the frame-fundamental wavenumber. A tuple overrides them and None
                disables a gate.
            stitch_taper (float): Cosine edge-taper width in octaves blending the
                overlapping aperture bands. Zero or None gives a hard in-band
                indicator.
            solve_eta (np.ndarray, optional): Field or array, the same shape as
                the elevation, used only for the wavevector solve. The power and
                variance still come from the elevation. Pass a de-pistoned field
                to keep the uniform long wave out of the cross-staff phases.
            depiston (bool): Built-in de-piston. When True (and `solve_eta` is
                not given), the per-frame spatial mean is removed from a copy of
                the field used only for the wavevector solve, suppressing the
                residual low-wavenumber bias from waves longer than the footprint
                (which the reliability gate cannot catch, since their phase
                differences are small but consistent). Default False; enable for
                field-of-view-limited data where the dominant wave approaches or
                exceeds the footprint.
            antialias_gate (bool): When True, an aperture whose ceiling is below
                `antialias_mult` times the spectral peak wavenumber deposits onto
                the k and nu grids only the frequencies whose dispersion
                wavenumber is below the ceiling. S(f) and F(f, theta) are
                unaffected. Enable this when the widest baseline exceeds the
                dominant wavelength.
            antialias_mult (float): Multiple of the peak wavenumber above which an
                aperture is considered safe from aliasing.
            reliability_gate (float or None): The default low-k guard. Drops
                every (frequency, time) sample whose normalised least-squares
                phase misfit (see :func:`solve_wavevectors`) exceeds this
                threshold before it is binned, so short waves aliased by an
                over-large aperture do not deposit their spurious low wavenumber.
                The misfit separates cleanly-resolved samples (-> 0) from aliased
                or noisy ones (-> 1), so the gate removes only the unreliable
                samples. It makes no dispersion assumption, so it also handles
                bound waves or a Doppler-shifted sea. Default ``0.6``; set to None
                to disable. S(f)/F(f, theta) are unaffected.
            sign_anchor (tuple, optional): External (frequency, direction, R)
                estimate resolving the 180-degree ambiguity, for example from a
                3D FFT. It overrides the Longuet-Higgins fallback.
            sign_anchor_rmin (float): Reliability floor for the anchor.
            sign_coh_min (float): Coherence gate for the below-peak sign unwrap.
            return_apertures (bool): Attach the per-aperture diagnostics on an
                `aperture` dimension.

        Returns:
            xr.Dataset: Dataset containing the omnidirectional and directional
            spectra and the directional diagnostics.
        """

        # default radial grids when none are supplied
        if freqs is None or k_grid is None or nu_grid is None:
            df, dk, dn = default_grids(self.dx)
            freqs = df if freqs is None else freqs
            k_grid = dk if k_grid is None else k_grid
            nu_grid = dn if nu_grid is None else nu_grid
        freqs = np.asarray(freqs, float)
        k_grid = np.asarray(k_grid, float)
        nu_grid = np.asarray(nu_grid, float)

        # built-in de-piston: when requested and no explicit solve field is
        # given, remove the per-frame spatial mean (the uniform "piston" mode of
        # waves longer than the footprint) from a copy used only for the
        # wavevector solve. targets the residual low-wavenumber bias the
        # reliability gate cannot reach (a sub-footprint wave gives small but
        # consistent phase differences with |k| biased towards zero).
        if depiston and solve_eta is None:
            if self.mode == "field":
                f = self.field
                if self.valid is not None:
                    masked = np.where(self.valid[:, :, None], f, np.nan)
                    solve_eta = f - np.nanmean(masked, axis=(0, 1), keepdims=True)
                else:
                    solve_eta = f - f.mean(axis=(0, 1), keepdims=True)
            else:
                elev = self.dataset["surface_elevation"].transpose(
                    "element", "time").values.astype(float)
                solve_eta = elev - elev.mean(0, keepdims=True)

        # build the per-aperture solutions for the chosen input mode
        if self.mode == "field":
            if apertures is None:
                ny, nx = self.field.shape[:2]
                apertures = auto_apertures(ny, nx, self.dx)
            ap, lh, var_eta, time = self._build_apertures_field(
                freqs, apertures, n_staff, seed, omega0, solve_eta,
                lo_frac, lo_frac_broad, hi_frac, flip_dir)
        else:
            if apertures is not None:
                logger.warning(
                    "`apertures` is ignored for a sparse array; the "
                    "sub-apertures are grouped by baseline."
                )
            ap, lh, var_eta, time = self._build_apertures_sparse(
                freqs, omega0, solve_eta, lo_frac, lo_frac_broad, hi_frac,
                flip_dir)
            if len(ap) == 1:
                logger.warning(
                    "Sparse array yielded a single aperture, so the result "
                    "degenerates to a single-aperture solve."
                )

        # run the common estimator and wrap the result into a dataset
        result = self._run_multiaperture(
            ap, lh, var_eta, time, freqs, k_grid, nu_grid, dd, kappa,
            rel_bandwidth, radial_bandwidth_mode, power_weighted, nu_lo_broad,
            nu_f_lim, nu_k_lim, stitch_taper, antialias_gate, antialias_mult,
            sign_anchor, sign_anchor_rmin, sign_coh_min, reliability_gate)

        return self._to_dataset(
            result, return_apertures=return_apertures, n_staff=n_staff,
            seed=seed, omega0=omega0, kappa=kappa, dd=dd, lo_frac=lo_frac,
            lo_frac_broad=lo_frac_broad, hi_frac=hi_frac,
            stitch_taper=stitch_taper, antialias_gate=antialias_gate,
            antialias_mult=antialias_mult, power_weighted=power_weighted,
            radial_bandwidth_mode=radial_bandwidth_mode,
            rel_bandwidth=rel_bandwidth)


    def _to_dataset(self, r, *, return_apertures, **attrs):
        """Wrap the estimator result dictionary into a labelled dataset."""

        # the engine works in degrees clockwise from North. when the math angle
        # convention is requested the direction coordinate and the per-frequency
        # direction diagnostics are remapped, and the dataset is later sorted
        direction = r['theta'].astype(float)
        freq_dir = (r['thbar'], r['sigma'], r['lh'], r['sign_ref'])
        if self.direction_convention == 'math':
            direction = _cw_from_N_to_math_angle(direction)
            freq_dir = tuple(_cw_from_N_to_math_angle(a) for a in freq_dir)
        thbar, sigma, lh, sign_ref = freq_dir

        # assemble the data variables and coordinates
        data_vars = {
            "frequency_spectrum": (["frequency"], r['Sf']),
            "wavenumber_spectrum": (["wavenumber"], r['Fk']),
            "nu_spectrum": (["nu"], r['Qn']),
            "directional_spectrum_f": (["frequency", "direction"], r['Fft']),
            "directional_spectrum_k": (["wavenumber", "direction"], r['Fkd']),
            "directional_spectrum_nu": (["nu", "direction"], r['Qnd']),
            "mean_direction": (["frequency"], thbar),
            "directional_spread": (["frequency"], sigma),
            "lh_direction": (["frequency"], lh),
            "sign_reference": (["frequency"], sign_ref),
            "var_eta": ((), r['var_eta']),
        }
        coords = {
            "frequency": r['freqs'],
            "wavenumber": r['k'],
            "nu": r['nu'],
            "direction": direction,
        }
        ds = xr.Dataset(data_vars=data_vars, coords=coords)
        if self.direction_convention == 'math':
            ds = ds.sortby("direction")

        # attach the variable metadata
        for coord in ("frequency", "wavenumber", "nu", "direction"):
            if coord in VARIABLE_NAMES:
                ds[coord].attrs = VARIABLE_NAMES[coord]
        for var in ds.data_vars:
            if var in VARIABLE_NAMES:
                ds[var].attrs = VARIABLE_NAMES[var]

        # optional per-aperture diagnostics on an aperture dimension
        if return_apertures:
            nap = len(r['ap_names'])
            ap_fk = np.stack([o for o, _ in r['ap_ok_omni']])
            ds = ds.assign_coords(aperture=np.arange(nap))
            ds["aperture_name"] = ("aperture", np.array(r['ap_names'], dtype="U8"))
            ds["aperture_bmax"] = ("aperture", np.asarray(r['ap_bmax']))
            ds["aperture_klo"] = ("aperture",
                                  np.array([b[0] for b in r['ap_bands']]))
            ds["aperture_khi"] = ("aperture",
                                  np.array([b[1] for b in r['ap_bands']]))
            ds["aperture_Fk"] = (["aperture", "wavenumber"], ap_fk)
            ds["stitch_weight_k"] = (["wavenumber"], r['ck'])

        # global attributes, casting booleans and None to netcdf-safe scalars
        safe = {"method": "multiaperture EWDM",
                "mode": self.mode,
                "depth": self.depth,
                "dx": float(self.dx),
                "fs": float(self.fs),
                "direction_convention": self.direction_convention}
        for key, value in attrs.items():
            if isinstance(value, bool):
                value = int(value)
            elif value is None:
                value = "None"
            safe[key] = value
        ds.attrs = safe
        return ds
# }}}

"""Helper functions for the multi-aperture gallery example
(:ref:`sphx_glr_auto_gallery_plot_multiaperture.py`).

Gathers the pieces that would otherwise clutter the example: a power-normalised
circular Tukey window and the per-method spectral estimators (a sparse pentagon
via `ewdm.Arrays`, the direct 3-D FFT spectra and a spreading normaliser) used to
build the comparison figures.
"""
import numpy as np
import xarray as xr
from scipy.signal.windows import tukey
from scipy.ndimage import convolve1d
import ewdm
from ewdm.multiaperture import k_dispersion
from ewdm.parameters import GRAV


# -- power-normalised circular Tukey window (radial space + temporal) ----------
def _radial_tukey(s1, s2, tw):
    md = min(s1, s2)
    y, x = np.mgrid[1:s1 + 1, 1:s2 + 1].astype(float)
    x -= x.mean(); y -= y.mean()
    r = np.sqrt(x ** 2 + y ** 2) / (md / 2)
    cp, rs = tw * 2, 1 - tw
    w = (np.cos(2 * np.pi / cp * (r - rs)) + 1) / 2
    w[r < (1 - tw)] = 1.0
    w[r > 1] = 0.0
    return w


def circular_tukey(a, taper_width=0.2, normalization="power", temporal_alpha=0.0):
    """Apply a radial spatial Tukey (and optional temporal Tukey) to a field
    ``a`` of shape (ny, nx) or (ny, nx, time), power-normalised by default."""
    a = np.asarray(a, float)
    s1, s2, s3 = a.shape if a.ndim == 3 else (a.shape[0], a.shape[1], 1)
    w2d = _radial_tukey(s1, s2, taper_width)
    w_t = tukey(s3, temporal_alpha) if (temporal_alpha > 0 and s3 > 1) else np.ones(s3)
    if normalization == "power":
        C = 1.0 / np.sqrt(np.mean(w2d ** 2) * np.mean(w_t ** 2))
    elif normalization == "spectral":
        C = np.sqrt(s1 * s2) / np.linalg.norm(w2d, 2)
    else:
        C = 1.0
    if a.ndim == 3:
        return C * w2d[:, :, None] * w_t[None, None, :] * a
    return C * w2d * a


# -- per-method spectral estimators -------------------------------------------
def pentagon_spectra(field, X, Y, dx, fs, t, indices):
    """Run `ewdm.Arrays` on a set of pixels (one fixed array). Returns a dict
    with the frequency and wavenumber spectra, the variance-calibration factor
    and the anti-alias ceiling k_max = pi / b_max."""
    px = np.array([X[i, j] for i, j in indices])
    py = np.array([Y[i, j] for i, j in indices])
    eta = np.array([field[i, j, :] for i, j in indices]).T
    bmax = float(np.hypot(px[:, None] - px[None, :], py[:, None] - py[None, :]).max())
    vpix = float(np.var(eta, axis=0).mean())
    arr = ewdm.Arrays.from_numpy(time=t, surface_elevation=eta,
                                 position_x=px, position_y=py, fs=fs)
    Sf = arr.compute(cross_wavelet=True, kappa=36, coordinate="frequency")
    Fk = arr.compute(cross_wavelet=True, kappa=36, coordinate="wavenumber")
    fac = vpix / np.trapezoid(Sf["frequency_spectrum"].values, Sf["frequency"].values)
    return dict(Sf=Sf, Fk=Fk, fac=fac, kmax=np.pi / bmax)


def _log_edges(c):
    e = np.empty(len(c) + 1)
    e[1:-1] = np.sqrt(c[:-1] * c[1:])
    e[0] = c[0] ** 2 / e[1]; e[-1] = c[-1] ** 2 / e[-2]
    return e


def direct_spectra(field, dx, fs, var_eta, freqs, k_grid, pad=2, ddir=4.0):
    """Windowed 3-D FFT of the whole field. Returns the omnidirectional S(f) and
    F(k) and the spoke-corrected directional D(f, theta), D(k, theta).

    The omnidirectional S(f) and F(k) come from a sharply windowed, unpadded FFT
    (spatial Tukey 0.2, temporal 0.25) that keeps the full spatial resolution, so
    the swell/wind-sea saddle in F(k) is resolved rather than smeared into a flat
    shelf. The 2-D directional fields, which are noisier and benefit from
    smoothing, use a wider taper (0.5) and a `pad`-fold zero-pad."""
    ny, nx, T = field.shape
    eta0 = field - field.mean(2, keepdims=True)
    fN = np.fft.rfftfreq(T, 1 / fs)
    fed = _log_edges(freqs)
    fin = np.digitize(fN, fed) - 1
    fok = (fin >= 0) & (fin < len(freqs))

    # omnidirectional S(f), F(k) from a sharp, unpadded FFT (full resolution)
    A = np.fft.rfftn(circular_tukey(eta0, 0.2, "power", 0.25), axes=(0, 1, 2))
    P = np.abs(A) ** 2; del A
    nf = P.shape[2]; P[:, :, 1:nf - 1] *= 2.0; P *= var_eta / P.sum()
    kmag = np.hypot(*np.meshgrid(2 * np.pi * np.fft.fftfreq(nx, dx),
                                 2 * np.pi * np.fft.fftfreq(ny, dx))).ravel()
    Sf_d = np.zeros(len(freqs))
    np.add.at(Sf_d, fin[fok], P.sum((0, 1))[fok])
    Sf_d /= np.diff(fed)
    dkf = 2 * np.pi / (nx * dx)
    klin = np.arange(dkf, np.pi / dx, dkf)
    Fk_d = np.histogram(kmag, bins=np.r_[klin - dkf / 2, klin[-1] + dkf / 2],
                        weights=P.sum(2).ravel())[0] / dkf
    del P

    # directional spectra from a wider-tapered, zero-padded FFT (smoother 2-D
    # fields), with a spoke-correction (divide by Cartesian cell count per bin)
    # and angular + radial smoothing of the hard-binned result
    A = np.fft.rfftn(circular_tukey(eta0, 0.5, "power", 0.5),
                     s=(pad * ny, pad * nx, T), axes=(0, 1, 2))
    P = np.abs(A) ** 2; del A
    P[:, :, 1:nf - 1] *= 2.0; P *= var_eta / P.sum()
    kxv = 2 * np.pi * np.fft.fftfreq(pad * nx, dx)
    kyv = 2 * np.pi * np.fft.fftfreq(pad * ny, dx)
    KX, KY = np.meshgrid(kxv, kyv)
    kmag = np.hypot(KX, KY).ravel()
    ndir = int(360 // ddir); dbins = np.arange(-180, 180, ddir)
    didx = (np.rint((np.degrees(np.arctan2(KY, KX)) + 180) / ddir).astype(int) % ndir).ravel()
    ked = _log_edges(k_grid)
    kri = np.digitize(kmag, ked) - 1
    kok = (kri >= 0) & (kri < len(k_grid))
    cnt_th = np.bincount(didx, minlength=ndir).astype(float)
    cnt_kth = np.zeros((len(k_grid), ndir))
    np.add.at(cnt_kth, (kri[kok], didx[kok]), 1.0)
    Fkd = np.zeros((len(k_grid), ndir))
    np.add.at(Fkd, (kri[kok], didx[kok]), P.sum(2).ravel()[kok])
    Fft = np.zeros((len(freqs), ndir))
    for fi in range(1, nf):
        ft = fin[fi]
        if 0 <= ft < len(freqs):
            Fft[ft] += np.bincount(didx, weights=P[:, :, fi].ravel(), minlength=ndir)
    del P
    ang = np.hanning(7); ang /= ang.sum()
    rad = np.hanning(7); rad /= rad.sum()
    Fft = convolve1d(Fft / np.maximum(cnt_th, 1.0)[None, :], ang, axis=1, mode="wrap")
    sm = lambda Z: convolve1d(convolve1d(Z, ang, axis=1, mode="wrap"), rad, axis=0, mode="nearest")
    Fkd = sm(Fkd) / np.maximum(sm(cnt_kth), 1e-9)
    Dfft = xr.DataArray(Fft, dims=["frequency", "direction"],
                        coords={"frequency": freqs, "direction": dbins})
    Dfkd = xr.DataArray(Fkd, dims=["wavenumber", "direction"],
                        coords={"wavenumber": k_grid, "direction": dbins})
    return dict(Sf_d=Sf_d, klin=klin, Fk_d=Fk_d, dbins=dbins, Dfft=Dfft, Dfkd=Dfkd)


def spreading(da, Xname, lo, hi, blank_lowE=True):
    """Row-normalise a directional spectrum to the directional distribution
    D(X, theta) (int D dtheta = 1). Returns (theta, X[slice], D[slice])."""
    th = da["direction"].values
    dth = np.radians(np.median(np.diff(th)))
    v = da.transpose(Xname, "direction").values.astype(float)
    norm = np.nansum(v, axis=1, keepdims=True) * dth
    D = np.where(norm > 0, v / norm, np.nan)
    if blank_lowE:
        nn = norm[:, 0]
        D[nn < 1e-3 * np.nanmax(nn)] = np.nan
    Xv = da[Xname].values
    s = (Xv >= lo) & (Xv <= hi)
    return th, Xv[s], D[s]

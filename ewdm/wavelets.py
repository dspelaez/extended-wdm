#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

import numpy as np
import xarray as xr

from typing import Union

from .helpers import get_sampling_frequency


# definition of mother wavelets {{{
class Morlet(object):
    """Implements the Morlet wavelet class.

    Note that the input parameters f and f0 are angular frequencies.
    f0 should be more than 0.8 for this function to be correct, its
    default value is f0 = 6.

    Adapted from: @regeirk/pycwt
    """

    def __init__(self, f0=6):
        self._set_f0(f0)
        self.name = 'Morlet'

    def _set_f0(self, f0):
        """
        Sets the Morlet wave number, the degrees of freedom and the
        empirically derived factors for the wavelet bases 
        """
        self.f0 = f0             # Wave number
        self.dofmin = 2          # Minimum degrees of freedom
        if self.f0 == 6:
            self.cdelta = 0.776  # Reconstruction factor
            self.gamma = 2.32    # Decorrelation factor for time averaging
            self.deltaj0 = 0.60  # Factor for scale averaging
        else:
            self.cdelta = -1
            self.gamma = -1
            self.deltaj0 = -1

    @property
    def flambda(self):
        """Fourier wavelength as of Torrence and Compo (1998)."""
        return (4 * np.pi) / (self.f0 + np.sqrt(2 + self.f0 ** 2))

    @property
    def coi(self):
        """e-Folding Time as of Torrence and Compo (1998)."""
        return 1. / np.sqrt(2)

    @property
    def sup(self):
        """Wavelet support defined by the e-Folding time."""
        return 1. / self.coi

    def psi_ft(self, f):
        """Fourier transform of the approximate Morlet wavelet."""
        return (np.pi ** -0.25) * np.exp(-0.5 * (f - self.f0) ** 2)

    def psi(self, t):
        """Morlet wavelet as described in Torrence and Compo (1998)."""
        return (np.pi ** -0.25) * np.exp(1j * self.f0 * t - t ** 2 / 2)

    def cone_of_influence(self, ntime, fs):
        """Determines the cone-of-influence.

        Note that it is returned as a function of time in Fourier periods. Uses
        triangualr Bartlett window with non-zero end-points.
        """
        _coi = (ntime / 2 - np.abs(np.arange(0, ntime) - (ntime - 1) / 2))
        return  self.flambda * self.coi * _coi / fs

# }}}

#  continuous wavelet transform {{{
def cwt(
        data: Union[np.ndarray, xr.DataArray],
        freqs: np.ndarray,
        fs: float = None,
        mother: Morlet = Morlet(6.),
        normalise: bool = False
    ) -> xr.DataArray:
    """This function computes the continuous wavelet transform.

    Args:
        data (p.ndarray, xr.DataArray): Input signal array.
        freqs (np.ndarray): Frequency array.
        fs (float): Sampling frequency. Default is 2 * freqs[-1].
        mother (Morlet): Mother wavelet class. Default is Morlet(6.).
        normalise (bool): Whether to normalise the signal variance.

    Returns:
        xr.DataArray: Continuous wavelet transform of the input signal.
    """

    assert isinstance(data, xr.DataArray) or isinstance(data, np.ndarray)

    ntime = len(data)

    # sampling frequency
    if fs is None:
        fs = get_sampling_frequency(data["time"])

    # wavelet scales using given frequencies.
    scales = 1 / (mother.flambda * freqs)

    # number scales
    nscale = len(scales)

    # normalise with signal variance
    if normalise:
        signal_std = data.std().item()
        data = (data - data.mean()) / data.std()
            
    # fourier frequencies and fourier transform of signal
    fft = np.fft.fft(data)
    omega = 2*np.pi * np.fft.fftfreq(ntime, 1./fs)

    # loop for fill the window and scales of wavelet
    s = scales[:,None]
    w = np.sqrt(s * omega[1] * ntime) * mother.psi_ft(s * omega)

    # convolve window and transformed series
    # fac = np.sqrt(2 / fs / mother.flambda)
    W = np.fft.ifft(fft[None,:] * w, ntime)

    # finally try to compute time array
    try:
        time = data["time"].values
    except:
        time = np.arange(ntime) / fs

    return xr.DataArray(
        data=W, coords=[freqs, time], dims=["frequency", "time"],
        attrs={
            "sampling_rate": fs,
            "normalised": normalise
        }
    )
# }}}

#  cross wavelet transform {{{
def xwt(
        data_x: Union[np.ndarray, list, xr.DataArray],
        data_y: Union[np.ndarray, list, xr.DataArray],
        freqs: np.ndarray,
        fs: float = None,
        mother: Morlet = Morlet(6.)
    ) -> xr.DataArray:
    """
    This function compute the continuous wavelet transform

    Args:
        data_x, data_y (numpy.ndarray, list, xr.DataArray): Input signal arrays.
        freqs (numpy.ndarray): Frequencies
        fs (float): Sampling frequency. Default fs=1
        mother: instance of Mother wavelet class. Default is Morlet(6.).

    Returns:
        Wxy: Cross-wavelet coefficients
    """

    # number of times
    if len(data_x) != len(data_y):
        raise Exception("Time series do not have the same lenght")
    else:
        ntime = len(data_x)

    if fs is None:
        fs = get_sampling_frequency(data_x["time"])

    # compute continuous wavelet transform for each signal
    Wx = cwt(data_x, fs=fs, freqs=freqs, mother=mother)
    Wy = cwt(data_y, fs=fs, freqs=freqs, mother=mother)
    
    # cross wavelet coefficients
    return Wx * Wy.conj()
# }}}

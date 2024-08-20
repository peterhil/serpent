"""DSP module."""

from __future__ import annotations

import numpy as np
from numpy import fft


def fft_spectra(decoded, n=None, *, norm='forward', **kwargs):
	"""FFT spectrum of the data.

	Returns positive frequencies and the spectrum using Numpy FFT module.

	See Numpy FFT docs:
	https://numpy.org/doc/stable/reference/routines.fft.html#module-numpy.fft
	"""
	n = n or len(decoded)

	dft = fft.rfft(decoded, n=n, norm=norm, **kwargs)
	# Could use fftfreq(n, d=1/n) to directly get the freqs
	norm_freqs = fft.rfftfreq(n)

	# Limit to positive frequencies (0 is DC, n / 2 is the Nyquist frequency)
	positive = slice(1, n // 2)
	dft = dft[positive]
	norm_freqs = norm_freqs[positive]

	freqs = norm_freqs * n
	spectra = np.abs(dft)  # TODO Use power spectrum? = abs(dft) ** 2

	return [freqs, spectra]

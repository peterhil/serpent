"""Statistical function."""
from __future__ import annotations

from collections import Counter, OrderedDict

import numpy as np


def autocorrelate(seq, index):
	"""Autocorrelate sequence at index."""
	return np.corrcoef(seq[:-index], seq[index:])[0, 1]


def autocorrelogram(seq, length=256):
	"""Autocorrelate sequence upto length."""
	# TODO Use signal processing correlogram to get the indices directly
	# See: https://stackoverflow.com/a/7981132/470560
	length = min(len(seq) - 1, length)
	return np.array([1] + [
		autocorrelate(seq, index)
		for index in range(1, length)
	])


def ac_peaks(ac, limit=0.05):
	"""Find peaks (=repeat lengths) from autocorrelation sequence."""
	peaks = np.arange(len(ac))[ac > limit]
	values = ac[peaks]

	return OrderedDict(zip(peaks, values))


def count_sorted(items):
	"""Count items and return as sorted array."""
	counts = Counter(items)
	return np.array(sorted(counts.items())).T

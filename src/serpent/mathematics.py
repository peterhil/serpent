"""Mathematical utilities."""

from __future__ import annotations

import numpy as np


def logn(number, base=np.e):
	"""Logarithm of number on some base."""
	return np.log2(number) / np.log2(base)


def magnitude(data, base=64):
	"""Logarithmic magnitude of data in some base."""
	return int(np.ceil(logn(np.max(data), base)))


def normalise(seq):
	"""Normalise data max to be 1."""
	return np.asanyarray(seq) / np.amax(seq)

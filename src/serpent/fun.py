"""Functional utils."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np


def map_array(function, arr, dtype=None):
	"""Naive map on Numpy arrays."""
	# TODO Use np.apply_along_axis or something faster than this
	return np.array(list(map(function, arr)), dtype=dtype)


def str_join(seq: Iterable, joiner='') -> str:
	"""Join a sequence into a string."""
	return joiner.join(seq)

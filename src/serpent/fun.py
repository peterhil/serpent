"""Functional utils."""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import ItemsView, Iterable, Mapping, Sequence

import numpy as np


def is_not_none(x):
	"""Filter helper for when you ONLY want to filter out the None values."""
	return x is not None


def inverse_od(mapping: Mapping) -> OrderedDict:
	"""Swap keys and values on an OrderedDict."""
	return OrderedDict([(v, k) for k, v in mapping.items()])


def map_array(function, arr, dtype=None):
	"""Naive map on Numpy arrays."""
	# TODO Use np.apply_along_axis or something faster than this
	return np.array(list(map(function, arr)), dtype=dtype)


def second(seq: Sequence):
	"""Return second item of a sequence."""
	return seq[1]



def sort_values(items: ItemsView, reverse=False):
	"""Sort mapping items by values."""
	if isinstance(items, Mapping):
		items = items.items()

	return sorted(items, key=second, reverse=reverse)


def str_join(seq: Iterable, joiner='') -> str:
	"""Join a sequence into a string."""
	return joiner.join(seq)

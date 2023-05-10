"""Information theoretic functions."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from serpent.math.basic import logn


def entropy(data: Sequence[int], n=2, base=2):
	"""Shannon entropy of the data.

	Arguments:
	---------
	data: sequence of decoded integers
	n: number of symbols in the data
	base: base of encoding system (default: binary)
	"""
	assert np.min(data) >= 0, 'Expected positive integers'
	assert np.max(data) < n, f'Expected data to fit in {n} symbols'
	propabilities = np.ones(len(data)) / n  # Equal propabilites

	return -np.sum(logn(propabilities, base))

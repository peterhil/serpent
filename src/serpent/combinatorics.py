"""Combinatorics module."""

from __future__ import annotations

import numpy as np


def spread(n: int=16, offset=0):
	"""Spread (or deal) a sequence evenly across N piles and combine them.

	It's exactly like dealing a deck of cards across N players.

	Example:
	-------
	>>> spread(16, offset=1)
	array([15,  0,  4,  8, 12,  1,  5,  9, 13,  2,  6, 10, 14,  3,  7, 11])
	"""
	seq = np.arange(n)
	magn = int(np.ceil(np.log2(np.max(seq))))
	max_n = 2 ** magn

	offsetted = (seq - offset) % max_n
	rounds = np.array(divmod((offsetted * magn), max_n))
	spread = np.apply_along_axis(np.sum, 0, rounds)

	return spread

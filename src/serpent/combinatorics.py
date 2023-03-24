"""Combinatorics module."""

from __future__ import annotations

from collections.abc import Sequence

import more_itertools as mitr
import numpy as np


def spread(seq: Sequence, n: int, *, offset: int=0):
	"""Spread (or deal) a sequence evenly across N piles and combine them.

	It's exactly like dealing a deck of cards across N players.
	The offset parameter is like cutting the deck.

	Example:
	-------
	>>> spread(np.arange(16), n=4, offset=1)
	array([15,  0,  4,  8, 12,  1,  5,  9, 13,  2,  6, 10, 14,  3,  7, 11])

	>>> spread(np.arange(7), n=2, offset=0)
	array([0, 2, 4, 6, 1, 3, 5])
	"""
	hands = mitr.chunked(seq, n)
	shuffled = [list(hand) for hand in mitr.unzip(hands)]
	spread = np.roll(np.concatenate(shuffled), offset)

	return spread

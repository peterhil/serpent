"""Combinatorics module."""

from __future__ import annotations

import itertools as itr
from collections.abc import Sequence

import more_itertools as mit
import numpy as np

from serpent.fun import is_not_none


def spread(seq: Sequence, n: int, offset: int=0):
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
	if len(seq) == 0:
		return seq
	assert 1 <= n <= len(seq), 'Parameter N should be strictly between 1 and len(seq)'

	hands = mit.chunked(seq, n)
	shuffled = [list(hand) for hand in mit.unzip(hands)]
	cut = np.roll(np.concatenate(shuffled), offset)

	return cut


def unspread(seq: Sequence, n: int, offset: int=0):
	"""Unspread a sequence, see spread.

	Example:
	-------
	>>> seq = np.array([15,  0,  4,  8, 12,  1,  5,  9, 13,  2,  6, 10, 14,  3,  7, 11])
	>>> unspread(seq, n=4, offset=1)
	array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15])

	>>> unspread(spread(np.arange(7), n=2, offset=2), n=2, offset=2)
	array([0, 1, 2, 3, 4, 5, 6])

	"""
	if len(seq) == 0:
		return seq
	assert 1 <= n <= len(seq), 'Parameter N should be strictly between 1 and len(seq)'

	uncut = np.roll(seq, -offset)
	# FIXME: Oh, why is Python so hostile to functional programming?
	# Fix this on a better day... and make it clean. (just compare to spread FFS!)
	unshuffled = [*itr.zip_longest(*mit.chunked(uncut, int(np.ceil(len(seq) / n))))]
	piles = [list(filter(is_not_none, pile)) for pile in unshuffled]
	deck = np.concatenate(piles)

	return deck

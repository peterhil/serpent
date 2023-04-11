"""Statistical function."""
from __future__ import annotations

from collections import Counter, OrderedDict, defaultdict
from collections.abc import Callable, Iterable, Sequence

import numpy as np

from serpent.mathematics import logn


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


def symbol_run_lengths(data: Iterable[str]):
	"""Statistics of repeat lengths for each symbol of a sequence."""
	last = OrderedDict()
	stat = defaultdict(list)
	max_length = 0

	for index, sym in enumerate(data):
		previous = last.get(sym, 0)
		last[sym] = index
		run = index - previous
		max_length = max(run, max_length)
		stat[sym].append(run)

	return stat, max_length


def quasar(data: Iterable[str], fn: Callable):
	"""Map a callback function through symbol run lengths."""
	stat, max_length = symbol_run_lengths(data)
	result = OrderedDict(sorted([
		(k, fn(np.array(s), max_length))
		for k, s in stat.items()
		if len(s)
	]))

	return result

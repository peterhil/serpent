"""Statistical function."""
from __future__ import annotations

from collections import Counter, OrderedDict, defaultdict
from collections.abc import Callable, Iterable, Sequence

import numpy as np

from serpent.mathematics import logn
from serpent.padding import pad_end


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


def gc_content(counts: Counter) -> float:
	"""Get GC content from nucleotide counts."""
	return counts['G'] + counts['C']


def pulse_repetition_intervals(data: Iterable[str]):
	"""Pulse repetition intervals for each symbol of data.

	See: https://en.wikipedia.org/wiki/Pulse_width
	"""
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


def map_pulses(data: Iterable[str], fn: Callable, key=None):
	"""Map a callback function through symbol run lengths."""
	stat, max_length = pulse_repetition_intervals(data)
	result = OrderedDict(sorted([
		(k, fn(np.array(s), max_length))
		for k, s in stat.items()
		if len(s)
	], key=key))

	return result


def quasar_pulses(data, cumulative=False, key=None):
	[stat, scale] = pulse_repetition_intervals(data)
	height = np.amax(np.array([len(s) for s in stat.values()]))

	def to_pulses(s):
		rhythm = np.cumsum(s) if cumulative else s
		return pad_end(rhythm, 0, n=height)

	if key is not None:
		pulses = OrderedDict([
			(k, to_pulses(np.array(stat.get(k, [-1]))))
			for k in key
		])

		# Check that key includes all symbols in the data
		if set(stat.keys()) > set(key):
			extra = set(stat.keys()) - set(key)
			err_msg = f'Insufficient key, needs to include: {extra}'
			raise ValueError(err_msg)
	else:
		pulses = OrderedDict(sorted([
			(k, to_pulses(np.array(s)))
			for k, s in stat.items()
			if len(s)
		]))

	return pulses, height, scale

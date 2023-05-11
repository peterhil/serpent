"""Information theoretic functions."""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Sequence

import numpy as np

from serpent.fun import str_join
from serpent.math.basic import logn, percent


def ensure_counter(data: Counter | Sequence) -> Counter:
	return data if isinstance(data, Counter) else Counter(data)


def abs_rate(data: Counter | Sequence) -> int:
	"""Absolute rate: unique symbols in the data."""
	stats = ensure_counter(data)

	return len(stats)


def entropy(data: Counter | Sequence, base: float=2) -> float:
	"""Shannon entropy of the data.

	Arguments:
	---------
	data: Counter or Sequence
	base: unit of information (common choices: 2, e, 10...)
	"""
	stats = ensure_counter(data)
	counts = np.array([v for v in stats.values() if v != 0])

	total = np.sum(counts)
	probabilities = counts / total
	entropy = np.sum(-probabilities * logn(probabilities, base))

	return entropy


def abs_redundancy(data: Counter | Sequence, base: float=2) -> float:
	"""Absolute redundancy."""
	counts = ensure_counter(data)
	rate = abs_rate(counts)

	return rate - entropy(counts, base)


def rel_redundancy(data: Counter | Sequence, base: float=2) -> float:
	"""Relative redundancy."""
	counts = ensure_counter(data)
	rate = abs_rate(counts)

	return abs_redundancy(counts, base) / rate


def efficiency(data: Counter | Sequence, base: float=2) -> float:
	"""Efficiency."""
	counts = ensure_counter(data)
	rate = abs_rate(counts)
	eff2 = entropy(counts, base) / rate
	efficiency = 1.0 - rel_redundancy(counts, base)

	assert efficiency == eff2, 'Definitions not equal'
	return efficiency


def max_compr_ratio(data: Counter | Sequence, base: float=2) -> float:
	"""Maximum compression ratio."""
	counts = ensure_counter(data)
	rate = abs_rate(counts)

	return rate / entropy(counts, base)


def statistics_header():
	return str_join([
		'Rate',
		'Entropy',
		'Perplex',
		'Inform.',
		'Length',
		'Abs.rdn',
		'Rel.rdn',
		'Eff.',
		'Compmax',
	], '\t')


def statistics(data: Counter | Sequence, base: float=2) -> float:
	"""Information statistics."""
	c = ensure_counter(data)

	rate = abs_rate(c)
	entr = entropy(c, base)
	total = np.sum([*c.values()])
	info = entr * total

	abs_red = rate - entr
	rel_red = abs_red / rate
	eff = 1.0 - rel_red
	max_compr = rate / entr
	perplex = base ** entr

	# TODO Split into format_info_stats
	return str_join([
		f'{rate}',
		f'{entr :.2f}',
		f'{perplex :.2f}',
		f'{math.ceil(info)}',
		f'{total}',
		f'{abs_red :.2f}',
		f'{percent(rel_red, 2) :.2f}%',
		f'{percent(eff, 2) :.2f}%',
		f'{max_compr :.2f}',
	], '\t')

"""Information theoretic functions."""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence

import numpy as np

from serpent.math.basic import logn


def abs_rate(data: Counter | Sequence) -> int:
	"""Absolute rate: unique symbols in the data."""
	stats = data if isinstance(data, Counter) else Counter(data)

	return len(stats)


def entropy(data: Counter | Sequence, base: float=2) -> float:
	"""Shannon entropy of the data.

	Arguments:
	---------
	data: Counter or Sequence
	base: unit of information (common choices: 2, e, 10...)
	"""
	stats = data if isinstance(data, Counter) else Counter(data)
	counts = np.array([v for v in stats.values() if v != 0])
	total = np.sum(counts)
	probabilities = counts / total

	entropy = np.sum(-probabilities * logn(probabilities, base))

	return entropy

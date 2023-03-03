"""Mathematical utilities."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Union

import numpy as np
from numpy import dtype, ndarray

LogBase = Union[int, float]
Numeric = Union[int, float, complex]
NpNumeric = ndarray[Numeric, dtype[Any]]
NumericSeq = Union[Sequence[Numeric], NpNumeric]


def logn(number: Numeric, base: LogBase=np.e) -> float:
	"""Logarithm of number on some base."""
	return np.log2(number) / np.log2(base)


def magnitude(data: NumericSeq, base: LogBase=64) -> int:
	"""Logarithmic magnitude of data in some base."""
	return int(np.ceil(logn(np.max(data), base)))


def normalise(seq: NumericSeq) -> NumericSeq:
	"""Normalise data max to be 1."""
	return np.asanyarray(seq) / np.amax(seq)

"""Mathematical utilities."""

from __future__ import annotations

from typing import TypeVar, Union

import numpy as np
from numpy import dtype, ndarray
from numpy.typing import ArrayLike, NBitBase

T = TypeVar('T', bound=NBitBase)


LogBase = Union[int, float]
NumericSeq = ArrayLike


def logn(number: NumericSeq, base: LogBase=np.e) -> float:
	"""Logarithm of number on some base."""
	return float(np.log2(number) / np.log2(base))


def magnitude(data: NumericSeq, base: LogBase=64) -> int:
	"""Logarithmic magnitude of data in some base."""
	return int(np.ceil(logn(np.max(data), base)))


def normalise(seq: np.floating[T]) -> ndarray[T, dtype[np.floating[T]]]:
	"""Normalise data max to be 1."""
	return np.asanyarray(seq) / np.amax(seq)

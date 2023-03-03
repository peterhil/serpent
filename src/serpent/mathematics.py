"""Mathematical utilities."""

from __future__ import annotations

from typing import TypeVar, Union

import numpy as np
import numpy.typing as npt
from numpy import dtype, ndarray

T = TypeVar('T', bound=npt.NBitBase)


LogBase = Union[int, float]
Numeric = Union[int, float, complex]
NumericSeq = npt.ArrayLike


def logn(number: Numeric, base: LogBase=np.e) -> float:
	"""Logarithm of number on some base."""
	return float(np.log2(number) / np.log2(base))


def magnitude(data: NumericSeq, base: LogBase=64) -> int:
	"""Logarithmic magnitude of data in some base."""
	return int(np.ceil(logn(np.max(data), base)))


def normalise(seq: np.floating[T]) -> ndarray[T, dtype[np.floating[T]]]:
	"""Normalise data max to be 1."""
	return np.asanyarray(seq) / np.amax(seq)

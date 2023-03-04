"""Mathematical utilities."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Union

import numpy as np
from numpy.typing import ArrayLike, NDArray

LogBase = Union[int, float]
Numeric = Union[int, float, complex]
NumericSeq = Union[Sequence[Numeric], ArrayLike]


def logn(number: NumericSeq, base: LogBase=np.e) -> float:
	"""Logarithm of number on some base."""
	return float(np.log2(number) / np.log2(base))


def magnitude(data: NumericSeq, base: LogBase=64) -> int:
	"""Logarithmic magnitude of data in some base."""
	return int(np.ceil(logn(np.max(data), base)))


def normalise(seq: NumericSeq) -> NDArray[np.float64]:
	"""Normalise data max to be 1."""
	return np.asanyarray(seq, dtype=np.float64) / np.amax(seq)

"""Mathematical utilities."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Union

import numpy as np
from numpy.typing import ArrayLike, NDArray

LogBase = Union[int, float]
Numeric = Union[int, float, complex]
NumericSeq = Union[Sequence[Numeric], ArrayLike]


phi = (1 + np.sqrt(5)) / 2
phi_large = 1 / phi
phi_small = 1 - phi_large


def autowidth(n, base=64, aspect=phi):
	"""Automatic image width rounded to a multiple of some base width.

	Aspect ratio should be a positive float greater than 0, and
	can be used to give landscape (>1) or portrait (<1) images.

	Width will be at least once the base.
	"""
	unrounded_width = np.sqrt(n * aspect)
	multiple = max(1, int(np.round(unrounded_width / base)))

	return int(base * multiple)


def logn(number: NumericSeq, base: LogBase=np.e) -> float:
	"""Logarithm of number on some base."""
	return float(np.log2(number) / np.log2(base))


def magnitude(data: NumericSeq, base: LogBase=64) -> int:
	"""Logarithmic magnitude of data in some base."""
	return int(np.ceil(logn(np.max(data), base)))


def normalise(seq: NumericSeq) -> NDArray[np.float64]:
	"""Normalise data max to be 1."""
	return np.asanyarray(seq, dtype=np.float64) / np.amax(seq)

"""Padding utilities."""
from __future__ import annotations

from typing import TypeVar

import numpy as np
from numpy.typing import NDArray

T = TypeVar('T', bound=np.generic)


def get_padding(data: NDArray[T], fill: T, n: int=3) -> NDArray[T]:
	"""Get padding for the divisor.

	The length of data and padding will be evenly divisible by the divisor.
	"""
	padding = []
	rem = len(data) % n
	if rem != 0:
		pad_length = n - rem
		padding = pad_length * [fill]

	return np.array(padding)


def pad_to_left(data: NDArray[T], fill: T, *, n: int=3) -> NDArray[T]:
	"""Pad data with the fill characters on the end."""
	padding = get_padding(data, fill, n)

	return np.concatenate([data, padding])


def pad_to_right(data: NDArray[T], fill: T, *, n: int=3) -> NDArray[T]:
	"""Pad data with the fill character on the beginning."""
	padding = get_padding(data, fill, n)

	return np.concatenate([padding, data])

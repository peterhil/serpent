"""Padding utilities."""
from __future__ import annotations

from typing import TypeVar

T = TypeVar('T', str, int)


def get_padding(data: list[T], fill: T, n: int=3) -> list[T]:
	"""Get padding for the divisor.

	The length of data and padding will be evenly divisible by the divisor.
	"""
	padding = []
	rem = len(data) % n
	if rem != 0:
		pad_length = n - rem
		padding = pad_length * [fill]

	return padding


def pad_to_left(data: list[T], fill: T, *, n: int=3) -> list[T]:
	"""Pad data with the fill characters on the end."""
	padding: list[T] = get_padding(data, fill, n)

	return list(data) + padding


def pad_to_right(data: list[T], fill: T, *, n: int=3) -> list[T]:
	"""Pad data with the fill character on the beginning."""
	padding: list[T] = get_padding(data, fill, n)

	return padding + list(data)

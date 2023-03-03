"""Padding utilities."""
from __future__ import annotations


def get_padding(data, divisor=3, fill="A"):
	"""Get padding for the divisor.

	The length of data and padding wil be evenly divisible by the divisor.
	"""
	padding = []
	rem = len(data) % divisor
	if rem != 0:
		pad_length = divisor - rem
		padding = pad_length * [fill]

	return padding


def pad_to_left(data, divisor=3, fill="A"):
	"""Pad data with the fill characters on the end."""
	padding = get_padding(data, divisor, fill)

	return list(data) + padding


def pad_to_right(data, divisor=3, fill="A"):
	"""Pad data with the fill character on the beginning."""
	padding = get_padding(data, divisor, fill)

	return padding + list(data)

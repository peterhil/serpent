#!/usr/bin/env python


def get_padding(data, divisor=3, fill="A"):
	"""Return suitable padding in order to pad the data into a length
	evenly divisible by the divisor."""
	padding = []
	rem = len(data) % divisor
	if rem != 0:
		pad_length = divisor - rem
		padding = pad_length * [fill]

	return padding


def pad_to_left(data, divisor=3, fill="A"):
	"""Pad data to a length divisible by the divisor with the fill
	characters on the end"""
	padding = get_padding(data, divisor, fill)

	return list(data) + padding


def pad_to_right(data, divisor=3, fill="A"):
	"""Pad data to a length divisible by the divisor with the fill
	character on the beginning"""
	padding = get_padding(data, divisor, fill)

	return padding + list(data)

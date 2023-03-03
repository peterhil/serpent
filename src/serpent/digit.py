"""Digit conversions into various numerical bases."""
from __future__ import annotations


def digits_to_number(seq, base=64):
	"""Convert digits in sequence into a number in given base."""
	number = 0
	for i, digit in enumerate(reversed(seq)):
		number = number + base ** i * digit

	return number


def number_to_digits(number, base=64):
	"""Convert a number in given base into digits in a sequence."""
	result = []
	remainder = number

	while True:
		[multiplier, remainder] = divmod(remainder, base)
		result.append(remainder)
		remainder = multiplier
		if multiplier == 0:
			break

	return list(reversed(result))

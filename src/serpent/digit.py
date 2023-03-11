"""Digit conversions into various numerical bases."""
from __future__ import annotations

__all__ = [
	'digits_to_num',
	'num_to_digits',
]

def digits_to_num(seq, base=64):
	"""Convert digits sequence into a number in given base."""
	number = 0
	for i, digit in enumerate(reversed(seq)):
		number = number + base ** i * digit

	return number


def num_to_digits(number, base=64):
	"""Convert a number in given base into digits sequence."""
	result = []
	remainder = number

	while True:
		[multiplier, remainder] = divmod(remainder, base)
		result.append(remainder)
		remainder = multiplier
		if multiplier == 0:
			break

	return list(reversed(result))

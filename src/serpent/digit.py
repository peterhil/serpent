from __future__ import annotations


def digits_to_number(seq, base=64):
	value = 0
	for i, n in enumerate(reversed(seq)):
		value = value + base**i * n

	return value


def number_to_digits(number, base=64):
	result = []
	rem = number

	while True:
		[mul, rem] = divmod(rem, base)
		result.append(rem)
		rem = mul
		if mul == 0:
			break

	return list(reversed(result))

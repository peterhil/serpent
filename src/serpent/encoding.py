from __future__ import annotations

from itertools import combinations

from serpent.fun import map_array


def char_range(char1, char2):
	"""Generate range of characters from `char1` to `char2`, inclusive."""
	for ch in range(ord(char1), ord(char2) + 1):
		yield chr(ch)


# Base 64 alphabet variant.
# This is the most readable of the different options,
# so do not change easily.
base64 = (
	list(char_range("A", "Z")) +
	list(char_range("a", "z")) +
	list(char_range("0", "9")) +
	["+", "."]
)

alphabet64 = dict({(i, char) for (i, char) in enumerate(base64)})
combos = map_array(lambda cm: "".join(cm), combinations(base64, 2))

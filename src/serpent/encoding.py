"""Base64 variant and other encodings of DNA data."""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Iterable

from serpent.fun import inverse_od, str_join


def char_range(start: str, end: str) -> Iterable[str]:
	"""Generate range of characters from `start` to `end`, inclusive."""
	for char in range(ord(start), ord(end) + 1):
		yield chr(char)


# Base 64 alphabet variant.
# This is the most readable of the different options,
# so do not change easily.
BASE64 = str_join(
	list(char_range("A", "Z")) +
	list(char_range("a", "z")) +
	list(char_range("0", "9")) +
	["+", "."]
)

num_to_base64 = OrderedDict(enumerate(BASE64))
base64_to_num = inverse_od(num_to_base64)

"""Convert data between numbers and base64 encoding."""

from __future__ import annotations

from collections import OrderedDict

from serpent.encoding import BASE64
from serpent.fun import inverse_od

__all__ = [
	'base64_to_num',
	'num_to_base64',
]

num_to_base64 = OrderedDict(enumerate(BASE64))
base64_to_num = inverse_od(num_to_base64)

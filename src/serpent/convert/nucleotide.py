"""Convert data between nucleotide bases and numbers."""

from __future__ import annotations

from collections import OrderedDict

from serpent.config import BASE_ORDER
from serpent.fun import inverse_od

__all__ = [
	'nt_to_num',
	'num_to_nt',
]

# TODO: Mapping order and numbering could be distinct if using a gray code like
# scheme. For example: With GACT/GACU order the numbering could be 2013, which
# would be the same as ACGT in linear ordering.
#
# ACGT bases:
# A 00 0
# C 01 1
# G 10 2
# T 11 3
#
# GACT bases:
# G 10 2
# A 00 0
# C 01 1
# T 11 3

num_to_nt = OrderedDict(enumerate(BASE_ORDER))
nt_to_num = inverse_od(num_to_nt)

"""Convert data between codons and numbers."""

from __future__ import annotations

import itertools as itr
from collections import OrderedDict

import numpy as np

from serpent.config import BASE_ORDER
from serpent.convert.nucleotide import nt_to_num
from serpent.fun import inverse_od, str_join


def codons_array(bases=BASE_ORDER):
	"""Sequence of the 64 codons using the given order."""
	return np.fromiter(map(str_join, itr.product(bases, repeat=3)), dtype='<U3')


codons_inverse = OrderedDict([*enumerate(codons_array(BASE_ORDER))])
codons = inverse_od(codons_inverse)


def decode_codon(codon: str) -> int:
	"""Decode a codon string into a a number between 0 and 63."""
	# TODO Remove if and when codon_to_num works
	result = 0
	for num, char in enumerate(reversed(codon)):
		result += nt_to_num[char] << num * 2

	return result


def codon_to_num(codon: str) -> int:
	"""Decode a single codon (three letter string) into a number between 0 and 63."""
	return codons[codon]


def num_to_codon(code: int) -> str:
	"""Encode a number between 0 and 63 into a codon."""
	return codons_inverse[code]

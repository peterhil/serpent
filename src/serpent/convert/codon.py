"""Convert data between codons and numbers."""

from __future__ import annotations

import itertools as itr
from collections import OrderedDict

import numpy as np

from serpent.fun import inverse_od, str_join
from serpent.settings import BASE_ORDER

CODONS_LEN = 64


def base_orders():
	"""All permutations of the four nucleotide bases (ACGT)."""
	return [*map(str_join, itr.permutations('ACGT', 4))]


def valid_base_order(bases):
	"""Check that given base order is valid."""
	return ''.join(sorted(bases)) == 'ACGT'


def codons_array(bases=BASE_ORDER):
	"""Sequence of the 64 codons using the given order."""
	assert valid_base_order(bases), 'Invalid nucleotide base order.'

	return np.fromiter(map(str_join, itr.product(bases, repeat=3)), dtype='<U3')


# TODO Maybe do all the base orders like with genetic code?
codons_inverse = OrderedDict([*enumerate(codons_array(BASE_ORDER))])
codons = inverse_od(codons_inverse)


def codon_to_num(codon: str) -> int:
	"""Decode a single codon (three letter string) into a number between 0 and 63."""
	return codons[codon]


def num_to_codon(code: int) -> str:
	"""Encode a number between 0 and 63 into a codon."""
	return codons_inverse[code]

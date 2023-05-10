"""Degenerate nucleotide conversions."""

from __future__ import annotations

from collections import OrderedDict
from functools import reduce
from operator import or_
from typing import Union

from serpent.fun import inverse_od, second, str_join
from serpent.math.bits import popcount
from serpent.math.combinatorics import unspread
from serpent.settings import BASE_ORDER

Dnt = str
DntBits = int
DntSet = Union[set, str]

def is_degenerate(nts: str, bases: str=BASE_ORDER) -> bool:
	"""Return true if any nucleotide is degenerate in the given string.

	Accepts a string of nucleotides or single codon.
	"""
	return not (set(nts) <= set(bases))


def binomial_dnt(bases=BASE_ORDER):
	# TODO Ignores base order, change non-generate bases to use powers of two
	# (1, 2, 4, 8) and use them here
	binomial_order = [
		# 0 bases (1)
		'Z', #: 0,  # Zero
		# 1 base (4)
		'G', #: G,
		'A', #: A,
		'C', #: C,
		'T', #: T,
		# 2 bases (6)
		'R', #: G | A,  # puRine
		'S', #: G | C,  # Strong
		'K', #: G | T,  # Keto
		'M', #: A | C,  # aMino
		'W', #: A | T,  # Weak
		'Y', #: C | T,  # pYrimidine
		# 3 bases (4)
		'V', #: G | A | C | 0,  # Not T
		'D', #: G | A | 0 | T,  # Not C
		'B', #: G | 0 | C | T,  # Not A
		'H', #: 0 | A | C | T,  # Not G
		# 4 bases (1)
		'N', #: G | A | C | T,  # Any base
	]
	# Spread bases into positions 0, 4, 8, 12 (multiples of base length)
	evenly_spread_bases = unspread(binomial_order, len(bases), 1)
	mapping = OrderedDict([*enumerate(evenly_spread_bases)])

	return mapping


def inverse_exp_dnt(bases=BASE_ORDER):
	# TODO Does not match binomial_dnt
	idnt = OrderedDict([(base, 2 ** n) for n, base in enumerate(bases)])
	[G, A, C, T] = [idnt['G'], idnt['A'], idnt['C'], idnt['T']]

	idnt.update({
		'Z': 0,  # Zero
		# 2 bases
		'R': G | A,  # puRine
		'Y': C | T,  # pYrimidine
		'S': G | C,  # Strong
		'W': A | T,  # Weak
		'M': A | C,  # aMino
		'K': G | T,  # Keto
		# 3 bases
		'H': 0 | A | C | T,  # Not G
		'B': G | 0 | C | T,  # Not A
		'D': G | A | 0 | T,  # Not C
		'V': G | A | C | 0,  # Not T
		# 4 bases
		'N': G | A | C | T,  # Any base
	})

	idnt = OrderedDict([*sorted(idnt.items(), key=second)])
	return idnt


dnt_binomial = binomial_dnt(BASE_ORDER)
inv_dnt_binomial = inverse_od(dnt_binomial)

dnt_to_bits = inverse_exp_dnt(BASE_ORDER)
bits_to_dnt = inverse_od(dnt_to_bits)


def num_to_dnt(num: int) -> str:
	"""Convert number (<16) to degenerate nucleotide base."""
	return dnt_binomial[num]


def dnt_to_num(dnt: str) -> int:
	"""Convert degenerate nucleotide base to number (<16)."""
	return inv_dnt_binomial[dnt]


def dntset_to_bits(dntset: DntSet) -> DntBits:
	"""Convert a set of nucleotides to an integer representing four bits.

	Returns the highest DNT bits representing all the DNTs in given DNT set.

	The bits represent the bases in BASE_ORDER, for example:
	{G: 1, A: 2, C: 4, T: 8}

	# Few examples
	>>> dntset_to_bits('GACT')
	15
	>>> dntset_to_bits('GAC')
	7
	>>> dntset_to_bits('R')
	3

	# Test powerset of bases
	>>> from serpent.settings import BASE_ORDER
	>>> import more_itertools as mit
	>>>
	>>> sets = [set(nts) for nts in mit.powerset(BASE_ORDER)]
	>>> [dntset_to_bits(s) for s in sets]
	[0, 1, 2, 4, 8, 3, 5, 9, 6, 10, 12, 7, 11, 13, 14, 15]
	"""
	return reduce(or_, (dnt_to_bits[dnt] for dnt in set(dntset)), 0)


def compress_dntset(dntset: DntSet) -> Dnt:
	"""Contract a set of DNTs into more generic DNT that covers the whole set.

	Complementary to decompress_dnt.
	"""
	return bits_to_dnt.get(dntset_to_bits(dntset))


def decompress_dnt(dnt: Dnt, bases=BASE_ORDER) -> set:
	"""Expand DNT to nucleotide bases (G, A, C, T).

	Complementary to compress_dntset.
	"""
	assert dnt in dnt_to_bits, f'DNT not found: {dnt}'
	return {str_join(nt) for nt in bases if dnt_include(nt, dnt or 'Z')}


def bits_include(bits: DntBits, mask: DntBits) -> bool:
	return ~mask & bits == 0


def dnt_include(dnt: Dnt, mask: Dnt) -> bool:
	bits = dnt_to_bits[dnt]
	mask = dnt_to_bits[mask]
	return bits_include(bits, mask)


def dnt_degree(dnt: str) -> int:
	"""Degree of degeneracy for a degenerate nucleotide symbol."""
	bits = dnt_to_bits[dnt]
	return popcount(bits)

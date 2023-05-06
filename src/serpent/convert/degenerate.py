"""Convert data between degenerate nucleotide bases, codons, amino acids and numbers."""

from __future__ import annotations

import itertools as itr
from collections import OrderedDict

import numpy as np

from serpent.fun import inverse_od, second, str_join
from serpent.math.combinatorics import unspread
from serpent.settings import BASE_ORDER


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
	# TODO Does not match create_dnt currently, but kept as a reference
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


def inverse_quaternary_dnt(bases=BASE_ORDER):
	"""Degenerate nucleotide mapping in arbitrary quaternary base.

	For GACT base order the mapping is simply:
	00-03: 'GACT'
	04-07: 'RWSY'
	08-11: 'HBDV'
	12-15: 'ZMKN'
	"""
	# TODO Enable larger bases?
	idnt = inverse_od(OrderedDict(enumerate(bases[:4])))

	for high, group in zip([4, 8, 12], ['RSWY', 'HBDV', 'ZMKN']):
		idnt.update([
			(sym, high + idnt.get(nt))
			for sym, nt in zip(list(group), list('GACT'))
		])
	idnt = OrderedDict([*sorted(idnt.items(), key=second)])

	return idnt


# Bases

degenerate = binomial_dnt(BASE_ORDER)
inv_degenerate = inverse_od(degenerate)
inv_degenerate_exp = inverse_exp_dnt(BASE_ORDER)
degenerate_exp = inverse_od(inv_degenerate_exp)
# inv_degenerate = inverse_quaternary_dnt(BASE_ORDER)
# degenerate = inverse_od(inv_degenerate)


def num_to_dnt(num: int) -> str:
	"""Convert number (<16) to degenerate nucleotide base."""
	return degenerate[num]


def dnt_to_num(dnt: str) -> int:
	"""Convert degenerate nucleotide base to number (<16)."""
	return inv_degenerate[dnt]


def degenerate_codons():
	symbols = str_join(degenerate.values())
	degenerates = map(str_join, itr.product(symbols, repeat=3))

	return np.array([*degenerates], dtype='U3')


def is_degenerate(codon: str, bases: str=BASE_ORDER) -> bool:
	return not (set(codon) < set(bases))


# Codons

inv_degen_codons = OrderedDict([*enumerate(degenerate_codons())])
degen_codons = inverse_od(inv_degen_codons)

CODON_LEN = 3
DEGEN_MAX = 4096

def degen_to_num(codon: str) -> int:
	"""Decode a single degenerate codon (three letter string) into a number (<4096)."""
	assert len(codon) == CODON_LEN, 'Invalid codon length.'
	return degen_codons[codon]


def num_to_degen(code: int) -> str:
	"""Encode a number between (<4096) into a degenerate codon."""
	assert code < DEGEN_MAX, 'Invalid degenerate codon number.'
	return inv_degen_codons[code]

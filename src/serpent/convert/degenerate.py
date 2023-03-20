"""Convert data between degenerate nucleotide bases, codons, amino acids and numbers."""

from __future__ import annotations

import itertools as itr
from collections import OrderedDict

import numpy as np

from serpent.fun import inverse_od, second, str_join
from serpent.settings import BASE_ORDER


def create_inverse_dnt(bases=BASE_ORDER):
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

# Bases

inv_degenerate = create_inverse_dnt(BASE_ORDER)
degenerate = inverse_od(inv_degenerate)


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

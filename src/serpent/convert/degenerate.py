"""Degenerate codon conversions."""

from __future__ import annotations

import itertools as itr
from collections import OrderedDict

import numpy as np

from serpent.convert.dnt import degenerate
from serpent.fun import inverse_od, str_join

CODON_LEN = 3
DEGEN_MAX = 4096


def degenerate_codons():
	symbols = str_join(degenerate.values())
	degenerates = map(str_join, itr.product(symbols, repeat=3))

	return np.array([*degenerates], dtype='U3')


inv_degen_codons = OrderedDict([*enumerate(degenerate_codons())])
degen_codons = inverse_od(inv_degen_codons)


def degen_to_num(codon: str) -> int:
	"""Decode a single degenerate codon (three letter string) into a number (<4096)."""
	assert len(codon) == CODON_LEN, 'Invalid codon length.'
	return degen_codons[codon]


def num_to_degen(code: int) -> str:
	"""Encode a number between (<4096) into a degenerate codon."""
	assert code < DEGEN_MAX, 'Invalid degenerate codon number.'
	return inv_degen_codons[code]

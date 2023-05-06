"""Degenerate codon conversions."""

from __future__ import annotations

import itertools as itr
from collections import OrderedDict
from operator import itemgetter

import more_itertools as mit
import numpy as np

from serpent.convert.dnt import compress_dntset, dnt_binomial
from serpent.convert.genetic_code import STANDARD_CODONS, STANDARD_TABLES
from serpent.fun import inverse_od, str_join

CODON_LEN = 3
DEGEN_MAX = 4096


def degenerate_codons():
	symbols = str_join(dnt_binomial.values())
	degenerates = map(str_join, itr.product(symbols, repeat=3))

	return np.array([*degenerates], dtype='U3')


inv_degen_codons = OrderedDict([*enumerate(degenerate_codons())])
degen_codons = inverse_od(inv_degen_codons)


def snp_to_degen(codons: list[str]) -> str:
	"""Convert a list of single nucleotide molymorphic (SNP) codons to a degenerate codon.

	>>> snp_to_degen(['ACG', 'ACA', 'ACC', 'ACT'])
	'ACN'
	>>> snp_to_degen(['ATA', 'ATC', 'ATT'])
	'ATH'
	"""
	# TODO Assert SNP property of codons in the list
	transposed = [list(g) for g in mit.unzip(codons)]
	return str_join([compress_dntset(set(bases)) for bases in transposed])


def create_degen_to_amino_map(table: int=1) -> str:
	"""Get a mapping for degenerate codon to an amino acid conversion."""
	# TODO Memoize or precalculcate the mappings for each table
	aminos = STANDARD_TABLES[table].raw
	codons = STANDARD_CODONS
	iterable = zip(aminos, codons)
	grouper = mit.groupby_transform(iterable, itemgetter(0), itemgetter(1), snp_to_degen)

	# TODO Handle degenerate B, J, and Z amino acid IUPAC codes:
	# Update degenerate_amino on convert.amino to handle all tables
	mapping = OrderedDict([(degen, amino) for amino, degen in grouper])

	return mapping


def degen_to_amino(degen: str, table: int=1) -> str:
	"""Convert a degenerate codon to an amino acid."""
	mapping = create_degen_to_amino_map(table)
	# FIXME Use trie with a custom key, as now this only finds the exact degenerate codons
	return mapping.get(degen, 'X')


def degen_to_num(codon: str) -> int:
	"""Decode a single degenerate codon (three letter string) into a number (<4096)."""
	assert len(codon) == CODON_LEN, 'Invalid codon length.'
	return degen_codons[codon]


def num_to_degen(code: int) -> str:
	"""Encode a number between (<4096) into a degenerate codon."""
	assert code < DEGEN_MAX, 'Invalid degenerate codon number.'
	return inv_degen_codons[code]

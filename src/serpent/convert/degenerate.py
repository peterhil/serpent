"""Degenerate codon conversions."""

from __future__ import annotations

import itertools as itr
from collections import OrderedDict
from functools import cache, lru_cache
from operator import itemgetter

import more_itertools as mit
import numpy as np
from pytrie import SortedStringTrie as Trie

from serpent.convert.dnt import (
	compress_dntset,
	decompress_dnt,
	dnt_binomial,
	dnt_degree,
	dnt_include,
)
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
	"""Convert a list of SNP* codons to a degenerate codon.

	* = single nucleotide polymorphic

	>>> snp_to_degen(['ACG', 'ACA', 'ACC', 'ACT'])
	'ACN'
	>>> snp_to_degen(['ATA', 'ATC', 'ATT'])
	'ATH'

	"""
	transposed = [str_join(list(g)) for g in mit.unzip(codons)]

	# Assert SNP property of the codons
	length_check = set(sorted([len(set(g)) for g in transposed])[:-1])
	err_msg = f'Codons should differ in only one position, got: {codons}'
	assert length_check == {1}, err_msg

	return str_join([compress_dntset(set(bases)) for bases in transposed])


@cache
def create_degen_to_amino_map(table: int=1) -> str:
	"""Get a mapping for degenerate codon to an amino acid conversion."""
	aminos = STANDARD_TABLES[table].raw
	codons = STANDARD_CODONS
	iterable = zip(aminos, codons, strict=False)
	grouper = mit.groupby_transform(
		iterable,
		itemgetter(0),
		itemgetter(1),
		snp_to_degen
	)

	# TODO Handle degenerate B, J, and Z amino acid IUPAC codes:
	# Update degenerate_amino on convert.amino to handle all tables
	mapping = OrderedDict([(degen, amino) for amino, degen in grouper])

	return mapping


def degen_to_amino(degen: str, table: int=1) -> str:
	"""Convert a degenerate codon to an amino acid."""
	mapping = create_degen_to_amino_map(table)

	t = Trie(mapping)
	items = t.items(degen[:2])

	matching = [amino for dg, amino in items if dnt_include(degen[2], dg[2])]
	amino = matching[0] if matching else 'X'

	return amino


@lru_cache(maxsize=3375)
def degen_to_aminoset(degen: str, table: int=1) -> set:
	"""Return set of amino acids possibly represented by a degenerate codon."""
	nt_sets = [decompress_dnt(dnt) for dnt in degen]
	codons = [str_join(nts) for nts in itr.product(*nt_sets)]

	return {degen_to_amino(dgn, table) for dgn in codons}


def degen_to_num(codon: str) -> int:
	"""Decode a single degenerate codon (three letter string) into a number (<4096)."""
	assert len(codon) == CODON_LEN, 'Invalid codon length.'
	return degen_codons[codon]


def num_to_degen(code: int) -> str:
	"""Encode a number between (<4096) into a degenerate codon."""
	assert code < DEGEN_MAX, 'Invalid degenerate codon number.'
	return inv_degen_codons[code]


def degeneracy(degen: str):
	"""Degree of degeneracy for a degenerate codon."""
	degree = np.prod([dnt_degree(dnt) for dnt in degen])
	return degree

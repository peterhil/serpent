from __future__ import annotations

from collections import OrderedDict

from serpent.convert.codon import codon_to_num, num_to_codon
from serpent.convert.genetic_code import (
	GENETIC_CODE,
	genetic_code,
	genetic_code_inverse,
)
from serpent.fun import inverse_od, map_array, str_join

aa_tables = [*GENETIC_CODE.keys()]

# TODO: Handle degenerate data
degenerate_amino = OrderedDict({
	"RAY": "B",  # D | N
	"MUH": "J",  # L | I
	"SAR": "Z",  # E | Q
	"NNN": "X",  # Any
	"TAY": "*",  # Termination (stop codon)
	"ZZZ": "-",  # Gap of indeterminate length
})

degenerate_codon = inverse_od(degenerate_amino)


def codon_to_amino(codon, table=1):
	return degenerate_amino.get(codon) or genetic_code[table][codon]

# ruff: noqa: F601
# multi-value-repeated-key-literal
# TODO: Use degenerate nucleotide codes
def amino_to_codon(amino, table=1):
	return degenerate_codon.get(amino) or genetic_code_inverse[table][amino]


def amino_to_num(amino: str, table: int=1) -> int:
	"""Decode an amino acid IUPAC string into a number between 0 and 63."""
	try:
		codon = amino_to_codon(amino, table)
	except KeyError:
		return 0  # TODO Handle degenerate data
	return codon_to_num(codon)


def num_to_amino(number: int, table: int=1) -> str:
	"""Encode a number between 0 and 63 into an amino acid IUPAC string."""
	codon = num_to_codon(number)
	return codon_to_amino(codon, table)


aminos = OrderedDict([
	(num, amino) for (num, amino)
	in enumerate(genetic_code_inverse[1])
])

aminos_inverse = inverse_od(aminos)


def decode_aminos(
	dna,
	table=1
):
	"""Decode amino acid data.

	Uses given nucleotide base order and genetic code translation table.
	"""
	# TODO Use tables?
	table
	# return map_array(lambda d: amino_to_num(d, table), dna)
	return map_array(lambda a: aminos_inverse.get(a, 0), dna)


def aminos_for_table(table: int=1):
	return str_join(OrderedDict([
		(amino, index)
		for index, amino
		in enumerate(reversed(GENETIC_CODE[table]))
	]).keys())

from __future__ import annotations

from collections import OrderedDict

from serpent.fun import inverse_od, map_array

from .codon import CODONS_LEN, codon_to_num, codons_array, num_to_codon
from .genetic_code import GENETIC_CODE

CODONS_NCBI = codons_array('TCAG')
TRANSL_TABLE = 1


def create_genetic_table(table: str) -> OrderedDict[str, str]:
	assert len(table) == CODONS_LEN, 'Table should have 64 characters.'
	return OrderedDict(reversed(list(zip(
		CODONS_NCBI,
		table,
	))))


genetic_code = OrderedDict([
	(i, create_genetic_table(table))
	for i, table in GENETIC_CODE.items()
])

genetic_code_inverse = OrderedDict([
	(i, inverse_od(table))
	for i, table in genetic_code.items()
])


# TODO: Handle degenerate data
degenerate_amino = OrderedDict({
	"RAY": "B",  # D | N
	"MUH": "J",  # L | I
	"SAR": "Z",  # E | Q
	"NNN": "X",  # Any
	"ZZZ": "-",  # Gap of indeterminate length
})

degenerate_codon = inverse_od(degenerate_amino)


def codon_to_amino(codon, table=1):
	return degenerate_amino.get(codon) or genetic_code[table][codon]

# ruff: noqa: F601 # multi-value-repeated-key-literal
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


def decode_aminos(dna, table=1):
	"""Decode amino acid data.

	Uses given nucleotide base order and genetic code translation table.
	"""
	return map_array(lambda d: amino_to_num(d, table), dna)

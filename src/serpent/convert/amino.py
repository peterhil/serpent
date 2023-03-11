from __future__ import annotations

from collections import OrderedDict

from serpent.fun import inverse_od

from .codon import codon_to_num, codons_array, num_to_codon

# Condensed translation table for the Standard Genetic Code
#
# From the NCBI Taxonomy webpage (which also has 25 alternative tables):
# https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1

AMINO_CODE = {
	1:       "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	# Start: "---M------**--*----M---------------M----------------------------"
	# Base1: "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
	# Base2: "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
	# Base3: "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
}

CODONS_NCBI = codons_array('TCAG')

codon_to_amino = OrderedDict(reversed(list(zip(CODONS_NCBI, AMINO_CODE[1]))))

# TODO: Handle degenerate data
# codon_to_amino.update({
# 	"RAY": "B",  # D | N
# 	"MUH": "J",  # L | I
# 	"SAR": "Z",  # E | Q
# 	"NNN": "X",  # Any
# 	"ZZZ": "-",  # Gap of indeterminate length
# })

# ruff: noqa: F601 multi-value-repeated-key-literal
# TODO: Use degenerate nucleotide codes
amino_to_codon = inverse_od(codon_to_amino)


def amino_to_num(amino: str) -> int:
	"""Decode an amino acid IUPAC string into a number between 0 and 63."""
	try:
		codon = amino_to_codon[amino]
	except KeyError:
		return 0  # TODO Handle degenerate data
	return codon_to_num(codon)


def num_to_amino(number: int) -> str:
	"""Encode a number between 0 and 63 into an amino acid IUPAC string."""
	codon = num_to_codon(number)
	return codon_to_amino[codon]

"""DNA and codons data handling."""
from __future__ import annotations

import itertools as itr
import re

import numpy as np
from more_itertools import grouper

from serpent.amino import aminos, aminos_inverse
from serpent.digit import digits_to_number
from serpent.fasta import AMINO, BASE
from serpent.fun import map_array, str_join

bases = {
	"A": 0b00,
	"C": 0b01,
	"G": 0b10,
	"T": 0b11,
	"U": 0b11,
}


bases_inverse = {
	0: "A",
	1: "C",
	2: "G",
	3: "T"
}


def decode(dna, amino=False):
	"""Return codons or amino acids from DNA decoded into numbers 0..63."""
	# TODO Handle degenerate DNA data properly
	dna = clean_non_dna(dna, amino)  # TODO Handle degenerate data better
	if amino:
		return map_array(decode_amino, dna)
	else:
		codons = get_codons(dna)
		return map_array(decode_codon, codons)


def decode_amino(amino: str) -> int:
	"""Decode an amino acid IUPAC string into a number between 0 and 63."""
	return aminos[amino]


def encode_amino(code: int) -> str:
	"""Decode an amino acid IUPAC string into a number between 0 and 63."""
	return aminos_inverse.get(code, '')


def decode_codon(codon: str) -> int:
	"""Decode a codon string into a a number between 0 and 63."""
	result = 0
	for num, char in enumerate(reversed(codon)):
		result += bases[char] << num * 2

	return result


def clean_non_dna(data, amino=False):
	"""Clean up non DNA or RNA data. Warns if there are residual characters."""
	CODES = AMINO if amino else BASE
	cleaned = str_join(re.sub(fr"[^{CODES}]{6,}", "", data).split("\n"))
	residual = str_join(re.findall(fr"[^\n{CODES}]", data))

	if len(residual) > 0:
		# TODO Use logger.warn with warnings.warn?
		print("Residual characters:", residual)

	return cleaned


def get_codons(data, fill="A"):
	"""Get codons from data as Numpy array."""
	codons_list = list(grouper(data, 3, incomplete="fill", fillvalue=fill))
	codons = map_array(str_join, codons_list, dtype="U3")

	return codons


def codon_sequences(decoded, n=4, fill=0):
	"""Chunk data into length N sequences of codons.

	Count the occurences of different kmers as numbers between 0..64**n.
	Return index and counts.
	"""
	sequences = list(grouper(decoded, n, incomplete="fill", fillvalue=fill))
	numbers = np.apply_along_axis(digits_to_number, 1, sequences)

	return numbers


def codons():
	"""Sequence of the 64 codons using order ACGT."""
	return np.fromiter(map(str_join, itr.product('ACGT', repeat=3)), dtype='<U3')


def oligopeptides(length):
	"""Sequence of all oligopeptide combinations of requested length."""
	oligos = map(str_join, itr.product(codons(), repeat=length))

	return np.fromiter(oligos, dtype=f'<U{3*length}')

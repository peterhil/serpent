"""DNA and codons data handling."""
from __future__ import annotations

import itertools as itr
import re

import numpy as np
from more_itertools import grouper

from serpent.config import BASE_ORDER
from serpent.convert.amino import amino_to_num
from serpent.convert.digits import digits_to_num
from serpent.convert.nucleotide import nt_to_num
from serpent.fasta import AMINO, BASE
from serpent.fun import map_array, str_join


def decode(dna, amino=False):
	"""Return codons or amino acids from DNA decoded into numbers 0..63."""
	# TODO Handle degenerate DNA data properly
	dna = clean_non_dna(dna, amino)  # TODO Handle degenerate data better
	if amino:
		return map_array(amino_to_num, dna)
	else:
		codons = get_codons(dna)
		return map_array(decode_codon, codons)


def decode_codon(codon: str) -> int:
	"""Decode a codon string into a a number between 0 and 63."""
	result = 0
	for num, char in enumerate(reversed(codon)):
		result += nt_to_num[char] << num * 2

	return result


def clean_non_dna(data, amino=False):
	"""Clean up non DNA or RNA data. Warns if there are residual characters."""
	# TODO Convert RNA data into DNA, so everything can be handled in base 4 or
	# base 64, and convert back if necessary wehen printing.
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

	Count the occurences of different k-mers as numbers between 0..64**n.
	Return index and counts.
	"""
	sequences = list(grouper(decoded, n, incomplete="fill", fillvalue=fill))
	numbers = np.apply_along_axis(digits_to_num, 1, sequences)

	return numbers


def codons(bases=BASE_ORDER):
	"""Sequence of the 64 codons using the given order."""
	return np.fromiter(map(str_join, itr.product(bases, repeat=3)), dtype='<U3')


def oligopeptides(length):
	"""Sequence of all oligopeptide combinations of requested length.

	Wikipedia: https://en.wikipedia.org/wiki/Oligopeptide
	"""
	oligos = map(str_join, itr.product(codons(), repeat=length))

	return np.fromiter(oligos, dtype=f'<U{3*length}')

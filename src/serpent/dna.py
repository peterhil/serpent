"""DNA and codons data handling."""
from __future__ import annotations

import itertools as itr
import re
from collections.abc import Sequence

import numpy as np
from more_itertools import grouper

from serpent.convert.amino import decode_aminos
from serpent.convert.base64 import num_to_base64
from serpent.convert.codon import codon_to_num, codons_array, num_to_codon
from serpent.convert.degenerate import decode_degenerate
from serpent.convert.digits import digits_to_num
from serpent.fasta import AMINO, BASE, DEGENERATE
from serpent.fun import map_array, str_join
from serpent.settings import BASE_ORDER


def encode(decoded: Sequence[int], fmt: str = 'base64'):
	"""Encode decoded data into base64 or codon format."""
	if fmt in ['b', 'base64']:
		encoded = (num_to_base64.get(num, " ") for num in decoded)
	elif fmt in ['c', 'codon']:
		encoded = (num_to_codon(num) for num in decoded)
	else:
		raise ValueError('Unknown format: ' + fmt)

	return encoded


def file_extension(fmt: str = 'base64'):
	"""Get file extension for encoded data."""
	if fmt in ['b', 'base64']:
		extension = 'ser64'
	elif fmt in ['c', 'codon']:
		extension = 'codon.fasta'
	else:
		raise ValueError('Unknown format: ' + fmt)

	return extension


def decode(dna, amino=False, table=1, degen=False):
	"""Return codons or amino acids from DNA decoded into numbers 0..63.

	Warns if there are residual characters.
	"""
	# TODO: Handle degenerate amino acid data properly
	# TODO: Check data against amino option and warn if used incorrectly?
	[dna, residual] = clean_non_dna(dna, amino, degen)
	if len(residual) > 0:
		# TODO Use logger.warn with warnings.warn?
		print("Residual characters:", residual)

	if amino:
		return decode_aminos(dna, table)
	else:
		codons = get_codons(dna)
		if degen:
			return decode_degenerate(codons)
		else:
			return map_array(codon_to_num, codons)


def clean_non_dna(data, amino=False, degen=False):
	"""Clean up non DNA or RNA data."""
	# TODO Convert RNA data into DNA, so everything can be handled in base 4 or
	# base 64, and convert back when printing if necessary.
	CODES = AMINO if amino else BASE
	if degen and not amino:
		CODES += DEGENERATE
	cleaned = str_join(re.sub(fr"[^{CODES}]{6,}", "", data).split("\n"))
	residual = str_join(re.findall(fr"[^\n{CODES}]", data))

	return [cleaned, residual]


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


def oligopeptides(length, bases=BASE_ORDER):
	"""Sequence of all oligopeptide combinations of requested length.

	Wikipedia: https://en.wikipedia.org/wiki/Oligopeptide
	"""
	oligos = map(str_join, itr.product(codons_array(bases), repeat=length))

	return np.fromiter(oligos, dtype=f'<U{3*length}')

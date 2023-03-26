"""DNA and codons data handling."""
from __future__ import annotations

import sys
from collections.abc import Iterable, Iterator

import numpy as np
from more_itertools import grouper, partition
from numpy.typing import NDArray

from serpent.convert.amino import decode_aminos
from serpent.convert.base64 import num_to_base64
from serpent.convert.codon import codon_to_num, num_to_codon
from serpent.convert.degenerate import degen_to_num
from serpent.convert.digits import digits_to_num
from serpent.fasta import AMINO, BASE, DEGENERATE, RE_WHITESPACE
from serpent.fun import str_join
from serpent.io import err


def encode(decoded: Iterable[int], fmt: str = 'base64') -> Iterable[str]:
	"""Encode decoded data into base64 or codon format."""
	if fmt in ['b', 'base64']:
		encoded = (num_to_base64.get(num, " ") for num in decoded)
	elif fmt in ['c', 'codon']:
		encoded = (num_to_codon(num) for num in decoded)
	else:
		raise ValueError('Unknown format: ' + fmt)

	return encoded


def decode(
	dna: Iterable,
	amino: bool=False,
	table: int=1,
	degen: bool=False
) -> NDArray:
	"""Decode codons or amino acids into array from DNA into numbers.

	The scale of numbers is:
	- degen = False: 0..63     (dtype=np.int8)
	- degen = True:  0...4095  (dtype=np.int32)
	"""
	decoded = decode_iter(dna, amino, table, degen)
	array = _decoded_array(decoded, degen)

	return array


def _decoded_array(decoded: Iterable, degen: bool=False) -> NDArray:
	"""Return iterable decoded data as Numpy array."""
	dtype = np.int32 if degen else np.int8
	decoded = np.fromiter(decoded, dtype)

	return decoded


def decode_iter(
	dna: Iterable,
	amino: bool=False,
	table: int=1,
	degen: bool=False
) -> Iterator:
	"""Decode codons or amino acids iteratively from DNA into numbers 0..63.

	Warns if there are residual characters.
	"""
	# TODO: Handle degenerate amino acid data properly
	# TODO: Check data against amino option and warn if used incorrectly?
	[dna, residual] = clean_non_dna(dna, amino, degen)
	if len(residual) > 0:
		err(f'Residual characters: {residual}')
		if not degen:
			err('Try again with the --degen / -g option.')
			sys.exit(1)

	if amino:
		return decode_aminos(dna, table)
	else:
		codons = get_codons_iter(dna)
		if degen:
			# TODO Handle degenerate data as individual symbols in base 16 instead of codons
			# and convert scale back with np.log2?
			return map(degen_to_num, codons)
		else:
			return map(codon_to_num, codons)


def clean_non_dna(
	data: Iterable, amino: bool=False, degen: bool=False,
) -> tuple[str, str]:
	"""Clean up non DNA or RNA data."""
	# TODO Handle non-coding DNA marked with lowercase symbols.
	# TODO Convert RNA data into DNA, so everything can be handled in base 4 or
	# base 64, and convert back when printing if necessary.
	CODES = AMINO if amino else BASE
	if degen and not amino:
		CODES += DEGENERATE

	[residual, cleaned] = partition(lambda c: c in CODES, data)
	# Filter out whitespace from residual
	[residual, _] = partition(RE_WHITESPACE.match, residual)

	return (str_join([*cleaned]), str_join([*residual]))


def get_codons(data, fill="A"):
	"""Get codons from data as Numpy array."""
	codons = np.fromiter(get_codons_iter(data, fill), dtype="U3")

	return codons


def get_codons_iter(data, fill="A"):
	"""Get codons from data iteratively."""
	codons = grouper(data, 3, incomplete="fill", fillvalue=fill)

	return map(str_join, codons)


def codon_sequences(decoded, n=4, fill=0):
	"""Reinterpret codon data as oligopeptide sequences.

	Split the decoded codon data into sequences of length n and
	intepret the sequences as base 64 numbers.

	Resulting data will have the range of 0..64**n.
	"""
	sequences = list(grouper(decoded, n, incomplete="fill", fillvalue=fill))
	numbers = np.apply_along_axis(digits_to_num, 1, sequences)

	return numbers

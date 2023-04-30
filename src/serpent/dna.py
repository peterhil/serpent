"""DNA and codons data handling."""
from __future__ import annotations

import sys
from collections.abc import Callable, Iterable, Iterator, Sequence

import more_itertools as mit
import numpy as np
from numpy.typing import NDArray

from serpent.convert.amino import codon_to_amino, decode_aminos
from serpent.convert.base64 import num_to_base64
from serpent.convert.codon import codon_to_num, num_to_codon
from serpent.convert.degenerate import dnt_to_num
from serpent.convert.digits import change_base
from serpent.fasta import (
	AMINO,
	AMINO_DEGENERATE,
	BASE,
	BASE_NONCODING,
	DEGENERATE,
	RE_WHITESPACE,
	FastaToken,
	descriptions_and_data,
)
from serpent.fun import str_join
from serpent.io import err


def encode(decoded: Iterable[int], fmt: str = 'base64') -> Iterable[str]:
	"""Encode decoded data into base64 or codon format."""
	if fmt in ['b', 'base64']:
		encoded = (num_to_base64.get(num) for num in decoded)
	elif fmt in ['c', 'codon']:
		encoded = (num_to_codon(num) for num in decoded)
	else:
		raise ValueError('Unknown format: ' + fmt)

	return encoded


def to_amino(
	dna: Iterable,
	amino: bool=False,
	table: int=1,
	degen: bool=False
):
	"""Encode data as amino acids."""
	[cleaned, _] = clean_non_dna(dna, amino, degen)

	if amino:
		# Just clean the data
		yield from cleaned
	else:
		# Convert from codons
		if degen:
			err_msg ="Can't map degenerate data into amino acids."
			raise NotImplementedError(err_msg)
		codons = get_codons_iter(cleaned)
		yield from (codon_to_amino(c, table) for c in codons)


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

	if amino:
		return decode_aminos(dna, table)

	if degen:
		# Handle degenerate data as individual symbols in base 16
		# and combine to degenerate codons
		b16 = map(dnt_to_num, dna)
		return change_base(b16, base=16, n=3)
	else:
		codons = get_codons_iter(dna)
		return map(codon_to_num, codons)


def decode_seq(
	seq: Sequence[FastaToken],
	amino: bool=False,
	table: int=1,
	degen: bool=False,
	decoder: Callable=decode,
):
	"""Decode FASTA token sequence."""
	[descriptions, data] = descriptions_and_data(seq)
	decoded = decoder(data, amino, table, degen)

	return decoded, descriptions


def cleaned_and_residual(
	data: Iterable, amino: bool=False, degen: bool=False,
) -> tuple[Iterable[str], Iterable[str]]:
	"""Get cleaned and residual data from DNA or protein sequence."""
	# TODO Handle non-coding DNA marked with lowercase symbols.
	# TODO Convert RNA data into DNA, so everything can be handled in base 4 or
	# base 64, and convert back when printing if necessary.
	CODES = AMINO if amino else BASE + BASE_NONCODING
	if degen:
		if amino:
			CODES += AMINO_DEGENERATE
		else:
			CODES += DEGENERATE

	(residual, cleaned) = mit.partition(lambda c: c in CODES, data)
	# Filter out whitespace from residual
	(residual, _) = mit.partition(RE_WHITESPACE.match, residual)

	return (cleaned, residual)


def clean_non_dna(
	data: Iterable, amino: bool=False, degen: bool=False,
) -> tuple[str, str]:
	"""Clean up non DNA or RNA data."""
	(cleaned, residual) = map(mit.peekable, cleaned_and_residual(data, amino, degen))

	if residual.peek(''):
		err(f'Residual characters: {str_join([*residual])}')
		if not degen:
			err('Try again with the --degen / -g option.')
			sys.exit(1)

	return (cleaned, residual)


def get_codons(data, fill="A"):
	"""Get codons from data as Numpy array."""
	codons = np.fromiter(get_codons_iter(data, fill), dtype="U3")

	return codons


def get_codons_iter(data, fill="A"):
	"""Get codons from data iteratively."""
	codons = mit.grouper(data, 3, incomplete="fill", fillvalue=fill)

	return map(str_join, codons)


def codon_sequences(decoded, n=4, fill=0):
	"""Reinterpret codon data as oligopeptide sequences.

	Split the decoded codon data into sequences of length n and
	intepret the sequences as base 64 numbers.

	Resulting data will have the range of 0..64**n.
	"""
	return change_base(decoded, base=64, n=n, fill=fill)

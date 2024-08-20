"""Convert data between numbers and encodings."""
# flake8: noqa: F401

from __future__ import annotations

from collections.abc import Callable, Sequence

from serpent.fun import str_join

from .amino import amino_to_num, num_to_amino
from .digits import digits_to_num, num_to_digits
from .nucleotide import nt_to_num, num_to_nt


def convert(func: Callable, seq: Sequence):
	"""Convert a sequence with a translation function to different format.

	For example:
	>>> from .amino import amino_to_codon
	>>> convert(amino_to_codon, 'WPRPQIPP')
	'TGGCCTCGTCCTCAAATTCCTCCT'
	"""
	return str_join([*map(func, seq)])

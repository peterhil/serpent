from __future__ import annotations

from collections import OrderedDict

complement_dnt = OrderedDict({
	'Z': 'Z',
	'G': 'C',
	'A': 'T',
	'C': 'G',
	'T': 'A',
	'R': 'Y',
	'S': 'S',
	'K': 'M',
	'M': 'K',
	'W': 'W',
	'Y': 'R',
	'V': 'B',
	'D': 'H',
	'B': 'V',
	'H': 'D',
	'N': 'N',
})


def complement(seq):
	"""Complement a nucleotide sequence."""
	return (complement_dnt.get(dnt) for dnt in seq)

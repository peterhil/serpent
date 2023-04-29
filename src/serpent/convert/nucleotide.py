"""Convert data between nucleotide bases and numbers."""

from __future__ import annotations

from collections import OrderedDict

import more_itertools as mit

from serpent.fun import inverse_od
from serpent.settings import BASE_ORDER

__all__ = [
	'nt_to_num',
	'num_to_nt',
]

# TODO: Mapping order and numbering could be distinct if using a gray code like
# scheme. For example: With GACT/GACU order the numbering could be 2013, which
# would be the same as ACGT in linear ordering.
#
# ACGT bases:
# A 00 0
# C 01 1
# G 10 2
# T 11 3
#
# GACT bases:
# G 10 2
# A 00 0
# C 01 1
# T 11 3

num_to_nt = OrderedDict(enumerate(BASE_ORDER))
nt_to_num = inverse_od(num_to_nt)


def split_nucleotides(
	codons,
	# Hard coded table 11, TODO Use genetic code tables
	start='GTG ATG ATA ATC ATT CTG TTG',
	stop='TGA TAG TAA',
	split='n'):
	"""Split nucleotide acid sequence by start and stop codons."""
	# TODO Get all start and stop codons from the genetic code tables!
	start = start.split(' ')
	stop = stop.split(' ')

	def start_or_stop(a, b):
		# 'sentences': Split when a start codon follows a stop codon
		if split == 'n':
			return (a in stop and b in start)
		# 'words': Split on start or stop codons
		elif split == 'r':
			return (a in stop or b in start) and a != b

	codons = mit.split_when(codons, start_or_stop)
	return codons

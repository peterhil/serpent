from __future__ import annotations

import more_itertools as mit


def start_or_stop(a, b, start, stop, split='n'):
	# 'sentences': Split when a start codon follows a stop codon
	if split == 'n':
		return (a in stop and b in start)
	# 'words': Split on start or stop codons
	elif split == 'r':
		return (a in stop or b in start) and a != b


def split_aminos(aminos, start='M', stop='*', split='r'):
	"""Split amino acid sequence by start and stop codons."""
	# TODO Get all start and stop codons from the genetic code tables!
	def splitter(a, b):
		return start_or_stop(a, b, start, stop, split)
	aminos = mit.split_when(aminos, splitter)

	yield from aminos


def split_nucleotides(
	codons,
	# Hard coded table 11, TODO Use genetic code tables
	start='GTG ATG ATA ATC ATT CTG TTG',
	stop='TGA TAG TAA',
	split='n'
):
	"""Split nucleotide sequence by start and stop codons."""
	start = start.split(' ')
	stop = stop.split(' ')

	def splitter(a, b):
		return start_or_stop(a, b, start, stop, split)
	codons = mit.split_when(codons, splitter)

	yield from codons


def split_encoded(encoded, fmt, split='n'):
	if fmt in ['a', 'amino']:
		regions = split_aminos(encoded, split=split)
	else:
		regions = split_nucleotides(encoded, split=split)

	yield from regions

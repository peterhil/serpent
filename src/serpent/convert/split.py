from __future__ import annotations

import more_itertools as mit


def split_aminos(aminos, start='M', stop='*', split='r'):
	"""Split amino acid sequence by start and stop codons."""
	# TODO Get all start and stop codons from the genetic code tables!
	def start_or_stop(a, b):
		# 'sentences': Split when a start codon follows a stop codon
		if split == 'n':
			return (a in stop and b in start)
		# 'words': Split on start or stop codons
		elif split == 'r':
			return (a in stop or b in start) and a != b
	aminos = mit.split_when(aminos, start_or_stop)

	return aminos


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

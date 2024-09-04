from __future__ import annotations

import more_itertools as mit

from serpent.convert.genetic_code import STANDARD_TABLES


def splitter_for(start, stop, split='n'):
	started = False

	def splitter(a, b):
		if split == 'n':
			# sentences: Split when a start codon follows a stop codon
			return (a in stop and b in start)

		elif split == 'r':
			# words: Split on start or stop codons
			return (a in stop or b in start) and a != b

		elif split == 'f':
			# frames: Split by open reading frames
			nonlocal started
			if not started and b in start:
				started = True
				return True
			elif started and b in stop:
				started = False
				return True
			else:
				return False
		else:
			err_msg = 'Unknown split type'
			raise ValueError(err_msg)

	return splitter


def split_aminos(aminos, start='M', stop='*', split='f'):
	"""Split amino acid sequence by start and stop codons."""
	# TODO Get all start and stop codons from the genetic code tables!
	splitter = splitter_for(start, stop, split)
	aminos = mit.split_when(aminos, splitter)

	yield from aminos


def split_nucleotides(
	codons,
	table=1,
	split='f'
):
	"""Split nucleotide sequence by start and stop codons."""
	start = STANDARD_TABLES[table].start
	stop = STANDARD_TABLES[table].stop

	splitter = splitter_for(start, stop, split)
	codons = mit.split_when(codons, splitter)

	yield from codons


def split_encoded(encoded, fmt, table, split='f'):
	if fmt in ['a', 'amino']:
		regions = split_aminos(encoded, split=split)
	else:
		regions = split_nucleotides(encoded, table, split=split)

	yield from regions

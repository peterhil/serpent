from __future__ import annotations

from serpent import dna
from serpent.fasta import descriptions_and_data
from serpent.printing import reflow


def encode_data(data, fmt, amino=False, table=1, degen=False):
	"""Encode raw dna data into other formats."""
	if fmt in ['a', 'amino']:
		encoded = dna.to_amino(data, amino, table, degen)
	elif fmt in ['b', 'base64', 'c', 'codon']:
		decoded = dna.decode_iter(data, amino, table, degen)
		encoded = dna.encode(decoded, fmt)
	else:
		err_msg = f'Unknown format: {fmt}'
		raise ValueError(err_msg)

	yield from encoded



def encode_sequences(
	seqs, width, fmt,
	*, amino=False, table=1, degen=False
):
	for sequence in seqs:
		[descriptions, data] = descriptions_and_data(sequence)
		yield from descriptions

		encoded = encode_data(data, fmt, amino, table, degen)
		lines = reflow(encoded, width)
		yield from lines

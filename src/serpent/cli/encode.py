from __future__ import annotations

from serpent import dna


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

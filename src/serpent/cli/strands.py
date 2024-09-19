from __future__ import annotations

import more_itertools as mit

from serpent import dna
from serpent.fun import str_join


def to_strands(
	nucleotides, width,
	*,
	table=1, degen=False
):
	"""Convert nucleotide data to amino acid strands."""
	strands = mit.stagger(nucleotides, offsets=(0, 1, 2), fillvalue=' ', longest=True)
	chunks = mit.chunked(strands, width * 3)
	for chunk in chunks:
		yield from (
			str_join(dna.to_amino(strand, amino=False, table=table, degen=degen))
			for strand in mit.transpose(chunk)
		)
		yield '-' * width

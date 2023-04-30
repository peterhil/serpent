from __future__ import annotations

from serpent import dna
from serpent.convert.digits import num_to_digits
from serpent.fun import str_join
from serpent.io import echo
from serpent.printing import format_lines
from serpent.stats import count_sorted


def analyse_repeats(decoded, length=4, limit=2, fmt='codon'):
	"""Analyse codon sequence repeats."""
	[index, count] = count_sorted(dna.codon_sequences(decoded, length))
	repeats = index[count >= limit]

	echo("Repeated codon sequences:")
	b64 = (num_to_digits(num, base=64) for num in repeats)
	encoded = (str_join(dna.encode(codon, fmt)) for codon in b64)

	width = 32 if fmt in ['b', 'base64'] else 8
	yield from format_lines(encoded, width)

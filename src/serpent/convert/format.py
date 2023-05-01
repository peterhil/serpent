from __future__ import annotations

from serpent.convert.amino import aminos_for_table
from serpent.convert.codon import codons_array
from serpent.encoding import BASE64


def symbols_for(fmt, table=1):
	if fmt in ['a', 'amino']:
		symbols = aminos_for_table(table)
	elif fmt in ['b', 'base64']:
		symbols = BASE64
	elif fmt in ['c', 'codon']:
		symbols = codons_array()  # Use BASE_ORDER
	else:
		err_msg = f'Uknown format: {fmt}'
		raise NotImplementedError(err_msg)

	return symbols

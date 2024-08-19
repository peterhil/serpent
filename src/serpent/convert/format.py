from __future__ import annotations

from serpent.convert.amino import aminos_for_table
from serpent.convert.codon import codons_array


def symbols_for(fmt, table=1):
	if fmt in ['a', 'amino']:
		symbols = aminos_for_table(table)
	elif fmt in ['c', 'codon']:
		symbols = codons_array()  # Use BASE_ORDER
	else:
		err_msg = f'Uknown format: {fmt}'
		raise NotImplementedError(err_msg)

	return symbols

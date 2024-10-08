from __future__ import annotations

import itertools as itr

import more_itertools as mit
import numpy as np

from serpent.convert.format import symbols_for
from serpent.fun import str_join


def all_symbols_for(fmt, seql, table=1):
	item_size = seql * 3 if fmt in ['c', 'codon'] else seql
	dtype = f'U{item_size}'
	symbols = symbols_for(fmt, table)
	combos = np.fromiter(map(str_join, itr.product(symbols, repeat=seql)), dtype=dtype)

	return combos


def peptides_of_length(encoded, seql):
	"""Return peptides of sequence length from encoded data."""
	yield from (str_join(pep) for pep in mit.chunked(encoded, seql))

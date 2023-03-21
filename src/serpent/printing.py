from __future__ import annotations

from collections import Counter
from collections.abc import Iterable

import numpy as np
from more_itertools import chunked, take

from serpent.fun import str_join
from serpent.settings import COUNT_LIMIT


def format_lines(data, width=80, sep=' '):
	"""Format lines for printing."""
	def fmt(chunk):
		return str_join(chunk, sep)
	lines = [
		f"{i * width}:\t{fmt(chunk)}"
		for i, chunk
		in enumerate(chunked(data, width))
	]

	return lines


def format_counts(counts: Counter, limit=COUNT_LIMIT) -> Iterable[str]:
	yield from (
		f'{item}\t{count}'
		for (item, count)
		in take(limit, counts.most_common())
	)

def format_decoded(decoded):
	uniq = len(np.unique(decoded))
	print(f"Decoded ({len(decoded)}, unique: {uniq}):")
	strings = map(str, iter(decoded))
	lines = format_lines(strings, 32)

	return lines

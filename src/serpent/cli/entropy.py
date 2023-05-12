from __future__ import annotations

import more_itertools as mit

from serpent.math.information import statistics, statistics_header


def format_entropy(data, base=2, seql=None):
	if seql:
		yield statistics_header()
		for chunk in mit.chunked(data, seql):
			yield statistics(chunk, base)
	else:
		yield statistics_header()
		yield statistics(data, base)

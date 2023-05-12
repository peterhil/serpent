from __future__ import annotations

import more_itertools as mit

from serpent.math.information import statistics, statistics_header


def check_step_size(seql, step):
	if step is None:
		step = seql
	if not 0 <= step <= seql:
		err_msg = f'Step size should be between 0 and --seql/-q {seql} (inclusive)'
		raise ValueError(err_msg)

	return step


def format_entropy(data, base=2, seql=None, step=None):
	yield statistics_header()
	if seql:
		step = check_step_size(seql, step)
		for chunk in mit.windowed(data, seql, step=step):
			yield statistics(chunk, base)
	else:
		yield statistics(data, base)

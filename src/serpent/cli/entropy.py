from __future__ import annotations

import matplotlib.pyplot as plt
import more_itertools as mit
import numpy as np

from serpent.math.information import efficiency, statistics, statistics_header


def check_step_size(seql, step):
	if step is None:
		step = seql
	err_msg = f'Step size should be between 1 and --seql/-q {seql} (inclusive)'
	assert 1 <= step <= seql, err_msg

	return step


def format_entropy(data, base=2, seql=None, step=None):
	yield statistics_header()
	if seql:
		step = check_step_size(seql, step)
		for chunk in mit.windowed(data, seql, step=step):
			yield statistics(chunk, base)
	else:
		yield statistics(data, base)


def plot_entropy(data, base=2, seql=None, step=None, color=None):
	if seql is None:
		seql = 256

	step = check_step_size(seql, step)
	eff = np.array(
		[efficiency(w, base) for w in mit.windowed(data, seql, step=step)],
		dtype=np.float64
	)
	plt.plot(eff, color=color)

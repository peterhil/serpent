from __future__ import annotations

import more_itertools as mit

from serpent.cli.entropy import check_step_size
from serpent.math.statistic import symbol_frequencies


def sequence_probabilities(data, symbols, seql=64, step=None):
	step = check_step_size(seql, step)
	seqs = mit.windowed(data, seql, step=step, fillvalue='')

	probabilities = (symbol_frequencies(seq, symbols) for seq in seqs)

	yield from probabilities

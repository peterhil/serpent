from __future__ import annotations

from collections import Counter

import more_itertools as mit
import numpy as np

from serpent.cli.entropy import check_step_size


def symbol_frequencies(seq, symbols):
	counts = Counter(seq)
	tide = np.array([counts.get(symbol, 0) for symbol in symbols])

	total = np.sum(tide)
	probabilities = tide / total

	return probabilities


def tide_sequence(data, symbols, seql=64, step=None):
	step = check_step_size(seql, step)
	seqs = mit.windowed(data, seql, step=step, fillvalue='')

	probabilities = (symbol_frequencies(seq, symbols) for seq in seqs)

	yield from probabilities

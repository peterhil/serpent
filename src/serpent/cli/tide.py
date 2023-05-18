from __future__ import annotations

import itertools as itr
from collections import Counter, deque

import more_itertools as mit
import numpy as np

from serpent.cli.entropy import check_step_size
from serpent.fun import str_join


def sequence_probabilities(seq, symbols):
	counts = Counter(seq)

	tide = np.array([counts.get(symbol, 0) for symbol in symbols])
	total = np.sum(tide)

	return tide / total


def tide_sequence(data, symbols, seql=64, step=None):
	step = check_step_size(seql, step)

	data = itr.chain(iter(data))
	initial = mit.take(seql - step, data)
	rest = mit.windowed(data, step, step=step, fillvalue='')

	pool = deque(initial, seql)
	for wave in rest:
		# print('wave:', str_join(filter(None, wave)))
		for symbol in wave:
			pool.append(symbol)

		seq = str_join(pool)
		# yield seq
		yield sequence_probabilities(seq, symbols)

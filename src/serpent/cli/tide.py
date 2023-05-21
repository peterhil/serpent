from __future__ import annotations

import matplotlib.pyplot as plt
import more_itertools as mit
import numpy as np

from serpent.cli.entropy import check_step_size
from serpent.math.statistic import symbol_frequencies
from serpent.visual.palette import hex_spectrum


def sequence_probabilities(data, symbols, seql=64, step=None):
	step = check_step_size(seql, step)
	seqs = mit.windowed(data, seql, step=step, fillvalue='')

	probabilities = (symbol_frequencies(seq, symbols) for seq in seqs)

	yield from probabilities


def tide_sequence(data, symbols, seql, *, step=None, cumulative=False):
	step = check_step_size(seql, step)

	tides = sequence_probabilities(data, symbols, seql, step)
	tides = np.array([*tides])

	if cumulative:
		tides = np.cumsum(tides, axis=1)

	return tides


def plot_tides(tides, symbols):
	colours = hex_spectrum(len(symbols) + 1, sat=0.75, lightness=245, offset=-15/360)
	for tide, color in zip(tides.T, colours):
		plt.plot(tide, color)

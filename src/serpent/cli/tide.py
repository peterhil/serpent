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


def tide_normalise(tides, total):
	scale = total[-1]
	tides = tides / scale
	total = total / scale

	return tides, total


def tide_total(tides, norm=False, color='#000000'):
	total = np.sum(tides, axis=1)

	if norm:
		[tides, total] = tide_normalise(tides, total)

	plt.plot(total, color=color)

	return total


# ruff: noqa: PLR0913
def tide_sequence(
	data, symbols, seql,
	*,
	step=None,
	slopes=False,
	cumulative=False,
):
	step = check_step_size(seql, step)

	tides = sequence_probabilities(data, symbols, seql, step)
	tides = np.array([*tides])

	tides = np.repeat(tides, step, axis=0)

	if cumulative:
		tides = np.cumsum(tides, axis=1)

	if slopes:
		tides = np.cumsum(tides, axis=0)

	return tides


def plot_tides(tides, symbols, alpha=None):
	colours = hex_spectrum(len(symbols) + 1, sat=0.75, lightness=245, offset=-15/360)
	for tide, colour in zip(tides.T, colours, strict=False):
		color = colour + alpha if alpha else colour
		plt.plot(tide, color)

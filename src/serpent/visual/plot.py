"""Visualise data."""

from __future__ import annotations

import sys

import matplotlib.pyplot as plt
import numpy as np

from serpent import dna
from serpent.mapping.amino_spiral_cube import amino_spiral
from serpent.mathematics import magnitude
from serpent.stats import count_sorted

bin_choices = [
	'base', 'auto', 'fd', 'doane', 'scott', 'stone', 'rice', 'sturges', 'sqrt',
]


def interactive() -> None:
	"""Use interactive mode with pyplot."""
	plt.ion()
	plt.show()


def plot_amino_labels(ax, scale=1, color='black', size=10):
	for symbol, coords in amino_spiral.items():
		[y, x, z] = coords
		xyz = np.array([x, y, z]) * scale
		ax.text(*xyz, symbol, color=color, size=size)


def plot_directions(ax, dirs, projection='3d', title=None, **kwargs):
	if projection == '3d':
		plot_path_3d(dirs, ax, **kwargs)
	else:
		plot_path_2d(dirs, ax, **kwargs)

	if title:
		ax.set_title(title)

	return ax


def plot_histogram(
	data,
	*,
	bins='auto',
	cumulative=False,
	density=False,
	histtype='stepfilled',
	**kwargs
):
	"""Plot histograms using np.histogram and plt.hist.

	size:
	Histogram with N (=size) automatically sized bins.
	The np.histogram bins argument seems to have off-by-one error with linspace.

	histtype:
	step and stepfilled are significantly faster for >1000 bins
	default is bar and barstacked is also an option.
	"""
	hist_kwargs = kwargs.copy()
	hist_kwargs.pop('color')

	hist, bins = np.histogram(data, bins=bins, **hist_kwargs)

	plt.hist(
		bins[:-1],
		bins,
		weights=hist,
		histtype=histtype,
		cumulative=cumulative,
		density=density,
		**kwargs
	)

	return [hist, bins]


def plot_histogram_sized(data, *args, size='base', base=64, multi=1, **kwargs):
	"""Plot histograms with automatically sized bins.

	String size can be:
	- base (prefer multiples of some numeric base)
	- auto (fd or sturges)
	- fd, doane, scott, stone, rice, sturges, sqrt

	It makes sense to use base 64 with DNA data, for example.

	Note! Size can be one of the strings described here:
	https://numpy.org/doc/stable/reference/generated/numpy.histogram_bin_edges.html
	"""
	# Prevent overflow and other errors from histogram_bin_edges with long
	# sequences
	if np.min(data) > sys.maxsize or np.max(data) > sys.maxsize:
		size = max(1024, base * multi)

	if size == 'base':
		end = base ** magnitude(np.max(data))
		step = base * multi
		bins = np.linspace(0, end, step, endpoint=True)
	elif isinstance(size, int):
		bins = np.arange(size + 1)
	else:
		bins = np.histogram_bin_edges(data, bins=size)

	hist, bins = plot_histogram(data, *args, bins=bins, **kwargs)

	return [hist, bins]


def plot_path_3d(dirs, ax, **kwargs):
	[y, x, z] = dirs[::, :3].T
	lines = ax.plot3D(x, y, z, **kwargs)

	return lines


def plot_path_2d(dirs, ax, **kwargs):
	[y, x, z] = dirs[::, :3].T
	lines = ax.plot(x, y, **kwargs)

	return lines


def plot_sequence_counts(decoded, *args, n=4, **kwargs):
	"""Plot codon sequence counts of length N."""
	numbers = dna.codon_sequences(decoded, n)
	[index, count] = count_sorted(numbers)

	size = 64**n
	data = np.zeros(size, dtype=np.uint64)
	data[index] = count

	plt.plot(np.arange(size), data, *args, **kwargs)

	return [count, index]

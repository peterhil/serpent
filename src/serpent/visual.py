from __future__ import annotations

import sys

import matplotlib.pyplot as plt
import numpy as np
from numpy.fft import fft
from PIL import Image

from serpent import dna
from serpent.mathematics import magnitude, normalise
from serpent.padding import pad_to_left
from serpent.stats import count_sorted


def interactive():
	plt.interactive(True)
	plt.show()


def plot_fft(decoded, n=64, *args, **kwargs):
	ft = np.abs(fft(decoded, *args, n=n, norm='ortho', **kwargs))
	plt.plot(ft)

	return ft


def plot_histogram(
	data,
	bins='auto',
	cumulative=False,
	density=False,
	histtype='stepfilled',
	*args, **kwargs
):
	"""
	Plot histograms using np.histogram and plt.hist.

	size:
	Histogram with N (=size) automatically sized bins.
	The np.histogram bins argument seems to have off-by-one error with linspace.

	histtype:
	step and stepfilled are significantly faster for >1000 bins
	default is bar and barstacked is also an option.
	"""
	hist, bins = np.histogram(data, *args, bins=bins, **kwargs)

	plt.hist(
		bins[:-1],
		bins,
		*args,
		weights=hist,
		histtype=histtype,
		cumulative=cumulative,
		density=density,
		**kwargs
	)

	return [hist, bins]


def plot_histogram_sized(data, size='base', base=64, multi=1, *args, **kwargs):
	"""Plot histograms with automatically sized bins.

	String size can be:
	- base (prefer multiples of some numeric base)
	- auto (fd or sturges)
	- fd, doane, scott, stone, rice, sturges, sqrt

	It makes sense to use base 64 with DNA data, for example.

	Note! Size can be one of the strings described here:
	https://numpy.org/doc/stable/reference/generated/numpy.histogram_bin_edges.html
	"""
	# Prevent overflow and other errors from histogram_bin_edges with long sequences
	if np.min(data) > sys.maxsize or np.max(data) > sys.maxsize:
		size = max(1024, base * multi)

	if size == 'base':
		bins = np.linspace(0, base ** magnitude(np.max(data)), base * multi, endpoint=True)
	elif type(size) == 'int':
		bins = np.arange(size + 1)
	else:
		bins=np.histogram_bin_edges(data, bins=size)

	hist, bins = plot_histogram(data, *args, bins=bins, **kwargs)

	return [hist, bins]


def plot_sequence_counts(decoded, n=4, *args, **kwargs):
	numbers = dna.codon_sequences(decoded, n)
	[index, count] = count_sorted(numbers)

	size = 64**n
	data = np.zeros(size, dtype=np.uint64)
	data[index] = count

	plt.plot(np.arange(size), data, *args, **kwargs)

	return [count, index]


def show_image(decoded, width=64, fill=0, mode="RGB"):
	padded = np.array(pad_to_left(decoded, 3 * width, fill))
	norm = normalise(padded)
	channels = len(mode)

	rows = len(norm) / (channels * width)
	height = int(np.floor(rows))

	if channels > 1:
		rgb = norm.reshape(height, width, channels)
	else:
		rgb = norm.reshape(height, width)
	img = Image.fromarray(np.uint8(rgb * 255), mode=mode)
	img.show()

	return img

#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

from PIL import Image
from numpy.fft import fft

from serpent import dna
from serpent.mathematics import normalise
from serpent.padding import pad_to_left
from serpent.stats import count_sorted


def interactive():
	plt.interactive(True)
	plt.show()


def plot_fft(decoded, n=64, *args, **kwargs):
	ft = np.abs(fft(decoded, n=n, norm='ortho', *args, **kwargs))
	plt.plot(ft)

	return ft


def plot_histogram(data, *args, **kwargs):
	hist, bins = np.histogram(data, *args, **kwargs)

	plt.plot(bins[:-1], hist)

	return [hist, bins]


def plot_histogram_sized(data, size, *args, **kwargs):
	"""Histogram with N (=size) automatically sized bins.
	The np.histogram bins argument seems to have off-by-one error with linspace."""
	bins = np.arange(size + 1)
	hist, bins = plot_histogram(data, bins=bins, *args, **kwargs)

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

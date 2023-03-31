"""Visualise data."""

from __future__ import annotations

import sys

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from PIL import Image

from serpent import dna
from serpent.bitmap import height_for, num_to_pixel
from serpent.mathematics import magnitude
from serpent.padding import pad_end
from serpent.stats import count_sorted
from serpent.typing import CodonData

bin_choices = [
	'base', 'auto', 'fd', 'doane', 'scott', 'stone', 'rice', 'sturges', 'sqrt',
]


def interactive() -> None:
	"""Use interactive mode with pyplot."""
	plt.ion()
	plt.show()


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


def plot_sequence_counts(decoded, *args, n=4, **kwargs):
	"""Plot codon sequence counts of length N."""
	numbers = dna.codon_sequences(decoded, n)
	[index, count] = count_sorted(numbers)

	size = 64**n
	data = np.zeros(size, dtype=np.uint64)
	data[index] = count

	plt.plot(np.arange(size), data, *args, **kwargs)

	return [count, index]


def plot_path_3d(dirs, ax, **kwargs):
	[y, x, z] = dirs[::, :3].T
	lines = ax.plot3D(x, y, z, **kwargs)

	return lines


def plot_path_2d(dirs, ax, **kwargs):
	[y, x, z] = dirs[::, :3].T
	lines = ax.plot(x, y, **kwargs)

	return lines


def plot_directions(dirs, projection='3d', title='Serpent paths', **kwargs):
	# fig = plt.figure()
	ax = plt.axes(projection=projection)

	if projection == '3d':
		plot_path_3d(dirs, ax, **kwargs)
	else:
		plot_path_2d(dirs, ax, **kwargs)

	if title:
		ax.set_title(title)

	return ax


# ruff: noqa: PLR0913
def dna_image(
	decoded: CodonData, width=64, fill=0, mode="RGB",
	amino=False, degen=False,
) -> Image.Image:
	"""Get decoded DNA data as full colour image.

	See `dna_image_data` for details.
	"""
	rgb = dna_image_data(
		decoded,
		width=width, fill=fill, mode=mode,
		amino=amino, degen=degen,
	)
	img = Image.fromarray(rgb, mode=mode)

	return img


# ruff: noqa: PLR0913
def dna_image_data(
	decoded: CodonData, width=64, fill=0, mode="RGB",
	amino=False, degen=False,
) -> NDArray[np.uint8]:
	"""Convert decoded DNA data to full colour image.

	The codons are mapped to 64 ** 3 (=262144) RGB colours quite directly,
	so that: G: 1, A: 85, C: 169, T (or U): 253.

	Examples
	--------
	A (bluish) grey band of 4 colours from dark to light:

	>>> codons = dna.decode('GGGGGAGGCAAGAAAAATCCGCCACCCTTATTCTTT')
	>>> codons
	array([ 0,  1,  2, 20, 21, 23, 40, 41, 42, 61, 62, 63], dtype=int8)

	>>> im = dna_image(codons, width=4)
	>>> [*im.getdata()]
	[(1, 5, 9), (81, 85, 93), (161, 165, 169), (245, 249, 253)]
	"""
	# TODO decode data here, so gaps can be accomodated for requested width?
	channels: int = len(mode)

	if mode == 'P':
		fill = 255
		padded = np.array(pad_end(decoded, fill, n=width))
		uint8 = np.uint8(padded)
	else:
		padded = np.array(pad_end(decoded, fill, n=3 * width))
		uint8 = num_to_pixel(padded, amino, degen)

	height = height_for(padded, width, channels)

	if channels > 1:
		rgb = uint8.reshape(height, width, channels)
	else:
		rgb = uint8.reshape(height, width)

	return rgb


def dna_image_seq(
	seq,
	width=None, mode="RGB", fill=0,
	*,
	amino=False, degen=False, table=1,
):
	"""Get DNA data from a single sequence of tokens as full colour image data."""
	[decoded, _] = dna.decode_seq(seq, amino, table, degen)
	rgb = dna_image_data(
		decoded,
		width=width, fill=fill, mode=mode,
		amino=amino, degen=degen,
	)

	return rgb

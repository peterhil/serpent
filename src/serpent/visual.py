"""Visualise data."""

from __future__ import annotations

import sys

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from PIL import Image

from serpent import dna
from serpent.bitmap import height_for, num_to_pixel
from serpent.convert.quad import dna_to_quad, quads_to_rgb
from serpent.fasta import data_and_descriptions
from serpent.fun import str_join
from serpent.mapping.amino_spiral_cube import amino_spiral
from serpent.mathematics import logn, magnitude, normalise
from serpent.padding import pad_end
from serpent.printing import format_quasar
from serpent.settings import DEBUG
from serpent.stats import count_sorted, quasar_pulses
from serpent.typing import CodonData

RGB_MAX = 255

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
		fill = RGB_MAX
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


def dna_quad_image(
	seq,
	length=3, width=64,
	*,
	# amino=False,
	degen=False,
	# table=1,
):
	"""Convert DNA data to image using four base colours (YIQ colour space)."""
	(data, descriptions) = data_and_descriptions(seq)

	# TODO Remove degen option altogether and just use dnt conversion
	quads = dna_to_quad(data, length, degen)
	rgb = quads_to_rgb(quads, degen)

	# Compare this with dna_image_data!
	channels = 3
	fill = (0, 0, 0)

	padded = np.array(pad_end(rgb, fill, n=channels * width))
	uint8 = np.uint8(padded)

	height = height_for(padded, width, channels)
	pixels = uint8.reshape(height, width, channels)

	return pixels



def dna_quasar_seq(
	seq,
	amino=False, degen=False, table=1,
	cumulative=False, log=False, mod=0, test=False, key=None
):
	[aminos, descriptions] = dna.decode_seq(seq, amino, table, degen, dna.to_amino)
	pulses, height, scale = quasar_pulses(aminos, cumulative=cumulative, key=key)

	print(str_join(descriptions, '\n'))
	print(format_quasar(pulses.keys())[0])  # Print symbols

	rgb = pulses_to_rgb(pulses, scale, mod=mod, log=log, test=test)

	if DEBUG:
		print(rgb)
	print(f'scale: {scale}')

	return rgb


def pulses_to_rgb(pulses, scale, mod=0, log=False, test=False) -> NDArray[np.uint8]:
	"""Convert pulse repetition data to RGB values."""
	# yield from format_quasar_pulses(pulses, height)
	# OR same data as Numpy array:
	arr = np.vstack([*pulses.values()]).T

	if test:
		arr = np.arange(np.product(arr.shape)).reshape(arr.shape)

	if mod != 0:
		if not log:
			err_msg = f'Modulo needs to be between 1 and {RGB_MAX} when not using the log option'
			assert mod <= RGB_MAX, err_msg
		rgb = arr % mod

		# Enhance image
		normalised = logn(rgb + 1, base=mod) if log else rgb * np.floor(RGB_MAX / mod)
	else:  # Convert to uint8 pixels
		# normalised = normalise(np.log2(arr + 1)) if log else normalise(arr)
		normalised = logn(arr + 1, base=scale) if log else normalise(arr)

	rgb = RGB_MAX * normalised
	return np.uint8(rgb)

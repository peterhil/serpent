"""Bitmap utilities."""

from __future__ import annotations

import colorsys
from collections.abc import Iterable

import numpy as np

from serpent.mathematics import rescale
from serpent.palette import Rgb, amino_colour_map, codon_colour_map


def to_uint8(data, old=64, offset=0):
	"""Rescale data to 8 bit uints.

	>>> to_uint8(np.arange(4), 4, offset=1)
	array([  1,  65, 129, 193], dtype=uint8)
	"""
	return np.uint8(rescale(data, old, 256) + offset)


def height_for(data, width, channels=1):
	"""Return minimum heigth for data to fit in width with N channels."""
	rows: float = data.size / (width * channels)
	height: int = int(np.ceil(rows))

	return height


def num_to_pixel(decoded, amino=False, degen=False):
	"""Convert decoded data into pixel values."""
	# TODO Check degen amino max_value
	max_value = (22 if degen else 22) if amino else (4096 if degen else 64)
	return to_uint8(decoded, max_value, offset=1)  # 1, 5, 9, ..., 249, 253


def num_to_rgb(decoded, amino=False) -> Iterable[Rgb]:
	"""Convert decoded data into pixel values using a palette."""
	colour_map = amino_colour_map if amino else codon_colour_map
	default = (0, 0, 0)  # Black
	pixels = (colour_map.get(n, default) for n in decoded)

	return pixels


def decoded_to_pixels(
	decoded, mode, amino=False, degen=False
) -> Iterable[Rgb] | Iterable[int]:
	"""Convert decoded data into pixel values."""
	if mode == 'P':
		pixels = num_to_rgb(decoded, amino)
	else:
		pixels = num_to_pixel(decoded, amino, degen)

	return pixels


def yiq_to_rgb(yiq):
	"""Convert YIQ pixel values (+-1) to RGB pixels (0..255)."""
	return np.uint8(np.apply_along_axis(lambda c: colorsys.yiq_to_rgb(*c), -1, yiq) * 255)


def yiq_test_image(size: int=8, ymax: float=0.75, imax: float=0.75, qmax: float=0.75):
	"""YIQ colour space is used in the NTSC video format.

	It has Luminance (Y) and chrominance (range: -1...1) channels
	with orange-purple (I) and purple-green (Q) axes.

	Ranges:
	Y: 0...1
    I and Q channels: -1...1

	See: https://en.wikipedia.org/wiki/YIQ
	"""
	luma = np.ones((size, size)) * ymax
	i = np.linspace(-1, 1, size)  # blue-orange (I)
	q = np.linspace(-1, 1, size)  # green-purple (Q)

	iv, qv = np.meshgrid(i * imax, q * qmax)
	yiq = np.stack([luma, iv, qv], axis=-1)

	rgb = yiq_to_rgb(yiq)
	# im = Image.fromarray(rgb)
	# im.show()
	return rgb

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Sequence
from colorsys import hsv_to_rgb

import numpy as np
import PIL
from PIL.Image import Image
from PIL.ImagePalette import ImagePalette

Rgb = tuple[int, int, int]


def hue(hue: float=0, sat: float=1.0, lightness: int=255) -> Rgb:
	return hsv_to_rgb(hue, sat, lightness)


def rgb_to_hex(rgb: Rgb) -> str:
	hexa = ''.join([hex(ch)[2:].zfill(2) for ch in rgb])
	return f'#{hexa}'


def spectrum(n=64, sat=1.0, lightness=255, offset=0):
	"""Spectrum of n colours."""
	hues = (np.linspace(0, 1, n, endpoint=False) + offset) % 1
	colours = np.array([hue(v, sat, lightness) for v in hues], dtype=np.uint8)

	return colours


def hex_spectrum(n, **kwargs):
	return [rgb_to_hex(c) for c in spectrum(n, **kwargs)]


def spectrum_layers(n=9, layers=3, start=0.5, *, sat=0.75, offset=0):
	"""Spectrum colours with layers of lightness."""
	# TODO Rethink the amino acid to colour mapping?
	lightness = np.linspace(start, 1, layers, endpoint=True) * 256 - 1
	colours = np.concatenate([
		spectrum(n, sat, light, offset)
		for light in reversed(lightness)
	])

	return colours


def arr_to_palette(arr):
	length = len(arr)
	pad_length = min(768, 768 - 3 * length)
	padding = bytes(pad_length)

	return bytearray(arr.flat) + padding


def colours_to_palette(colours: Sequence[Rgb]) -> ImagePalette:
	"""Convert a sequence of RGB values into PIL ImagePalette."""
	return ImagePalette('RGB', palette=arr_to_palette(colours))


def set_palette(im: Image, palette: ImagePalette):
	im.palette = palette
	im.palette.dirty = True

	return im


def spectrum_layer_colours(amino: bool, degen: bool=False):
	if amino:
		# 27 colours
		layers = 3
		num_colours = 9
		start = 1/2
	else:
		# 64 colours (256 for degen)
		layers = 4 if not degen else 16
		num_colours = 16
		start = 1/2
	colours = spectrum_layers(num_colours, layers, start, offset=-10/360)

	return colours


def spectrum_layer_colours_for(key: str, amino: bool=False):
	colour_map = spectrum_layer_colours(amino)
	colours = OrderedDict([
		[aa, rgb_to_hex(rgb)]
		for aa, rgb in zip(key, colour_map, strict=False)
	])

	return colours


def spectrum_layer_colour_map(amino: bool, degen: bool=False):
	colours = spectrum_layer_colours(amino, degen)
	colour_map = OrderedDict((i, tuple(c)) for i, c in enumerate(colours))

	return colour_map


amino_colour_map = spectrum_layer_colour_map(amino=True)
codon_colour_map = spectrum_layer_colour_map(amino=False)
degen_colour_map = spectrum_layer_colour_map(amino=False, degen=True)


def apply_palette(img: Image, amino: bool=False) -> Image:
	colours = spectrum_layer_colours(amino)
	palette = colours_to_palette(colours)
	return set_palette(img, palette)


if __name__ == '__main__':
	# Test data
	data = np.arange(27, dtype=np.uint8).reshape(3, 9)

	spectra = PIL.Image.fromarray(data, mode='P')
	spectra = apply_palette(spectra, amino=True)
	spectra.show()

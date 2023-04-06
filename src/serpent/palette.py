from __future__ import annotations

from colorsys import hsv_to_rgb

import numpy as np
import PIL
from PIL.Image import Image
from PIL.ImagePalette import ImagePalette

saturation = 1.0
lightness = 255


def hue(hue: float=0, sat: float=1.0, lightness: int=255) -> tuple[int, int, int]:
	return hsv_to_rgb(hue, sat, lightness)


def spectrum(n=64, sat=1.0, lightness=255, offset=0):
	"""Spectrum of n colours."""
	hues = (np.linspace(0, 1, n, endpoint=False) + offset) % 1
	return np.array([hue(v, sat, lightness) for v in hues], dtype=np.uint8)


def spectrum_palette(n=64, sat=1.0, lightness=255, offset=0):
	"""Spectrum palette suitable for PIL Image instance."""
	colours = spectrum(n, sat, lightness, offset)
	return ImagePalette('RGB', palette=arr_to_palette(colours))


def spectrum_layers_palette(n=9, layers=3, start=0.5, *, sat=0.75, offset=0):
	"""Spectrum palette with layers of lightness.

	Suitable for PIL Image instance.
	"""
	lightness = np.linspace(start, 1, layers, endpoint=True) * 256 - 1
	colours = np.concatenate([
		spectrum(n, sat, light, offset)
		for light in reversed(lightness)
	])
	return ImagePalette('RGB', palette=arr_to_palette(colours))


def arr_to_palette(arr):
	length = len(arr)
	pad_length = min(768, 768 - 3 * length)
	padding = bytes(pad_length)

	return bytearray(arr.flat) + padding


def set_palette(im: Image, palette: ImagePalette):
	im.palette = palette
	im.palette.dirty = True

	return im


def apply_palette(img: Image, amino: bool=False) -> Image:
	# num_colours = 22 if amino else 64
	# set_palette(img, spectrum_palette(num_colours, 0.75))
	if amino:
		# 27 colours
		layers = 3
		num_colours = 9
		start = 0.25
	else:
		# 64 colours
		layers = 4
		num_colours = 16
		start = 0.5
	# TODO Rethink the mapping
	palette = spectrum_layers_palette(num_colours, layers, start, offset=-10/360)
	return set_palette(img, palette)


if __name__ == '__main__':
	# Test data
	data = np.arange(27, dtype=np.uint8).reshape(3, 9)

	spectra = PIL.Image.fromarray(data, mode='P')
	# new_palette = spectrum_palette(256, sat=0.75, offset=0/360)
	new_palette = spectrum_layers_palette(
		n=9, layers=3, start=0.25, sat=0.75, offset=-10/360
	)

	set_palette(spectra, new_palette)

	spectra.show()

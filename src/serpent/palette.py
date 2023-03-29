from __future__ import annotations

from colorsys import hsv_to_rgb

import numpy as np
import PIL
from PIL.Image import Image
from PIL.ImagePalette import ImagePalette

from serpent.combinatorics import spread

saturation = 1.0
lightness = 255


def hue(hue: float=0, sat: float=1.0, lightness: int=255) -> tuple[int, int, int]:
	return hsv_to_rgb(hue, sat, lightness)


def spectrum(n=64, sat=1.0, lightness=255):
	"""Spectrum of n colours."""
	hues = np.linspace(0, 1, n, endpoint=False)
	return np.array([hue(v, sat, lightness) for v in hues], dtype=np.uint8)


def spectrum_palette(n=64, sat=1.0, lightness=255):
	"""Spectrum palette suitable for PIL Image instance."""
	colours = spectrum(n, sat, lightness)
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


if __name__ == '__main__':
	# Test data
	a = np.arange(256 * 3)
	data = (spread(a, 64) % 256).reshape(3 * 16, 16).astype(np.uint8)

	spectra = PIL.Image.fromarray(data.T, mode='P')
	new_palette = spectrum_palette(256)

	set_palette(spectra, new_palette)

	spectra.show()

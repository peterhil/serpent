from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from PIL import Image

from serpent import dna
from serpent.bitmap import RGB_MAX, height_for, num_to_pixel
from serpent.convert.quad import dna_to_quad, quads_to_rgb
from serpent.fasta import descriptions_and_data
from serpent.padding import pad_end
from serpent.palette import apply_palette
from serpent.typing import CodonData


# ruff: noqa: PLR0913
def dna_image_data(
	decoded: CodonData, width=64, fill=0, mode="RGB",
	amino=False, degen=False,
) -> NDArray[np.uint8]:
	"""Convert decoded DNA data to full colour image.

	The codons are mapped to 64 ** 3 (=262144) RGB colours quite directly,
	so that: G: 1, A: 85, C: 169, T (or U): 253.
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
	[descriptions, data] = descriptions_and_data(seq)

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


def dna_image(
	seqs, width=None, mode="RGB",
	*,
	length=1,
	amino=False, degen=False, table=1,
) -> Image.Image:
	"""DNA data as full colour image in various modes."""
	if mode == 'Q':
		rgb = np.vstack([
			dna_quad_image(
				seq, length, width,
				# amino=amino,
				degen=degen,
				# table=table,
			) for seq in seqs
		])
		mode = 'RGB'
	else:
		rgb = np.vstack([
			dna_image_seq(
				seq, width, mode, fill=0,
				amino=amino, degen=degen, table=table)
			for seq in seqs
		])

	img = Image.fromarray(rgb, mode)

	if mode == 'P':
		img = apply_palette(img, amino)

	return img

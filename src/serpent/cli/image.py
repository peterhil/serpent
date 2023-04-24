from __future__ import annotations

import numpy as np
from PIL import Image

from serpent.palette import apply_palette
from serpent.visual import dna_image_seq, dna_quad_image


def dna_image(
	seqs, width=None, mode="RGB",
	*,
	length=1,
	amino=False, degen=False, table=1,
):
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

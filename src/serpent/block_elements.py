"""Unicode Block Element graphics for terminals.

See:
https://en.wikipedia.org/wiki/Block_Elements
"""
from __future__ import annotations

import itertools as itr
from collections.abc import Iterator, Sequence

from more_itertools import chunked, grouper

from serpent import ansi

HALF_BLOCK = '\u2580'


def pixels_to_blocks(pixels: Sequence, width: int) -> Iterator[str]:
	"""Convert a sequence of RGB pixel values to Unicode block element graphics."""
	rgb = grouper(pixels, 3, incomplete='fill', fillvalue=0)
	lines = chunked(rgb, width * 2)  # double the line width
	zero_pixel = (0, 0, 0)

	for line in lines:
		top_and_bottom = chunked(line, width)  # split in half
		columns = itr.zip_longest(*top_and_bottom, fillvalue=zero_pixel)

		blocks = ''
		for column in columns:
			colours = ansi.rgb(*column)
			blocks += colours + HALF_BLOCK

		yield blocks + ansi.RESET

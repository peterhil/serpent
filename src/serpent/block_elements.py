"""Unicode Block Element graphics for terminals.

See:
https://en.wikipedia.org/wiki/Block_Elements
"""
from __future__ import annotations

import itertools as itr
from collections.abc import Iterator, Sequence

import more_itertools as mit

from serpent import ansi

HALF_BLOCK = '\u2580'


def pixels_to_blocks(pixels: Sequence, width: int, *, mode: str='RGB') -> Iterator[str]:
	"""Convert a sequence of RGB pixel values to Unicode block element graphics."""
	if mode == 'RGB':
		rgb = mit.grouper(pixels, 3, incomplete='fill', fillvalue=0)
	else:
		rgb = pixels
	lines = mit.chunked(rgb, width * 2)  # double the line width
	zero_pixel = (0, 0, 0) if mode == 'RGB' else 0

	for line in lines:
		top_and_bottom = mit.chunked(line, width)  # split in half
		columns = itr.zip_longest(*top_and_bottom, fillvalue=zero_pixel)

		blocks = ''
		for column in columns:
			# FIXME Fixes the last line problem, that is caused by zero pixel bg colours
			# getting filtered out, probably because of the spread operators
			values = (column[0], zero_pixel) if len(column) == 1 else column
			colours = ansi.rgb(*values) if mode == 'RGB' else ansi.grey(*values)
			blocks += colours + HALF_BLOCK

		yield blocks + ansi.RESET

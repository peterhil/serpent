"""Flow command module."""

from __future__ import annotations

import itertools as itr

import more_itertools as mit

from serpent import dna
from serpent.visual.bitmap import decoded_to_pixels
from serpent.visual.block_elements import pixels_to_blocks, pixels_to_verbose_blocks


# ruff: noqa: PLR0913 # Too many arguments in function definition
def flow_blocks(
	data, width=64, mode='RGB',
	*,
	amino=False, degen=False, table=1,
	height: int | None=None,
):
	decoded = dna.decode(data, amino, table, degen)
	pixels = decoded_to_pixels(decoded, mode, amino, degen)
	yield from pixels_to_blocks(pixels, width, height=height, mode=mode)


# ruff: noqa: PLR0913 # Too many arguments in function definition
def verbose_flow_blocks(
	data, width=64, mode='RGB', fmt=None,
	*,
	amino=False, degen=False, table=1,
	height: int | None=None,
):
	(data1, data2) = itr.tee(data, 2)

	if not fmt:
		fmt = 'amino' if amino else 'codon'

	if fmt in ['a', 'amino']:
		text = dna.to_amino(data1, amino, table, degen)
	elif fmt in ['c', 'codon']:
		if not amino:
			text = data1
		else:
			decoded = dna.decode(data1, amino, table, degen)
			text = mit.flatten(dna.encode(decoded, fmt=fmt, degen=degen))
	else:
		err_msg = 'Unknown format'
		raise NotImplementedError(err_msg)

	decoded = dna.decode(data2, amino, table, degen)
	pixels = decoded_to_pixels(decoded, mode, amino, degen)
	# TODO Use three rows per codon for 1.5 times more data per screen?
	repeat = 3 * len(mode) if fmt in ['c', 'codon'] else len(mode)

	# TODO Decode data later after splitting to lines and
	# only after that convert to pixels
	yield from pixels_to_verbose_blocks(
		pixels, text, width, height=height, mode=mode, repeat=repeat
	)

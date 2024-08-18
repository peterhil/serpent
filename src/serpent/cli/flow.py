"""Flow command module."""

from __future__ import annotations

from serpent import dna
from serpent.fun import str_join
from serpent.visual.bitmap import decoded_to_pixels
from serpent.visual.block_elements import pixels_to_blocks, pixels_to_verbose_blocks


def flow_blocks(
	data, width=64, mode='RGB',
	*,
	amino=False, degen=False, table=1,
):
	decoded = dna.decode(data, amino, table, degen)
	pixels = decoded_to_pixels(decoded, mode, amino, degen)
	yield from pixels_to_blocks(pixels, width, mode=mode)


def verbose_flow_blocks(
	data, width=64, mode='RGB', fmt=None,
	*,
	amino=False, degen=False, table=1,
):
	decoded = dna.decode(data, amino, table, degen)

	if not fmt:
		fmt = 'amino' if amino else 'codon'

	if fmt in ['a', 'amino']:
		text = dna.to_amino(data, amino, table, degen)
	elif fmt in ['c', 'codon', 'b', 'base64']:
		text = str_join(dna.encode(decoded, fmt=fmt, degen=degen))
	else:
		err_msg = 'Unknown format'
		raise NotImplementedError(err_msg)

	pixels = decoded_to_pixels(decoded, mode, amino, degen)
	# TODO Use three rows per codon for 1.5 times more data per screen?
	repeat = 3 * len(mode) if fmt in ['c', 'codon'] else len(mode)

	yield from pixels_to_verbose_blocks(
		pixels, text, width, mode=mode, repeat=repeat
	)

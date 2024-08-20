"""Quasar command visualises symbol pulse repetition intervals as images."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from serpent import dna
from serpent.fun import str_join
from serpent.io.printing import format_quasar
from serpent.math.basic import logn, normalise
from serpent.math.statistic import quasar_pulses
from serpent.settings import DEBUG
from serpent.visual.bitmap import RGB_MAX


# ruff: noqa: PLR0913
def dna_quasar_seq(
	seq,
	*,
	cumulative=False, log=False, mod=0,
	test=False, key=None,
	amino=False, degen=False, table=1,
):
	[aminos, descriptions] = dna.decode_seq(seq, amino, table, degen, dna.to_amino)
	pulses, height, scale = quasar_pulses(aminos, cumulative=cumulative, key=key)

	print(str_join(descriptions, '\n'))
	print(format_quasar(pulses.keys())[0])  # Print symbols

	rgb = pulses_to_rgb(pulses, scale, mod=mod, log=log, test=test)

	if DEBUG:
		print(rgb)
	print(f'scale: {scale}')

	return rgb


def pulses_to_rgb(pulses, scale, mod=0, log=False, test=False) -> NDArray[np.uint8]:
	"""Convert pulse repetition data to RGB values."""
	# yield from format_quasar_pulses(pulses, height)
	# OR same data as Numpy array:
	arr = np.vstack([*pulses.values()]).T

	if test:
		arr = np.arange(np.prod(arr.shape)).reshape(arr.shape)

	if mod != 0:
		if not log:
			err_msg = f'Modulo needs to be 1...{RGB_MAX}, or use the log option'
			assert mod <= RGB_MAX, err_msg
		rgb = arr % mod

		# Enhance image
		normalised = logn(rgb + 1, base=mod) if log else rgb * np.floor(RGB_MAX / mod)
	else:  # Convert to uint8 pixels
		# normalised = normalise(np.log2(arr + 1)) if log else normalise(arr)
		normalised = logn(arr + 1, base=scale) if log else normalise(arr)

	rgb = RGB_MAX * normalised
	return np.uint8(rgb)

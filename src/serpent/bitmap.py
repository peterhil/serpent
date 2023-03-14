"""Bitmap utilities."""

from __future__ import annotations

import numpy as np

from serpent.mathematics import rescale


def to_uint8(data, old=64, offset=0):
	"""Rescale data to 8 bit uints.

	>>> to_uint8(np.arange(4), 4, offset=1)
	array([  1,  65, 129, 193], dtype=uint8)
	"""
	return np.uint8(rescale(data, old, 256) + offset)


def height_for(data, width, channels=1):
	"""Return minimum heigth for data to fit in width with N channels."""
	rows: float = len(data) / (width * channels)
	height: int = int(np.ceil(rows))

	return height

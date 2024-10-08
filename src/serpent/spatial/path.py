from __future__ import annotations

import itertools as itr

import numpy as np

from serpent.spatial.amino_spiral_cube import amino_spiral

directions = np.array(list(itr.product([-1, 0, 1], repeat=3)), dtype=np.int8)


def amino_to_dir(amino):
	return amino_spiral.get(amino)


def amino_path_3d(aminos):
	"""Make a 3d path with a directions mapping."""
	dirs = np.array(list(map(amino_to_dir, aminos)), dtype=np.int64)
	path = np.cumsum(dirs, axis=0)

	return path

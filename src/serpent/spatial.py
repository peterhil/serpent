from __future__ import annotations

import itertools as itr

import numpy as np

from serpent.mapping.amino_spiral_cube import amino_spiral

directions = np.array(list(itr.product([-1, 0, 1], repeat=3)), dtype=np.int8)


def amino_to_dir(amino):
	return amino_spiral.get(amino)

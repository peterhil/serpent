from __future__ import annotations

import numpy as np

from serpent.cli.entropy import check_step_size
from serpent.convert.quad import dna_to_quad
from serpent.math.basic import normalise


# ruff: noqa: PLR0913 # Too many arguments in function definition
def walk_sequence(
	data,
	seql,
	step=None,
	*,
	degen=False,
	norm=False,
	unit=False,
):
	step = check_step_size(seql, step)
	quads = dna_to_quad(data, length=seql, degen=degen, step=step)
	dirs = np.array(quads).T if unit else np.cumsum(quads, axis=0).T
	path = normalise(dirs) if norm else dirs if unit else dirs * step

	return path

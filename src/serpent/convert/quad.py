"""Quadrature colour conversions."""

from __future__ import annotations

import itertools as itr
from collections import OrderedDict
from collections.abc import Iterable

import more_itertools as mit
import numpy as np

from serpent.bitmap import yiq_to_rgb
from serpent.settings import BASE_ORDER

QUAD_ZERO = np.array((0, 0))

nt_to_quad = OrderedDict(zip(
	BASE_ORDER,
	[*itr.product([1, -1], repeat=2)]
))


def peptide_to_quad(peptide: Iterable[str]):
	if len(peptide) == 0:
		return QUAD_ZERO
	quads = [nt_to_quad.get(nt) for nt in peptide]
	return np.mean(quads, axis=0)


def dna_to_quad(dna: Iterable[str], length=3):
	# TODO Also try moving average and exponential smoothing!
	quads = np.array([peptide_to_quad(peptide) for peptide in mit.chunked(dna, length)])
	return quads


def quads_to_rgb(quads, lightness=0.75):
	# TODO Handle lightness on dna_to_quad or elsewhere
	luma = np.ones(len(quads)) * lightness
	yiq = np.vstack([luma, quads.T]).T

	rgb = yiq_to_rgb(yiq)

	return rgb

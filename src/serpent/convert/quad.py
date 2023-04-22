"""Quadrature colour conversions."""

from __future__ import annotations

import itertools as itr
from collections import OrderedDict
from collections.abc import Iterable

import more_itertools as mit
import numpy as np

from serpent.settings import BASE_ORDER

nt_to_quad = OrderedDict(zip(
	BASE_ORDER,
	[*itr.product([1, -1], repeat=2)]
))


def peptide_to_quad(peptide: Iterable[str]):
	quads = [nt_to_quad.get(nt) for nt in peptide]
	return np.mean(quads, axis=0)


def dna_to_quad(dna: Iterable[str], length=3):
	# TODO Try also moving average and exponential smoothing!
	quads = [peptide_to_quad(peptide) for peptide in mit.chunked(dna, length)]
	return quads

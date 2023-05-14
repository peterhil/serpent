"""Quadrature colour conversions."""

from __future__ import annotations

import itertools as itr
from collections import OrderedDict
from collections.abc import Iterable

import more_itertools as mit
import numpy as np

from serpent.convert.dnt import bits_to_dnt, dnt_to_bits
from serpent.fun import str_join
from serpent.settings import BASE_ORDER
from serpent.visual.bitmap import yiq_to_rgb

QUAD_ZERO = np.array((0, 0))

nt_to_quad = OrderedDict(zip(
	BASE_ORDER,
	[*itr.product([1, -1], repeat=2)]
))


def peptide_to_quad(peptide: Iterable[str]):
	if len(peptide) == 0:
		return QUAD_ZERO
	quads = [nt_to_quad[nt] for nt in peptide]
	return np.mean(quads, axis=0)


def create_dnt_to_quad(bases=BASE_ORDER):
	def quad_mapping(nts: str) -> tuple[str, tuple[float, float, float]]:
		"""Map string (combination) of nongenerate bases into YIQ quad value."""
		summed = np.sum([dnt_to_bits[nt] for nt in nts])  # TODO Use map with operator.or_?
		dnt = bits_to_dnt[summed]
		quad = np.concatenate([[1 - (len(nts) / len(bases))], peptide_to_quad(nts)])
		return (dnt, quad)

	return OrderedDict([
		quad_mapping(quad)
		for quad in mit.powerset(nt_to_quad)
	])


dnt_to_quad = create_dnt_to_quad()


def degen_to_quad(peptide: Iterable[str]):
	if len(peptide) == 0:
		peptide = 'Z'
	quads = [dnt_to_quad[nt] for nt in peptide]
	return np.mean(quads, axis=0)


def dna_to_quad(dna: Iterable[str], length=3, degen=False, step=1):
	convert = degen_to_quad if degen else peptide_to_quad
	# TODO Also try exponential moving average smoothing!
	chunks = [
		str_join(filter(None, chunk))
		for chunk in mit.windowed(dna, length, step=step)
	]
	quads = np.array([convert(peptide) for peptide in chunks])
	return quads


def quads_to_rgb(quads, degen=False):
	if degen:
		yiq = quads
	else:
		# TODO Handle lightness on dna_to_quad or elsewhere
		lightness = 0.75
		luma = np.ones(len(quads)) * lightness
		yiq = np.vstack([luma, quads.T]).T

	rgb = yiq_to_rgb(yiq)

	return rgb

from __future__ import annotations

import numpy as np
from hypothesis import given
from hypothesis import strategies as st
from numpy.testing import assert_array_equal

from serpent.convert.quad import QUAD_ZERO, dna_to_quad, nt_to_quad, peptide_to_quad

a_nucleotide = st.sampled_from('GACT')


def positive_integer(min_value=0, max_value=None):
	return st.integers(min_value=min_value, max_value=max_value)


def test_nt_to_quad():
	assert nt_to_quad.get('G') == (1, 1)
	assert nt_to_quad.get('A') == (1, -1)
	assert nt_to_quad.get('C') == (-1, 1)
	assert nt_to_quad.get('T') == (-1, -1)


@given(st.data())
def test_peptide_to_quad_repeats_are_same_as_one(data):
	nt = data.draw(a_nucleotide)
	length = data.draw(positive_integer(1, 12))
	quadrant = nt_to_quad.get(nt)

	assert_array_equal(peptide_to_quad(nt * length), quadrant)


@given(st.data())
def test_peptide_to_quad_dominant(data):
	dominant = data.draw(a_nucleotide)
	pep = data.draw(a_nucleotide)
	quadrant = nt_to_quad.get(dominant)

	# TODO Find a way to shuffle or draw a weighted sample
	actual = peptide_to_quad(dominant * 2 + pep)

	assert_array_equal(np.sign(actual), quadrant)


def test_peptide_to_quad_empty_is_zero():
	assert_array_equal(peptide_to_quad(''), QUAD_ZERO)


def test_dna_to_quad():
	G = (1, 1)
	assert_array_equal(dna_to_quad('GGG'), [G])

from __future__ import annotations

from numpy.testing import assert_array_equal

from serpent.convert.quad import dna_to_quad, nt_to_quad, peptide_to_quad


def test_nt_to_quad():
	assert nt_to_quad.get('G') == (1, 1)
	assert nt_to_quad.get('A') == (1, -1)
	assert nt_to_quad.get('C') == (-1, 1)
	assert nt_to_quad.get('T') == (-1, -1)


def test_peptide_to_quad():
	G = (1, 1)
	assert_array_equal(peptide_to_quad('GGG'), G)


def test_dna_to_quad():
	G = (1, 1)
	assert_array_equal(dna_to_quad('GGG'), [G])

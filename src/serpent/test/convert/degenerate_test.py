from __future__ import annotations

from numpy.testing import assert_array_equal

from serpent import dna


def test_dna_decode_is_same_for_degenerate():
	data = 'GGGAAACCCTTT'
	scale = 52

	assert_array_equal(
		dna.decode(data, degen=True) / scale,
		dna.decode(data),
	)

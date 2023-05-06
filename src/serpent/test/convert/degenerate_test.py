from __future__ import annotations

from numpy.testing import assert_array_equal

from serpent import dna
from serpent.convert.degenerate import degen_to_amino
from serpent.test.convert.fixed_degen import fixed_degen_to_amino


def test_degen_to_amino():
	expected_map = fixed_degen_to_amino[1]
	for degen in expected_map:
		assert degen_to_amino(degen) == expected_map[degen]


def test_dna_decode_is_same_for_degenerate():
	data = 'GGGAAACCCTTT'
	scale = 52

	assert_array_equal(
		dna.decode(data, degen=True) / scale,
		dna.decode(data),
	)

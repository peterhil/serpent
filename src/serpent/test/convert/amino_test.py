from __future__ import annotations

from serpent.convert.amino import num_to_amino
from serpent.convert.codon import codon_to_num
from serpent.test.convert.fixed_amino import fixed_codon_to_amino


def assert_equal(actual, expected):
	assert actual == expected, f'expected: {expected}, got: {actual}'


def codon_to_amino(codon):
	return num_to_amino(codon_to_num(codon))


def test_codon_to_amino():
	for codon, expected in fixed_codon_to_amino.items():
		assert_equal(
			codon_to_amino(codon),
			expected
		)

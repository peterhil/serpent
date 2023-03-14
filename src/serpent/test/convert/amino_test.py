from __future__ import annotations

from serpent.convert.amino import codon_to_amino
from serpent.test.convert.fixed_amino import fixed_codon_to_amino


def assert_equal(actual, expected):
	assert actual == expected, f'expected: {expected}, got: {actual}'


def test_codon_to_amino():
	for codon, expected in fixed_codon_to_amino.items():
		assert_equal(
			codon_to_amino(codon),
			expected
		)

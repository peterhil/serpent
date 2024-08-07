from __future__ import annotations

import math
from collections import Counter

import pytest

from serpent.math.information import abs_rate, entropy_shannon

abs_rates = [
	(-math.inf, ''),
	(0, 'AA'),
	(0.5, 'GA'),
	(0.792481250360578, 'GAC'),
	(1, 'GACT'),
	(1, 'ATGCCCTAG'),
]

@pytest.mark.parametrize(('expected', 'seq'), abs_rates)
def test_abs_rate(expected, seq):
	base = 4
	assert expected == abs_rate(seq, base) == abs_rate(Counter(seq), base)


bin_entropies = [
	(0.0, ''),
	(0.0, 'AAA'),
	(1.0, 'GAGA'),
	(1.8423709931771084, 'GATTACA'),
	(2.0, 'GACT'),
	(2.197159723424149, 'ATGCCCUGA'),
	(2.25, 'ATGCCUGA'),
]

quad_entropies = [
	(0.0, 'GG'),
	(0.4591479170272448, 'GAG'),
	(0.5, 'GA'),
	(0.792481250360578, 'ATG'),
	(1.0, 'GACT'),
]


@pytest.mark.parametrize(('expected', 'seq'), bin_entropies)
def test_entropy_shannon(expected, seq):
	assert expected == entropy_shannon(seq)


@pytest.mark.parametrize(('expected', 'seq'), bin_entropies)
def test_entropy_shannon_counter(expected, seq):
	assert expected == entropy_shannon(Counter(seq))


@pytest.mark.parametrize(('expected', 'seq'), quad_entropies)
def test_entropy_shannon_quad_base(expected, seq):
	assert expected == entropy_shannon(seq, base=4)


@pytest.mark.parametrize(('expected', 'seq'), quad_entropies)
def test_entropy_shannon_counter_quad_base(expected, seq):
	assert expected == entropy_shannon(Counter(seq), base=4)

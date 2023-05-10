from __future__ import annotations

from collections import Counter

import pytest

from serpent.math.information import entropy

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
def test_entropy(expected, seq):
	assert expected == entropy(seq)


@pytest.mark.parametrize(('expected', 'seq'), bin_entropies)
def test_entropy_counter(expected, seq):
	assert expected == entropy(Counter(seq))


@pytest.mark.parametrize(('expected', 'seq'), quad_entropies)
def test_entropy_quad_base(expected, seq):
	assert expected == entropy(seq, base=4)


@pytest.mark.parametrize(('expected', 'seq'), quad_entropies)
def test_entropy_counter_quad_base(expected, seq):
	assert expected == entropy(Counter(seq), base=4)

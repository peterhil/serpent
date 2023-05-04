from __future__ import annotations

from collections import OrderedDict

import pytest

from serpent.convert.codon import codons_array
from serpent.convert.genetic_code import (
	genetic_code,
	genetic_code_inverse,
	genetic_code_map,
)
from serpent.fun import inverse_od

fixture_codons = codons_array('GACT')
fixture_genetic_table = {
	1: 'GGGGEEDDAAAAVVVVRRSSKKNNTTTTMIIIRRRRQQHHPPPPLLLLW*CC**YYSSSSLLFF',
	2: 'GGGGEEDDAAAAVVVV**SSKKNNTTTTMMIIRRRRQQHHPPPPLLLLWWCC**YYSSSSLLFF',
	5: 'GGGGEEDDAAAAVVVVSSSSKKNNTTTTMMIIRRRRQQHHPPPPLLLLWWCC**YYSSSSLLFF',
}

fixture_code = OrderedDict({
	1: OrderedDict(zip(fixture_codons, fixture_genetic_table[1])),
	2: OrderedDict(zip(fixture_codons, fixture_genetic_table[2])),
	5: OrderedDict(zip(fixture_codons, fixture_genetic_table[5])),
})


def test_genetic_code_map():
	table = fixture_genetic_table[1]
	assert genetic_code_map(table, fixture_codons) == fixture_code[1]


@pytest.mark.parametrize(('n'), [1, 2, 5])
def test_genetic_code(n):
	err_msg = f'Genetic code table {n} differs.'
	assert genetic_code[n] == fixture_code[n], err_msg


@pytest.mark.parametrize(('n'), [1, 2, 5])
def test_genetic_code_inverse(n):
	err_msg = f'Inverse genetic code table {n} differs.'
	assert genetic_code_inverse[n] == inverse_od(fixture_code[n]), err_msg

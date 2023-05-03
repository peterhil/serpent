from __future__ import annotations

from collections import OrderedDict

import pytest

from serpent.convert.codon import codons_array
from serpent.convert.genetic_code import genetic_code

codons = codons_array('GACT')
fixture_genetic_table = {
	1: 'GGGGEEDDAAAAVVVVRRSSKKNNTTTTMIIIRRRRQQHHPPPPLLLLW*CC**YYSSSSLLFF',
	2: 'GGGGEEDDAAAAVVVV**SSKKNNTTTTMMIIRRRRQQHHPPPPLLLLWWCC**YYSSSSLLFF',
	5: 'GGGGEEDDAAAAVVVVSSSSKKNNTTTTMMIIRRRRQQHHPPPPLLLLWWCC**YYSSSSLLFF',
}

fixture_code = OrderedDict({
	1: OrderedDict(zip(codons, fixture_genetic_table[1])),
	2: OrderedDict(zip(codons, fixture_genetic_table[2])),
	5: OrderedDict(zip(codons, fixture_genetic_table[5])),
})


@pytest.mark.parametrize(('n'), [1, 2, 5])
def test_genetic_code(n):
	assert genetic_code[n] == fixture_code[n], f'Genetic code table {n} differs.'

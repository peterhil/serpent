from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from serpent import dna
from serpent.fasta import AMINO, BASE, DEGENERATE
from serpent.fun import str_join


def grey(dna: str) -> str:
	return str_join(np.array(list(dna)).repeat(3))


def test_decode():
	ac = dna.decode('GGGAAACCCTTTGTTATTCTTTTT')
	ex = np.array([ 0, 21, 42, 63, 15, 31, 47, 63])

	assert_array_equal(ac, ex)


def test_decode_degen():
	pass


@pytest.mark.parametrize(('data', 'expected'), [
	(BASE + DEGENERATE, (BASE, DEGENERATE)),
	('ACGTUWYS', ('ACGTU', 'WYS')),
])
def test_clean_non_dna_base(data, expected):
	actual = dna.clean_non_dna(data, amino=False, degen=False)
	assert actual == expected


@pytest.mark.parametrize(('data', 'expected'), [
	(BASE + DEGENERATE + 'J?', (BASE + DEGENERATE, 'J?')),
])
def test_clean_non_dna_base_degen(data, expected):
	actual = dna.clean_non_dna(data, amino=False, degen=True)
	assert actual == expected


@pytest.mark.parametrize(('data', 'expected'), [
	(AMINO + 'J?', (AMINO, 'J?')),
	('WAGMIGTSLSLIIRTELGNPS', ('WAGMIGTSLSLIIRTELGNPS', '')),
])
def test_clean_non_dna_amino(data, expected):
	actual = dna.clean_non_dna(data, amino=True, degen=False)
	assert actual == expected

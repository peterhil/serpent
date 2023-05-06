from __future__ import annotations

import more_itertools as mit
from numpy.testing import assert_array_equal

from serpent import dna
from serpent.convert.codon import codons_array
from serpent.convert.degenerate import degenerate_codons, is_degenerate
from serpent.fun import str_join


def test_dna_decode_is_same_for_degenerate():
	data = 'GGGAAACCCTTT'
	scale = 52

	assert_array_equal(
		dna.decode(data, degen=True) / scale,
		dna.decode(data),
	)


def test_is_degenerate():
	base_sets = [str_join(s) for s in mit.powerset('GACT')]
	for nts in base_sets:
		assert is_degenerate(nts) is False


def test_is_degenerate_for_codons():
	degens = degenerate_codons()
	codons = codons_array()
	for codon in degens:
		assert is_degenerate(codon) == (codon not in codons)

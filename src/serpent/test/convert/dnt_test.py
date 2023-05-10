from __future__ import annotations

import more_itertools as mit

from serpent.convert.codon import codons_array
from serpent.convert.degenerate import degenerate_codons
from serpent.convert.dnt import is_degenerate
from serpent.fun import str_join


def test_is_degenerate():
	base_sets = [str_join(s) for s in mit.powerset('GACT')]
	for nts in base_sets:
		assert is_degenerate(nts) is False


def test_is_degenerate_for_codons():
	degens = degenerate_codons()
	codons = codons_array()
	for codon in degens:
		assert is_degenerate(codon) == (codon not in codons)

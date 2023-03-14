from __future__ import annotations

import numpy as np
from numpy.testing import assert_array_equal

from serpent import dna
from serpent.fun import str_join


def grey(dna: str) -> str:
	return str_join(np.array(list(dna)).repeat(3))


def test_decode():
	ac = dna.decode('GGGAAACCCTTTGTTATTCTTTTT')
	ex = np.array([ 0, 21, 42, 63, 15, 31, 47, 63])

	assert_array_equal(ac, ex)


def test_decode_degen():
	pass

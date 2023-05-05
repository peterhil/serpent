from __future__ import annotations

import numpy as np
from numpy.testing import assert_array_equal

from serpent import dna
from serpent.visual import dna_image_data


def test_dna_image_data():
	# A (bluish) grey band of 4 colours from dark to light:
	codons = dna.decode('GGGGGAGGCAAGAAAAATCCGCCACCCTTATTCTTT')
	assert_array_equal(
		codons,
		np.array([ 0,  1,  2, 20, 21, 23, 40, 41, 42, 61, 62, 63], dtype=np.uint8)
	)

	actual = dna_image_data(codons, width=4)
	expected = np.array([[
		[  1,   5,   9],
		[ 81,  85,  93],
		[161, 165, 169],
		[245, 249, 253]
	]], dtype=np.uint8)

	assert_array_equal(actual, expected)

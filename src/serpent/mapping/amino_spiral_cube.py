from __future__ import annotations

from collections import OrderedDict

import numpy as np

from serpent.fun import str_join

# Spiraling order for each layer of R, N, and Y type of amino acids.
# It has some nice properties when all amino acids are placed on a 3x3x3 cube like this:
#
#       E           A           D
#      / .         / .         / .
#     R   L       T   L       S   F
#    / . / \     / . / \     / . / \
#   K   *   *   P   X   R   N   *   Y
#    \ . . /     \ . . /     \ . . /
#     M   W       S   V       I   C
#      \ /         \ /         \ /
#   /   Q           G           H   \
#  z                                 x
#      (R)         (N)         (Y)
#       y = -1      y = 0       y = 1
#
# The R/N/Y type is determined (loosely) by the degenerate type of the last codon.
# The spiral is also a gray code, except for the center piece at the end.
spiral = np.array([
	# x,  z
	[-1, -1],
	[-1,  0],
	[-1,  1],
	[ 0,  1],
	[ 1,  1],
	[ 1,  0],
	[ 1, -1],
	[ 0, -1],
	[ 0,  0],
], dtype=np.int8)


# TODO Handle stop codons and two groups of L/R/S separately
amino_spiral_order = str_join([
	'ERKMQW*L*',
	'ATPSGVRLX',
	'DSNIHCYF*',
])


coordinates = [
	(y, x, z)
	for y in [-1, 0, 1]
	for (x, z) in spiral
]


amino_spiral = OrderedDict(zip(amino_spiral_order, coordinates))

from __future__ import annotations

from collections import OrderedDict

# Condensed translation table for the Standard Genetic Code
#
# From the NCBI Taxonomy webpage (which also has 25 alternative tables):
# https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1

# AAs    = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
# Starts = ---M------**--*----M---------------M----------------------------
# Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
# Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
# Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG


# ruff: noqa: F601
aminos_inverse = OrderedDict({
	# A
	0:  "K",  # AAA
	1:  "K",  # AAG
	2:  "N",  # AAC
	3:  "N",  # AAU
	4:  "T",  # ACA
	5:  "T",  # ACG
	6:  "T",  # ACC
	7:  "T",  # ACU
	8:  "R",  # AGA
	9:  "R",  # AGG
	10: "S",  # AGC
	11: "S",  # AGU
	12: "I",  # AUA
	13: "M",  # AUG
	14: "I",  # AUC
	15: "I",  # AUU
	# C
	16: "Q",  # CAA
	17: "Q",  # CAG
	18: "H",  # CAC
	19: "H",  # CAU
	20: "P",  # CCA
	21: "P",  # CCG
	22: "P",  # CCC
	23: "P",  # CCU
	24: "R",  # CGA
	25: "R",  # CGG
	26: "R",  # CGC
	27: "R",  # CGU
	28: "L",  # CUA
	29: "L",  # CUG
	30: "L",  # CUC
	31: "L",  # CUU
	# G
	32: "E",  # GAA
	33: "E",  # GAG
	34: "D",  # GAC
	35: "D",  # GAU
	36: "A",  # GCA
	37: "A",  # GCG
	38: "A",  # GCC
	39: "A",  # GCU
	40: "G",  # GGA
	41: "G",  # GGG
	42: "G",  # GGC
	43: "G",  # GGU
	44: "V",  # GUA
	45: "V",  # GUG
	46: "V",  # GUC
	47: "V",  # GUU
	# U
	48: "*",  # UAA
	49: "*",  # UAG
	50: "Y",  # UAC
	51: "Y",  # UAU
	52: "S",  # UCA
	53: "S",  # UCG
	54: "S",  # UCC
	55: "S",  # UCU
	56: "*",  # UGA
	57: "W",  # UGG
	58: "C",  # UGC
	59: "C",  # UGU
	60: "L",  # UUA
	61: "L",  # UUG
	62: "F",  # UUC
	63: "F",  # UUU
	# B = D | N
	3:  "B",  # N
	35: "B",  # D
	# J = L | I
	15: "J",  # I
	61: "J",  # L
	# Z = E | Q
	17: "Z",  # Q
	33: "Z",  # E
	64: "X",  # Any
	65: "-",  # Gap of indeterminmate length
})


aminos = OrderedDict([(v, k) for k, v in aminos_inverse.items()])

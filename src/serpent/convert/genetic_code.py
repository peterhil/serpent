"""Condensed translation tables for the Genetic Codes.

From the NCBI Taxonomy webpage (which also has 25 alternative tables):
https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1

By default all transl_table in GenBank flatfiles are equal to id 1, and this
is not shown. When transl_table is not equal to id 1, it is shown as a
qualifier on the CDS feature.

 1. The Standard Code
 2. The Vertebrate Mitochondrial Code (transl_table=2)
 3. The Yeast Mitochondrial Code
 4. The Mold, Protozoan, and Coelenterate Mitochondrial Code
    and the Mycoplasma/Spiroplasma Code
 5. The Invertebrate Mitochondrial Code
 6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
 9. The Echinoderm and Flatworm Mitochondrial Code
10. The Euplotid Nuclear Code
11. The Bacterial, Archaeal and Plant Plastid Code
12. The Alternative Yeast Nuclear Code
13. The Ascidian Mitochondrial Code
14. The Alternative Flatworm Mitochondrial Code
16. Chlorophycean Mitochondrial Code
21. Trematode Mitochondrial Code
22. Scenedesmus obliquus Mitochondrial Code
23. Thraustochytrium Mitochondrial Code
24. Rhabdopleuridae Mitochondrial Code
25. Candidate Division SR1 and Gracilibacteria Code
26. Pachysolen tannophilus Nuclear Code
27. Karyorelict Nuclear Code
28. Condylostoma Nuclear Code
29. Mesodinium Nuclear Code
30. Peritrich Nuclear Code
31. Blastocrithidia Nuclear Code
33. Cephalodiscidae Mitochondrial UAA-Tyr Code
"""
from __future__ import annotations

from collections import OrderedDict

from serpent.convert.codon import CODONS_LEN, codons_array
from serpent.fun import inverse_od

CODONS_NCBI = codons_array('TCAG')

GENETIC_CODE = {
	# Base  1:  "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
	# Base  2:  "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
	# Base  3:  "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
	1:          "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	# Start 1:  "---M------**--*----M---------------M----------------------------",
	2:          "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
	# Start 2:  "----------**--------------------MMMM----------**---M------------",
	3:          "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	# Start 3:  "----------**----------------------MM---------------M------------",
	4:          "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	# Start 4:  "--MM------**-------M------------MMMM---------------M------------",
	5:          "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
	# Start 5:  "---M------**--------------------MMMM---------------M------------",
	6:          "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	# Start 6:  "--------------*--------------------M----------------------------",
	9:          "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
	# Start 9:  "----------**-----------------------M---------------M------------",
	10:         "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	# Start 10: "----------**-----------------------M----------------------------",
	11:         "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	# Start 11: "---M------**--*----M------------MMMM---------------M------------",
	12:         "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	# Start 12: "----------**--*----M---------------M----------------------------",
	13:         "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
	# Start 13: "---M------**----------------------MM---------------M------------",
	14:         "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
	# Start 14: "-----------*-----------------------M----------------------------",
	16:         "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	# Start 16: "----------*---*--------------------M----------------------------",
	21:         "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
	# Start 21: "----------**-----------------------M---------------M------------",
	22:         "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	# Start 22: "------*---*---*--------------------M----------------------------",
	23:         "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	# Start 23: "--*-------**--*-----------------M--M---------------M------------",
	24:         "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
	# Start 24: "---M------**-------M---------------M---------------M------------",
	25:         "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	# Start 25: "---M------**-----------------------M---------------M------------",
	26:         "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	# Start 26: "----------**--*----M---------------M----------------------------",
	27:         "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	# Start 27: "--------------*--------------------M----------------------------",
	28:         "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	# Start 28: "----------**--*--------------------M----------------------------",
	29:         "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	# Start 29: "--------------*--------------------M----------------------------",
	30:         "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	# Start 30: "--------------*--------------------M----------------------------",
	31:         "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	# Start 31: "----------**-----------------------M----------------------------",
	33:         "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
	# Start 33: "---M-------*-------M---------------M---------------M------------",
	# Base  1:  "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
	# Base  2:  "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
	# Base  3:  "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
}


def create_genetic_table(table: str) -> OrderedDict[str, str]:
	assert len(table) == CODONS_LEN, 'Table should have 64 characters.'
	return OrderedDict(reversed(list(zip(
		CODONS_NCBI,
		table,
	))))


genetic_code = OrderedDict([
	(i, create_genetic_table(table))
	for i, table in GENETIC_CODE.items()
])

genetic_code_inverse = OrderedDict([
	(i, inverse_od(table))
	for i, table in genetic_code.items()
])

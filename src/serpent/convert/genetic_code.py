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
from typing import NamedTuple

import numpy as np

from serpent.convert.codon import CODONS_LEN, codons_array
from serpent.fun import find_offsets, inverse_od, str_join

STANDARD_CODONS = codons_array('GACT')


class GeneticCode(NamedTuple):
	"""Standard genetic code table."""

	code: OrderedDict[str, str]
	start: set
	stop: set


def genetic_code_map(table: str, codons: str=STANDARD_CODONS) -> OrderedDict[str, str]:
	table = OrderedDict(zip(codons, table))
	assert len(table) == CODONS_LEN, 'Codons and table should have 64 characters.'
	return table


def special_set(codons, special: str, marker: str):
	return set(codons[[*find_offsets(special, marker)]])


def code_table(
	code: str, special: str, codons: str=STANDARD_CODONS
) -> OrderedDict[str, str]:
	err_msg = 'Codons and table should have 64 characters.'
	assert len(code) == len(special) == CODONS_LEN, err_msg

	code = genetic_code_map(code, codons)
	start = special_set(codons, special, 'M')
	stop = special_set(codons, special, '*')

	return GeneticCode(code, start, stop)


sgc = code_table


# flake8: noqa ET128
STANDARD_TABLES = OrderedDict({
	#		"GGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCTTTTTTTTTTTTTTTT"  # Base 1
	#		"GGGGAAAACCCCTTTTGGGGAAAACCCCTTTTGGGGAAAACCCCTTTTGGGGAAAACCCCTTTT"  # Base 2
	#		"GACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACT"  # Base 3
	1:  sgc("GGGGEEDDAAAAVVVVRRSSKKNNTTTTMIIIRRRRQQHHPPPPLLLLW*CC**YYSSSSLLFF",
			"----------------------------M---------------M----*--**------M---"),
	2:  sgc("GGGGEEDDAAAAVVVV**SSKKNNTTTTMMIIRRRRQQHHPPPPLLLLWWCC**YYSSSSLLFF",
			"------------M---**----------MMMM--------------------**----------"),
	3:  sgc("GGGGEEDDAAAAVVVVRRSSKKNNTTTTMMIIRRRRQQHHPPPPTTTTWWCC**YYSSSSLLFF",
			"------------M---------------MM----------------------**----------"),
	4:  sgc("GGGGEEDDAAAAVVVVRRSSKKNNTTTTMIIIRRRRQQHHPPPPLLLLWWCC**YYSSSSLLFF",
			"------------M---------------MMMM------------M-------**------MM--"),
	5:  sgc("GGGGEEDDAAAAVVVVSSSSKKNNTTTTMMIIRRRRQQHHPPPPLLLLWWCC**YYSSSSLLFF",
			"------------M---------------MMMM--------------------**------M---"),
	6:  sgc("GGGGEEDDAAAAVVVVRRSSKKNNTTTTMIIIRRRRQQHHPPPPLLLLW*CCQQYYSSSSLLFF",
			"----------------------------M--------------------*--------------"),
	9:  sgc("GGGGEEDDAAAAVVVVSSSSKNNNTTTTMIIIRRRRQQHHPPPPLLLLWWCC**YYSSSSLLFF",
			"------------M---------------M-----------------------**----------"),
	10: sgc("GGGGEEDDAAAAVVVVRRSSKKNNTTTTMIIIRRRRQQHHPPPPLLLLWCCC**YYSSSSLLFF",
			"----------------------------M-----------------------**----------"),
	11: sgc("GGGGEEDDAAAAVVVVRRSSKKNNTTTTMIIIRRRRQQHHPPPPLLLLW*CC**YYSSSSLLFF",
			"------------M---------------MMMM------------M----*--**------M---"),
	12: sgc("GGGGEEDDAAAAVVVVRRSSKKNNTTTTMIIIRRRRQQHHPPPPSLLLW*CC**YYSSSSLLFF",
			"----------------------------M---------------M----*--**----------"),
	13: sgc("GGGGEEDDAAAAVVVVGGSSKKNNTTTTMMIIRRRRQQHHPPPPLLLLWWCC**YYSSSSLLFF",
			"------------M---------------MM----------------------**------M---"),
	14: sgc("GGGGEEDDAAAAVVVVSSSSKNNNTTTTMIIIRRRRQQHHPPPPLLLLWWCC*YYYSSSSLLFF",
			"----------------------------M-----------------------*-----------"),
	16: sgc("GGGGEEDDAAAAVVVVRRSSKKNNTTTTMIIIRRRRQQHHPPPPLLLLW*CCL*YYSSSSLLFF",
			"----------------------------M--------------------*---*----------"),
	21: sgc("GGGGEEDDAAAAVVVVSSSSKNNNTTTTMMIIRRRRQQHHPPPPLLLLWWCC**YYSSSSLLFF",
			"------------M---------------M-----------------------**----------"),
	22: sgc("GGGGEEDDAAAAVVVVRRSSKKNNTTTTMIIIRRRRQQHHPPPPLLLLW*CCL*YYS*SSLLFF",
			"----------------------------M--------------------*---*---*------"),
	23: sgc("GGGGEEDDAAAAVVVVRRSSKKNNTTTTMIIIRRRRQQHHPPPPLLLLW*CC**YYSSSSL*FF",
			"------------M---------------M--M-----------------*--**-------*--"),
	24: sgc("GGGGEEDDAAAAVVVVKSSSKKNNTTTTMIIIRRRRQQHHPPPPLLLLWWCC**YYSSSSLLFF",
			"------------M---------------M---------------M-------**------M---"),
	25: sgc("GGGGEEDDAAAAVVVVRRSSKKNNTTTTMIIIRRRRQQHHPPPPLLLLWGCC**YYSSSSLLFF",
			"------------M---------------M-----------------------**------M---"),
	26: sgc("GGGGEEDDAAAAVVVVRRSSKKNNTTTTMIIIRRRRQQHHPPPPALLLW*CC**YYSSSSLLFF",
			"----------------------------M---------------M----*--**----------"),
	27: sgc("GGGGEEDDAAAAVVVVRRSSKKNNTTTTMIIIRRRRQQHHPPPPLLLLWWCCQQYYSSSSLLFF",
			"----------------------------M--------------------*--------------"),
	28: sgc("GGGGEEDDAAAAVVVVRRSSKKNNTTTTMIIIRRRRQQHHPPPPLLLLWWCCQQYYSSSSLLFF",
			"----------------------------M--------------------*--**----------"),
	29: sgc("GGGGEEDDAAAAVVVVRRSSKKNNTTTTMIIIRRRRQQHHPPPPLLLLW*CCYYYYSSSSLLFF",
			"----------------------------M--------------------*--------------"),
	30: sgc("GGGGEEDDAAAAVVVVRRSSKKNNTTTTMIIIRRRRQQHHPPPPLLLLW*CCEEYYSSSSLLFF",
			"----------------------------M--------------------*--------------"),
	31: sgc("GGGGEEDDAAAAVVVVRRSSKKNNTTTTMIIIRRRRQQHHPPPPLLLLWWCCEEYYSSSSLLFF",
			"----------------------------M-----------------------**----------"),
	33: sgc("GGGGEEDDAAAAVVVVKSSSKKNNTTTTMIIIRRRRQQHHPPPPLLLLWWCC*YYYSSSSLLFF",
			"------------M---------------M---------------M-------*-------M---"),
	#		"GGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCTTTTTTTTTTTTTTTT"  # Base 1
	#		"GGGGAAAACCCCTTTTGGGGAAAACCCCTTTTGGGGAAAACCCCTTTTGGGGAAAACCCCTTTT"  # Base 2
	#		"GACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACT"  # Base 3
})


GENETIC_CODE = OrderedDict([
	(key, table.code)
	for key, table in STANDARD_TABLES.items()
])

genetic_code = GENETIC_CODE
genetic_code_inverse = OrderedDict([
	(key, inverse_od(table))
	for key, table in genetic_code.items()
])

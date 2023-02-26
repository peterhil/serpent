#!/usr/bin/env python

import re

from typing import NamedTuple


DATA_TOKENS = ["AMINO", "BASE", "DEGENERATE"]


class Token(NamedTuple):
	type: str
	value: str
	line: int
	column: int


def tokenize(data, amino=False):
	"""Iterative FASTA sequence reader.
	See: https://docs.python.org/3/library/re.html#writing-a-tokenizer
	"""
	BASE = r"[ACGTU\n]"
	token_specification = [
		("DESCRIPTION", r">[^\n]*"),
		("BASE", BASE + r"+"),
		(
			"DEGENERATE",
			r"[WSMKRYBDHVNZ-]+?",
		),  # https://en.wikipedia.org/wiki/Nucleic_acid_sequence#Notation
		("NEWLINE", r"\n"),  # Line endings
		("SKIP", r"[ \t]+"),  # Skip over spaces and tabs
		("MISMATCH", r"."),  # Any other character
	]
	if amino:
		# https://en.wikipedia.org/wiki/Proteinogenic_amino_acid
		token_specification.insert(1, ("AMINO", r"[ARNDCQEGHILKMFPSTWYVUO]+"))

	tok_regex = "|".join("(?P<%s>%s)" % pair for pair in token_specification)
	line_num = 1
	line_start = 0
	for mo in re.finditer(tok_regex, data):
		kind = mo.lastgroup
		value = mo.group()
		column = mo.start() - line_start
		if kind == "NEWLINE":
			line_start = mo.end()
			line_num += 1
			continue
		elif kind == "SKIP":
			continue
		elif kind == "MISMATCH":
			raise RuntimeError(f"{value!r} unexpected on line {line_num}")
		yield Token(kind, value, line_num, column)

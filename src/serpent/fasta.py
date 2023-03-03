from __future__ import annotations

import re
from pathlib import Path
from typing import NamedTuple

DATA_TOKENS = ["AMINO", "BASE", "DEGENERATE"]


class Token(NamedTuple):
	type: str
	value: str
	line: int
	column: int

	@property
	def is_data(self):
		return self.type in DATA_TOKENS


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
			err_msg = f"{value!r} unexpected on line {line_num} column {column}"
			raise ValueError(err_msg)
		yield Token(kind, value, line_num, column)


def read(fn: str, amino: bool = False) -> str:
	data = []
	description_count = 0

	with Path(fn).open(encoding="UTF-8") as file:
		while (line := file.readline().rstrip()) and description_count < 2:
			for token in tokenize(line, amino):
				if token.type == "DESCRIPTION":
					description_count += 1
				if description_count == 2:
					break
				if token.is_data:
					data.append(token.value)

	# TODO Create a TUI or add CLI option to select sequences or
	# otherwise handle multiple sequences
	if description_count >= 2:
		print("Warning: File has more than one FASTA sequence!")

	return "\n".join(data)

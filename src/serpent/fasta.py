"""FASTA data io."""

from __future__ import annotations

import re
from collections.abc import Iterator
from pathlib import Path
from typing import NamedTuple

DATA_TOKENS = ["AMINO", "BASE", "DEGENERATE"]


class Token(NamedTuple):
	"""Token of parsed FASTA data."""

	type: str
	value: str
	line: int
	column: int

	@property
	def is_data(self):
		"""True if the token contains sequence data."""
		return self.type in DATA_TOKENS


def tokenize(data: str, amino: bool=False) -> Iterator[Token]:
	"""Read FASTA sequences iteratively.

	See: https://docs.python.org/3/library/re.html#writing-a-tokenizer.
	"""
	base = r"[ACGTU\n]"
	token_specification = [
		("DESCRIPTION", r">[^\n]*"),
		("BASE", base + r"+"),
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
		token_specification.insert(1, ("AMINO", r"[ARNDCQEGHILKMFPSTWYVUOX*-]+"))

	tok_regex = "|".join("(?P<%s>%s)" % pair for pair in token_specification)
	line_num: int = 1
	line_start: int = 0
	for matches in re.finditer(tok_regex, data):
		kind: str = matches.lastgroup or "MISMATCH"
		value: str = matches.group()
		column: int = matches.start() - line_start
		if kind == "NEWLINE":
			line_start = matches.end()
			line_num += 1
			continue
		elif kind == "SKIP":
			continue
		elif kind == "MISMATCH":
			err_msg = f"{value!r} unexpected on line {line_num} column {column}"
			raise ValueError(err_msg)
		yield Token(kind, value, line_num, column)


def read(filename: str, amino: bool = False) -> str:
	"""Read FASTA data.

	Amino:
		false = input contains nucleotide bases
		true  = input contains amino acids
	"""
	data = []
	description_count = 0
	max_count = 2000

	with Path(filename).open(encoding="UTF-8") as file:
		while (line := file.readline().rstrip()) and description_count < max_count:
			for token in tokenize(line, amino):
				if token.type == "DESCRIPTION":
					description_count += 1
				if description_count == max_count:
					break
				if token.is_data:
					data.append(token.value)

	# TODO Create a TUI or add CLI option to select sequences or
	# otherwise handle multiple sequences
	if description_count >= max_count:
		print("Warning: File has more than {max_count} FASTA sequences!")

	return "\n".join(data)

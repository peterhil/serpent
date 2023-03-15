"""FASTA data io."""

from __future__ import annotations

import re
from collections.abc import Iterator
from fileinput import FileInput
from pathlib import Path
from typing import NamedTuple

from termcolor import colored

from serpent.io import err
from serpent.settings import DEFAULT_TERM_COLOR

DATA_TOKENS = [
	"AMINO",
	"BASE",
	"DEGENERATE",
]

# IUPAC encodings
#
# Amino acids:
# https://en.wikipedia.org/wiki/Proteinogenic_amino_acid
#
# DNA/RNA with degenarate data:
# https://en.wikipedia.org/wiki/Nucleic_acid_sequence#Notation
#
# Note: Keep `-` at the end for regexp character ranges!
AMINO = "ARNDCQEGHILKMFPSTWYVUOX*-"  # TODO Handle 'X*-' as degenerate?
BASE = "ACGTU"
DEGENERATE = "WSMKRYBDHVNZ-"

RE_DESCRIPTION = r"^[>;](?P<description>.*)\n?"
RE_DEGENERATE = fr"[{DEGENERATE}]+?"


def get_description(string: str) -> str | None:
	re_description = re.compile(RE_DESCRIPTION)
	matches = re_description.match(string)

	return matches.group(1).strip() if matches else None


def find_fasta_files(fi: FileInput, debug=False) -> Iterator[str]:
	for line in fi:
		filename = fi.filename()
		if fi.isfirstline() and get_description(line):
			yield filename
		elif debug:
			err(f"{filename}: Not FASTA!")
		fi.nextfile()


def find_fasta_sequences(fi: FileInput, debug=False) -> Iterator[str]:
	"""Find FASTA files and print single and multiple sequences."""
	for line in fi:
		description = get_description(line)
		newfile = fi.isfirstline()
		if description:
			if newfile:
				filename = fi.filename()
				yield colored(filename, DEFAULT_TERM_COLOR)
			yield f">{description}"
		elif newfile:
			fi.nextfile()
		if not description and debug:  # TODO Use with cat
			yield line.rstrip()


class ParseError(Exception):
	pass


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


def tokenize(data: str, amino: bool=False, line: int=1) -> Iterator[Token]:
	"""Read FASTA sequences iteratively.

	See:
	https://en.wikipedia.org/wiki/FASTA_format
	https://docs.python.org/3/library/re.html#writing-a-tokenizer.
	"""
	token_specification = [
		("DESCRIPTION", RE_DESCRIPTION),
		("BASE", fr"[{BASE}]+"),
		("DEGENERATE", RE_DEGENERATE),
		("NEWLINE", r"\n"),  # Line endings
		("SKIP", r"[ \t]+"),  # Skip over spaces and tabs
		("MISMATCH", r"."),  # Any other character
	]
	if amino:
		token_specification.insert(1, ("AMINO", fr"[{AMINO}]+"))

	spec = "|".join("(?P<%s>%s)" % pair for pair in token_specification)
	# TODO Handle lowercase insertions better, see FASTA format
	rec_token = re.compile(spec, flags=re.I)
	# line: int = 1
	line_start: int = 0
	for matches in rec_token.finditer(data):
		kind: str = matches.lastgroup or "MISMATCH"
		value: str = matches.group()
		column: int = matches.start() - line_start
		if kind == "NEWLINE":
			line_start = matches.end()
			line += 1
			continue
		elif kind == "SKIP":
			continue
		elif kind == "MISMATCH":
			print(data)
			err_msg = f"{value!r} unexpected on line {line} column {column}"
			raise ParseError(err_msg)
		yield Token(kind, value, line, column)


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
		print(f"Warning: File has more than {max_count} FASTA sequences!")

	return "\n".join(data)

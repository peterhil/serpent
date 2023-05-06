"""FASTA data io."""

from __future__ import annotations

import re
from collections import OrderedDict
from collections.abc import Iterable, Iterator, Sequence
from fileinput import FileInput
from os import PathLike
from pathlib import Path
from typing import NamedTuple

import more_itertools as mit
from termcolor import colored

from serpent.fun import str_join
from serpent.io.files import err, info
from serpent.settings import DEBUG, DEFAULT_TERM_COLOR

DATA_TOKENS = [
	"AMINO",
	"AMINO_DEGENERATE",
	"BASE",
	"BASE_NONCODING",
	"DEGENERATE",
]

# IUPAC encodings
#
# Amino acids:
# https://en.wikipedia.org/wiki/Proteinogenic_amino_acid
#
# Degenerate amino acids notation (see third table: “Ambiguous amino acids”)
# https://en.wikipedia.org/wiki/Amino_acid#Table_of_standard_amino_acid_abbreviations_and_properties
#
# DNA/RNA with degenarate data:
# https://en.wikipedia.org/wiki/Nucleic_acid_sequence#Notation
#
# Note: Keep `-` at the end for regexp character ranges!
AMINO = "GEDAVRSKNTMIQHPLW*CYFUO"
AMINO_DEGENERATE = "BJZX-"
BASE = "ACGTU"
BASE_NONCODING = BASE.lower()
DEGENERATE = "WSMKRYBDHVNZ-"

RE_DESCRIPTION = r"^[>;@].*"
RE_WHITESPACE = re.compile(r'\s')


def auto_select_amino(filename: str, amino: bool) -> bool:
	if (match := re.search(r'\.f([an])a$', str(filename).lower())):
		amino = match.group(1) == 'a'
		info(f'info: {filename}: Auto selected amino = {amino}')
	return amino


def get_data(seq: Sequence[FastaToken]) -> str:
	data = str_join(token.data for token in seq if token.data)

	return data


def get_description(string: str) -> str | None:
	re_description = re.compile(fr"(?P<DESCRIPTION>{RE_DESCRIPTION})")
	matches = re_description.match(string)

	return matches.group(1).strip() if matches else None


def find_fasta_files(fi: FileInput, debug=DEBUG) -> Iterator[str]:
	for line in fi:
		filename = fi.filename()
		if fi.isfirstline() and get_description(line):
			yield filename
		elif debug:
			err(f"err: {filename}: Not FASTA!")
		fi.nextfile()


def find_fasta_sequences(fi: FileInput, num=False, debug=DEBUG) -> Iterator[str]:
	"""Find FASTA files and print single and multiple sequences."""
	for line in fi:
		description = get_description(line)
		newfile = fi.isfirstline()
		if description:
			if newfile:
				filename = fi.filename()
				yield colored(f'\t{filename}' if num else filename, DEFAULT_TERM_COLOR)
			if num:
				yield f"{fi.filelineno()}\t{description}"
			else:
				yield description
		elif newfile:
			fi.nextfile()
		if not description and debug:  # TODO Use with cat
			yield line.rstrip()


class ParseError(Exception):
	pass


class FastaToken(NamedTuple):
	"""Token of parsed FASTA data."""

	type: str
	line: int
	column: int
	data: str = ''
	value: str = ''

	@property
	def is_data(self):
		"""True if the token contains sequence data."""
		return self.type in DATA_TOKENS

	@property
	def is_description(self):
		"""True if the token is a description."""
		return self.type == 'DESCRIPTION'


def re_token(group):
	return fr"[{group}]+"


def get_token_specification(amino: bool=False):
	"""Get a regexp for parsing FASTA files.

	See:
	https://en.wikipedia.org/wiki/FASTA_format
	"""
	aminos = [
		("AMINO", re_token(AMINO)),
		("AMINO_DEGENERATE", re_token(AMINO_DEGENERATE)),
	]
	nucleotides = [
		("BASE", re_token(BASE)),
		("BASE_NONCODING", re_token(BASE_NONCODING)),
		("DEGENERATE", re_token(DEGENERATE)),
	]

	token_specification = OrderedDict(aminos if amino else nucleotides)
	token_specification.update([
		("DESCRIPTION", RE_DESCRIPTION),
		("NEWLINE", r"\n"),  # Line endings
		("SKIP", r"[ \t]+"),  # Skip over spaces and tabs
		("MISMATCH", r"."),  # Any other character
	])
	spec = "|".join((
		fr"(?P<{k}>{v})"
		for k, v in token_specification.items()
	))

	return re.compile(spec)


RE_AMINO = get_token_specification(amino=True)
RE_NUCLEOTIDE = get_token_specification(amino=False)


def tokenize(
	data: str, amino: bool=False, line: int=1
) -> Iterator[FastaToken]:
	"""Read FASTA sequences iteratively.

	See:
	https://docs.python.org/3/library/re.html#writing-a-tokenizer
	"""
	# TODO Handle lowercase insertions better, see FASTA format
	spec = RE_AMINO if amino else RE_NUCLEOTIDE
	# line: int = 1
	line_start: int = 0
	for matches in spec.finditer(data):
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
		elif kind in DATA_TOKENS:
			token = FastaToken(kind, line, column, data=value.upper())
			if DEBUG:
				print(token)
			yield token
		else:
			yield FastaToken(kind, line, column, value=value)


def read(filename: PathLike, amino: bool = False) -> str:
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
					data.append(token.data)

	# TODO Create a TUI or add CLI option to select sequences or
	# otherwise handle multiple sequences
	if description_count >= max_count:
		print(f"Warning: File has more than {max_count} FASTA sequences!")

	return "\n".join(data)


def read_tokens(filename: PathLike, amino: bool = False) -> Iterable[FastaToken]:
	"""Parse data as tokens from a FASTA file.

	Amino:
		false = input contains nucleotide bases
		true  = input contains amino acids
	"""
	lineno = 0
	with Path(filename).open(encoding="UTF-8") as file:
		while (line := file.readline()):
			lineno += 1
			for token in tokenize(line, amino, lineno):
				if token.type == "DESCRIPTION" or token.is_data:
					yield token
				else:
					print('Extra:', token)


def read_sequences(filename: PathLike, amino: bool=False) -> Iterable[list[FastaToken]]:
	"""Read sequences from a FASTA file."""
	tokens = read_tokens(filename, amino)
	sequences = mit.split_before(tokens, lambda token: token.is_description)

	yield from sequences


def descriptions_and_data(sequence):
	"""Partition a sequence into data and descriptions."""
	[tokens, descriptions] = mit.partition(lambda t: t.is_description, sequence)

	descriptions = (desc.value for desc in descriptions)
	# FIXME Read data iteratively by removing str_join (which breaks things)
	data = str_join(token.data for token in tokens if token.data)

	yield from (descriptions, data)

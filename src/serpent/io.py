from __future__ import annotations

from collections.abc import Iterator
from fileinput import FileInput, hook_encoded
from pathlib import Path
from sys import argv, stderr, stdout

from termcolor import colored

from serpent.config import DEFAULT_TERM_COLOR
from serpent.fasta import get_description

# TODO Check out hook_compressed and zzzip
openhook = hook_encoded("utf-8", "surrogateescape")


def echo(message):
	return stdout.write(f"{message}\n")


def err(message):
	return stderr.write(f"{message}\n")


def check_path(path: Path, recurse: bool=False) -> Iterator[Path]:
	if not path.exists():
		err(f"{argv[0]}: {path}: No such file or directory")
	if path.is_dir():
		if recurse:
			for file in sorted(path.iterdir()):
				yield from check_path(Path(file), recurse)
		else:
			err(f"{argv[0]}: {path}: Is directory")
	if path.is_file():
		yield path


def check_inputs(inputs: list[str], recurse: bool=False):
	for input in inputs:
		yield from check_path(Path(input), recurse)


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

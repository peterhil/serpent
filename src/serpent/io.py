from __future__ import annotations

from collections.abc import Iterator, Sequence
from fileinput import hook_encoded
from pathlib import Path
from sys import argv, stderr, stdout

# TODO Check out hook_compressed and zzzip
openhook = hook_encoded("utf-8", "surrogateescape")


def echo(message):
	return stdout.write(f"{message}\n")


def err(message):
	return stderr.write(f"{message}\n")


def check_path(input: Path | str, recurse: bool=False) -> Iterator[Path | str]:
	path = Path(input)

	if input == '-':
		yield '-'
		return
	elif not path.exists():
		return err(f"{argv[0]}: {path}: No such file or directory")

	if path.is_dir():
		if recurse:
			for file in sorted(path.iterdir()):
				yield from check_path(file, recurse)
		else:
			return err(f"{argv[0]}: {path}: Is directory")
	elif path.is_file():
		yield path
	else:
		return err(f'{path}: Is not a regular file')


def check_paths(inputs: Sequence[str], recurse: bool=False) -> Iterator[Path | str]:
	if len(inputs) == 0:
		yield '-'
	for input in inputs:
		yield from check_path(input, recurse)

from __future__ import annotations

from collections.abc import Iterator, Sequence
from fileinput import hook_encoded
from functools import partial
from pathlib import Path
from sys import argv, stdin

from serpent.io.printing import err
from serpent.settings import BASE_ORDER

# TODO Check out hook_compressed and zzzip
openhook = hook_encoded("utf-8", "surrogateescape")


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


def file_extension_for(fmt: str = 'codon'):
	"""Get file extension for encoded data format."""
	if fmt in ['a', 'amino']:
		extension = 'faa'
	elif fmt in ['c', 'codon']:
		extension = 'fna'
	else:
		raise ValueError('Unknown format: ' + fmt)

	return extension


# ruff: noqa: PLR0913 # Too many arguments in function definition
def image_name_for(
	filename, width=0, mode='RGB',
	*,
	amino=False, degen=False, table=1, length=1
):
	"""Add extra info for image variation to filename."""
	if mode == 'P':
		palette = ('Pa' if amino else 'Pn')
	elif mode == 'Q':
		palette = f'Q{length}'
	elif mode == 'L':
		palette = 'L'
	else:
		palette = ''
	code_table = f't{table}' if amino and table != 1 else ''

	outfile = '.'.join(filter(None, [
		filename,
		palette,
		f'w{width}' if width else '',
		f'{BASE_ORDER}g' if degen else BASE_ORDER,
		code_table,
		'png'
	]))

	return outfile


def write_iterable(lines, outfile):
	with Path(outfile).open("w", encoding="UTF-8") as f:
		writeln = partial(print, file=f)
		for line in lines:
			writeln(line)
	print(f"Wrote: {outfile}")


def openfile(filename):
	"""Open stdin or a file to be used in a context handler."""
	# ruff: noqa: SIM115
	return stdin if filename == '-' else Path(filename).open("r", encoding="UTF-8")


def readlines(filename):
	with openfile(filename) as f:
		while (line := f.readline()):
			yield line.rstrip()

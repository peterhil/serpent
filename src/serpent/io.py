from __future__ import annotations

from collections.abc import Iterator, Sequence
from fileinput import hook_encoded
from pathlib import Path
from sys import argv, stderr, stdout

from serpent.settings import BASE_ORDER

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


def file_extension_for(fmt: str = 'base64'):
	"""Get file extension for encoded data format."""
	if fmt in ['b', 'base64']:
		extension = 'ser64'
	elif fmt in ['c', 'codon']:
		extension = 'codon.fasta'
	else:
		raise ValueError('Unknown format: ' + fmt)

	return extension


def image_name_for(filename, width=0, mode='RGB', *, amino=False, table=1, length=1):
	"""Add extra info for image variation to filename."""
	if mode == 'P':
		palette = ('Pa' if amino else 'Pn')
	elif mode == 'Q':
		palette = f'Q{length}'
	else:
		palette = ''
	code_table = f't{table}' if amino and table != 1 else ''

	outfile = '.'.join(filter(None, [
		filename,
		palette,
		f'w{width}' if width else '',
		BASE_ORDER,
		code_table,
		'png'
	]))

	return outfile


def wait_user():
	"""Wait for user to enter any key.

	Useful for allowing the user time to explore the interactive plots.
	"""
	return input('Press ENTER when ready. ')

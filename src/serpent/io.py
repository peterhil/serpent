from __future__ import annotations

from pathlib import Path
from sys import argv, stderr


def err(message):
	return stderr.write(f"{message}\n")


def check_inputs(inputs):
	for input in inputs:
		path = Path(input)
		if not path.exists():
			err(f"{argv[0]}: {path}: No such file or directory")
		if path.is_dir():
			err(f"{argv[0]}: {path}: Is directory")
		if path.is_file():
			yield path

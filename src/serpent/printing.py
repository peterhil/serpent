from __future__ import annotations

from collections import Counter
from collections.abc import Iterable

import more_itertools as mit
import numpy as np

from serpent.fun import str_join
from serpent.settings import COUNT_LIMIT


def format_data(data, width=80, sep=' '):
	"""Format data for printing."""
	def fmt(chunk):
		return str_join(chunk, sep)
	lines = [
		fmt(chunk)
		for i, chunk
		in enumerate(mit.chunked(data, width))
	]

	return lines


def format_lines(data, width=80, sep=' '):
	"""Format lines for printing."""
	# TODO Just add line numbers here to lines
	def fmt(chunk):
		return str_join(chunk, sep)
	lines = [
		f"{i * width}:\t{fmt(chunk)}"
		for i, chunk
		in enumerate(mit.chunked(data, width))
	]

	return lines


def format_counts(counts: Counter, limit=COUNT_LIMIT) -> Iterable[str]:
	yield from (
		f'{item}\t{count}'
		for (item, count)
		in mit.take(limit, counts.most_common())
	)

def format_decoded(decoded, degen=False):
	uniq = len(np.unique(decoded))
	print(f"Decoded ({len(decoded)}, unique: {uniq}):")
	strings = (f'{d: >4}' if degen else f'{d: >2}' for d in iter(decoded))
	lines = format_lines(strings, 32)

	return lines


def format_quasar(data):
	strings = (f'{d: >3}' for d in iter(data))
	lines = format_data(strings, 32)

	return lines


def format_quasar_pulses(pulses, height):
	for row in range(height):
		yield from format_quasar([pulses[col][row] for col in pulses])


def format_split(regions, width=72, split='n'):
	for i, region in enumerate(regions):
		yield f'@split-{split}-{i}'
		lines = (str_join(line) for line in mit.chunked(str_join(region), width))
		yield from lines

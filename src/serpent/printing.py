from __future__ import annotations

import more_itertools as mit
import numpy as np

from serpent.fun import str_join
from serpent.mathematics import percent
from serpent.stats import gc_content


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


def reflow(data, width=80):
	"""Reflow text data into line width."""
	# TODO Replace str_join with an iterative solution
	return (str_join(line) for line in mit.chunked(str_join(data), width))


def format_counter(counts, show_gc=False, *, limit=None):
	yield from (f'{count}\t{symbol}' for symbol, count in counts.most_common(limit))
	if show_gc:
		gc_percent = percent(gc_content(counts), 3)
		yield f'GC content: {gc_percent}%'


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

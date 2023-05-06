from __future__ import annotations

import blessed
import more_itertools as mit
import numpy as np

from serpent.fun import str_join
from serpent.math.basic import percent
from serpent.math.statistic import gc_content


def auto_line_width(item_size, base=8, indent=0, sep=1):
	term_width = blessed.Terminal().width
	unrounded_width = (term_width + sep - indent) // (item_size + sep)
	width = base * max(1, (unrounded_width // base))

	return width


def auto_line_width_for(fmt, seql, base=8, indent=8, sep=1):
	item_size = seql * 3 if fmt in ['c', 'codon'] else seql
	width = auto_line_width(item_size, base, indent, sep)
	return width


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


def format_counter(counts, show_gc=False, *, limit=1, top=None, precision=2):
	total = np.sum([*counts.values()])

	yield 'symbol\tcount\tpercent'
	yield from (
		f'{symbol}\t{count}\t{percent(count / total, precision)}%'
		for symbol, count in counts.most_common(top)
		if count >= limit
	)
	if show_gc:
		gc = gc_content(counts)
		gc_percent = percent(gc / total, precision)
		yield f'GC\t{gc}\t{gc_percent}%'

	yield f'total\t{total}'


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

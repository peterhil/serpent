from __future__ import annotations

import colorsys
import math
from dataclasses import dataclass

import blessed
import numpy as np

from serpent import dna
from serpent.fun import str_join
from serpent.io.fasta import auto_select_amino, descriptions_and_data, read_sequences
from serpent.io.files import check_paths
from serpent.math.basic import rescale
from serpent.visual import ansi
from serpent.visual.bitmap import decoded_to_pixels
from serpent.visual.block_elements import HALF_BLOCK, pixels_to_blocks


@dataclass
class ZigzagState:
	"""Dataclass for keeping track of the zigzag state."""

	inputs: list[str]
	dirty: bool = True
	file_no: int = 0
	page_no: int = 0

	@property
	def current_input(self):
		return self.inputs[self.file_no]

	@property
	def total(self):
		return len(self.inputs) or 1

	def next_input(self):
		self.file_no = (self.file_no + 1) % self.total
		self.page_no = 0
		self.dirty = True

	def prev_input(self):
		self.file_no = (self.file_no - 1) % self.total
		self.page_no = 0
		self.dirty = True


def rgb_at_xy(term, x, y, t):
	h, w = term.height, term.width
	xw = x - w / 2
	yh = y - h / 2

	hue = 4.0 + (
		math.sin(x / 16)
		+ math.sin(y / 32)
		+ math.sin(math.sqrt(xw ** 2 + yh ** 2) / 8 + t * 3)
	) + math.sin(math.sqrt(x ** 2 + y ** 2) / 8)
	saturation = y / h
	lightness = x / w

	rgb = np.array(colorsys.hsv_to_rgb(hue / 8, saturation, lightness))
	rgb = rescale(rgb, 256, 65536).astype(np.int32)

	return rgb


def screen_page(term, render_fn, t):
	result = ''
	for y in range(term.height - 1):
		for x in range(term.width):
			# fg = term.color_rgb(*render_fn(term, x, y, t))
			# bg = term.on_color_rgb(*render_fn(term, x, y, t + 1))
			# result += fg + bg + HALF_BLOCK
			fg = render_fn(term, x, y, t)
			bg = render_fn(term, x, y, t + 1)
			colours = ansi.rgb(fg, bg)
			result += colours + HALF_BLOCK
	return result


# ruff: noqa: PLR0913 # Too many arguments in function definition
def page(
	term,
	state,
	*,
	width=64, mode='RGB',
	amino=False, degen=False, table=1,
):
	height = term.height - 1
	filename = state.current_input

	amino = auto_select_amino(filename, amino)
	seqs = read_sequences(filename, amino)

	for sequence in seqs:
		[descriptions, data] = descriptions_and_data(sequence)
		decoded = dna.decode(data, amino, table, degen)
		pixels = decoded_to_pixels(decoded, mode, amino, degen)
		yield from pixels_to_blocks(pixels, width, height=height, mode=mode)


def status(term, state):
	left_txt = f'file ({state.file_no + 1} / {state.total}): {state.current_input}'
	right_txt = f'width {term.width}; {term.number_of_colors} colors - ?: help'
	return (
		'\n' + term.normal +
		term.white_on_purple + term.clear_eol +
		left_txt +
		term.rjust(right_txt, term.width - len(left_txt))
	)


# ruff: noqa: PLR0913 # Too many arguments in function definition
def zigzag_blocks(
	inputs,
	*,
	width=64, mode='RGB',
	amino=False, degen=False, table=1,
):
	"""Browse DNA data as text paged into variable line widths."""
	term = blessed.Terminal()
	state = ZigzagState(
		inputs = [*map(str, check_paths(inputs))]
	)

	with term.cbreak(), term.hidden_cursor(), term.fullscreen():
		state.dirty = True
		while True:
			if state.dirty:
				outp = term.home
				outp += str_join(page(
					term,
					state,
					width=width,
					mode=mode,
					amino=amino,
					degen=degen,
					table=table,
				))
				outp += status(term, state)
				# print(outp, end='')
				# sys.stdout.flush()
				yield outp
				state.dirty = False

			key = term.inkey(timeout=None)
			if key == 'n':
				state.next_input()
			elif key == 'p':
				state.prev_input()
			elif key == 'q':
				break

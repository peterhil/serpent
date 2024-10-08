from __future__ import annotations

import colorsys
import fileinput
import math
import sys
from dataclasses import dataclass

import blessed
import numpy as np

from serpent.io.files import check_paths
from serpent.math.basic import rescale
from serpent.visual.block_elements import HALF_BLOCK


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
			fg = term.color_rgb(*render_fn(term, x, y, t))
			bg = term.on_color_rgb(*render_fn(term, x, y, t + 1))
			result += bg + fg + HALF_BLOCK
	return result


def status(term, state):
	left_txt = f'file ({state.file_no + 1} / {state.total}): {state.current_input}'
	right_txt = f'{term.number_of_colors} colors - ?: help'
	return (
		'\n' + term.normal +
		term.white_on_purple + term.clear_eol +
		left_txt +
		term.rjust(right_txt, term.width - len(left_txt))
	)


def zigzag_blocks(inputs):
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
				outp += screen_page(term, rgb_at_xy, state.file_no)
				outp += status(term, state)
				print(outp, end='')
				sys.stdout.flush()
				state.dirty = False

			key = term.inkey(timeout=None)
			if key == ' ':
				state.next_input()
			elif key == 'b':
				state.prev_input()
			elif key == 'q':
				break


def zigzag_text(inputs):
	"""Browse DNA data as text paged into variable line widths."""
	inputs = check_paths(inputs)

	# Read files in binary mode
	with fileinput.input(inputs, mode='rb') as fi:
		while (line := fi.readline()):
			if fi.isstdin():
				print('STDIN: ' , end='')
			print(str(line, encoding='UTF-8', errors='surrogateescape'), end='')

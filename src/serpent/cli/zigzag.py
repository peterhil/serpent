from __future__ import annotations

import sys
from dataclasses import dataclass

import blessed

from serpent import dna
from serpent.fun import str_join
from serpent.io.fasta import auto_select_amino, descriptions_and_data, read_sequences
from serpent.io.files import check_paths
from serpent.visual.bitmap import decoded_to_pixels
from serpent.visual.block_elements import pixels_to_blocks


@dataclass
class ZigzagState:
	"""Dataclass for keeping track of the zigzag state."""

	inputs: list[str]
	width: int = 80
	dirty: bool = True
	file_no: int = 0
	page_no: int = 0

	@property
	def current_input(self):
		return self.inputs[self.file_no]

	@property
	def total(self):
		return len(self.inputs) or 1

	def key(self, key):
		if key is None:
			return

		self.dirty = True
		match key:
			case 'n':
				self.next_input()
			case 'p':
				self.prev_input()
			case '-':
				self.width -= 1
			case '+':
				self.width += 1

	def next_input(self):
		self.file_no = (self.file_no + 1) % self.total
		self.page_no = 0
		self.dirty = True

	def prev_input(self):
		self.file_no = (self.file_no - 1) % self.total
		self.page_no = 0
		self.dirty = True


# ruff: noqa: PLR0913 # Too many arguments in function definition
def page(
	term,
	state,
	*,
	mode='RGB',
	amino=False, degen=False, table=1,
):
	width = state.width
	height = term.height - 1
	filename = state.current_input

	amino = auto_select_amino(filename, amino)
	seqs = read_sequences(filename, amino)

	for sequence in seqs:
		[descriptions, data] = descriptions_and_data(sequence)
		decoded = dna.decode(data, amino, table, degen)
		pixels = decoded_to_pixels(decoded, mode, amino, degen)
		yield from pixels_to_blocks(pixels, width=width, height=height, mode=mode)


def status(term, state):
	left_txt = f'file ({state.file_no + 1} / {state.total}): {state.current_input}'
	right_txt = str_join([
		f'term {term.width}x{term.height}',
		f'width {state.width}',
		# f'{term.number_of_colors} colors',
		# '?: help',
	], '; ')
	return (
		term.normal +
		term.white_on_purple + term.clear_eol +
		left_txt +
		term.rjust(right_txt, term.width - len(left_txt)) +
		term.normal
	)


# ruff: noqa: PLR0913 # Too many arguments in function definition
def zigzag_blocks(
	inputs,
	*,
	width=80, mode='RGB',
	amino=False, degen=False, table=1,
):
	"""Browse DNA data as text paged into variable line widths."""
	term = blessed.Terminal()
	state = ZigzagState(
		inputs = [*map(str, check_paths(inputs))],
		width=width
	)

	with term.cbreak(), term.hidden_cursor(), term.fullscreen():
		state.dirty = True
		while True:
			if state.dirty:
				yield term.clear
				yield from page(
					term,
					state,
					mode=mode,
					amino=amino,
					degen=degen,
					table=table,
				)
				yield status(term, state)
				sys.stdout.flush()
				state.dirty = False

			if key := term.inkey(timeout=None):
				if key == 'q':
					break
				state.key(key)

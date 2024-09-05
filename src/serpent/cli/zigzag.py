from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import more_itertools as mit

from serpent.cli.flow import flow_blocks, verbose_flow_blocks
from serpent.fun import str_join
from serpent.io.fasta import auto_select_amino, descriptions_and_data, read_sequences


@dataclass
class ZigzagState:
	"""Dataclass for keeping track of the zigzag state."""

	inputs: list[str]
	height: int = 80
	width: int = 80
	dirty: bool = True
	file_no: int = 0
	position: int = 0
	# last_key: str = ''
	verbose: bool = False

	@property
	def current_input(self):
		return self.inputs[self.file_no]

	@property
	def max_position(self):
		return Path(self.current_input).stat().st_size - self.screen_size

	@property
	def screen_size(self):
		# TODO Handle height with verbose
		return self.width * self.height

	@property
	def total(self):
		return len(self.inputs) or 1

	def change_position(self, change: int):
		self.position += change
		if change < 0:
			self.position = max(self.position, 0)
		else:
			self.position = min(self.position, self.max_position)

	# ruff: noqa: PLR0912: Too many branches
	def key(self, key):
		dirty_keys = ' bwsadADzxnpm-,.;:v'

		if key in dirty_keys:
			self.dirty = True

		match key:
			case ' ':
				self.change_position(self.screen_size * 3)
			case 'b':
				self.change_position(self.screen_size * -3)
			case 'w':
				self.change_position(self.width * 3)
			case 's':
				self.change_position(self.width * -3)
			case 'a':
				self.change_position(3)
			case 'd':
				self.change_position(-3)
			case 'A':
				self.change_position(self.width)
			case 'D':
				self.change_position(-self.width)
			# shift frames
			case 'z':
				self.change_position(1)
			case 'x':
				self.change_position(-1)
			case 'n':
				self.next_input()
			case 'p':
				self.prev_input()
			case ',':
				self.width -= 1
			case '.':
				self.width += 1
			case 'm':
				self.width -= 3
			case '-':
				self.width += 3
			case ';':
				self.width -= 9
			case ':':
				self.width += 9
			case 'v':
				self.verbose = not self.verbose

	def next_input(self):
		self.file_no = (self.file_no + 1) % self.total
		self.position = 0
		self.dirty = True

	def prev_input(self):
		self.file_no = (self.file_no - 1) % self.total
		self.position = 0
		self.dirty = True


# ruff: noqa: PLR0913 # Too many arguments in function definition
def zigzag_page(
	term,
	state,
	*,
	mode='RGB', fmt=None,
	amino=False, degen=False, table=1,
):
	width = state.width
	height = state.height = term.height
	filename = state.current_input

	amino = auto_select_amino(filename, amino)
	seqs = read_sequences(filename, amino)

	for sequence in seqs:
		[descriptions, data] = descriptions_and_data(sequence)

		data = mit.seekable(data)
		data.seek(state.position)
		data = mit.take(state.screen_size, data)

		if state.verbose:
			yield from descriptions
			yield from verbose_flow_blocks(
				data, width, mode, fmt,
				height=height,
				amino=amino, degen=degen, table=table)
		else:
			yield from flow_blocks(
				data, width, mode,
				height=height,
				amino=amino, degen=degen, table=table)


def zigzag_status(term, state):
	left_txt = f'file ({state.file_no + 1}/{state.total}): {state.current_input}'
	right_txt = str_join([
		f'term {term.width}x{term.height}',
		f'width {state.width}',
		f'height {state.height}',
		f'pos {state.position}',
		# f'key {state.last_key}',
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

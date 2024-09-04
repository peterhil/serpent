from __future__ import annotations

from dataclasses import dataclass

from serpent.cli.flow import flow_blocks, verbose_flow_blocks
from serpent.fun import str_join
from serpent.io.fasta import auto_select_amino, descriptions_and_data, read_sequences


@dataclass
class ZigzagState:
	"""Dataclass for keeping track of the zigzag state."""

	inputs: list[str]
	width: int = 80
	dirty: bool = True
	file_no: int = 0
	page_no: int = 0
	# last_key: str = ''
	verbose: bool = False

	@property
	def current_input(self):
		return self.inputs[self.file_no]

	@property
	def total(self):
		return len(self.inputs) or 1

	def key(self, key):
		if key is None:
			return
		# else:
		# 	self.last_key = key

		self.dirty = True
		match key:
			case 'n':
				self.next_input()
			case 'p':
				self.prev_input()
			case ',':
				self.width -= 1
			case '.':
				self.width += 1
			case ';':
				self.width -= 9
			case ':':
				self.width += 9
			case 'v':
				self.verbose = not self.verbose

	def next_input(self):
		self.file_no = (self.file_no + 1) % self.total
		self.page_no = 0
		self.dirty = True

	def prev_input(self):
		self.file_no = (self.file_no - 1) % self.total
		self.page_no = 0
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
	height = term.height - 1
	filename = state.current_input

	amino = auto_select_amino(filename, amino)
	seqs = read_sequences(filename, amino)

	for sequence in seqs:
		[descriptions, data] = descriptions_and_data(sequence)

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

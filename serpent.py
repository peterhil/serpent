import matplotlib.pyplot as plt
import numpy as np
import re
import sys

from PIL import Image
from collections import Counter
from itertools import combinations
from more_itertools import chunked
from pprint import pp
from scipy.fft import fft
from typing import NamedTuple
from warnings import warn


COUNT_LIMIT = 32
LINE_FEED = False
PLOT = True


def char_range(c1, c2):
	"""Generates the characters from `c1` to `c2`, inclusive."""
	for c in range(ord(c1), ord(c2)+1):
		yield chr(c)


def map_array(fn, arr, dtype=None):
	return np.array(list(map(fn, arr)), dtype=dtype)


def normalise(array):
	return array / np.amax(array)


bases = {'A': 0b00, 'C': 0b01, 'G': 0b10, 'T': 0b11, 'U': 0b11}
bases_inverse = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}

chars64= list(char_range('A', 'Z')) + \
	list(char_range('a', 'z')) + \
	list(char_range('0', '9')) + \
	(['+', '\n'] if LINE_FEED else ['+', '.'])

alphabet64 = dict({(i, char) for (i, char) in enumerate(chars64)})
combos = map_array(lambda cm: ''.join(cm), combinations(chars64, 2))


def decode(dna):
	"""Returns dna’s codons encoded into numbers 0..63"""
	return map_array(decode_codon, dna)


def decode_codon(codon):
	result = 0
	for num, char in enumerate(reversed(codon)):
		result += bases.get(char, 0) << num * 2 # Throw IndexError by using []?

	return result


def get_padding(data, divisor=3, fill='A'):
	"""Return suitable padding in order to pad the data into a length evenly divisible by the divisor."""
	padding = []
	rem = len(data) % divisor
	if rem != 0:
		pad_length = divisor - rem
		padding = pad_length * [fill]

	return padding


def pad_to_left(data, divisor=3, fill='A'):
	"""Pad data to a length divisible by the divisor with the fill characters on the end"""
	padding = get_padding(data, divisor, fill)

	return list(data) + padding


def pad_to_right(data, divisor=3, fill='A'):
	"""Pad data to a length divisible by the divisor with the fill character on the beginning"""
	padding = get_padding(data, divisor, fill)

	return padding + list(data)


def clean_non_dna(data):
	"""Clean up non DNA or RNA data. Warns if the data is in multiple parts."""
	cleaned = ''.join(re.sub(r'[^CGAT]{6,}', '', data).split('\n'))
	residual = re.sub(r'[CGAT]{6,}', '', cleaned)

	if len(residual) > 1:
		print('Residual characters:', residual)
		# raise UserWarning(f"Data has { len(residual) } extra characters, please check it carefully!")

	return cleaned


def get_codons(data):
	"""Get codons from data as Numpy array"""
	codons_list = list(chunked(data, 3, strict=True))
	codons = map_array(lambda c: ''.join(c), codons_list, dtype='U3')

	return codons


def digits_to_number(seq, base=64):
	value = 0
	for i, n in enumerate(reversed(seq)):
		value = value + base ** i * n

	return value


def number_to_digits(number, base=64):
	result = []
	rem = number

	while True:
		[mul, rem] = divmod(rem, base)
		result.append(rem)
		rem = mul
		if mul == 0:
			break

	return list(reversed(result))


def plot(data):
	plt.interactive(True)
	plt.show()

	return plt.plot(data)


def codon_sequences(decoded, n=4, fill=0):
	"""Chunk data into length N sequences of codons.
	Count the occurences of different kmers as numbers between 0..64**n.

	Return index and counts.
	"""
	padded = pad_to_left(list(decoded), n, fill)
	sequences = list(chunked(padded, n, strict=True))
	numbers = np.apply_along_axis(digits_to_number, 1, sequences)

	return numbers


def count_sorted(items):
	counts = Counter(items)
	return np.array(list(sorted(counts.items()))).T


def plot_sequence_counts(decoded, n=4, *args, **kwargs):
	numbers = codon_sequences(decoded, n)
	[index, count] = count_sorted(numbers)

	size = 64 ** n
	data = np.zeros(size, dtype=np.uint64)
	data[index] = count

	plt.plot(np.arange(size), data, *args, **kwargs)

	return [index, count]


def show_image(decoded, width=64, fill=0, mode='RGB'):
	padded = np.array(pad_to_left(decoded, 3 * width, fill))
	norm = padded / np.amax(padded)
	channels = len(mode)

	rows = len(norm) / (channels * width)
	height = int(np.floor(rows))

	if channels > 1:
		rgb = norm.reshape(height, width, channels)
	else:
		rgb = norm.reshape(height, width)
	img = Image.fromarray(np.uint8(rgb * 255), mode=mode)
	img.show()

	return img


def main(data, fn=None):
	data = clean_non_dna(data)
	data = pad_to_left(data, 3, 'A')

	codons = get_codons(data)
	print("Codons:\n", codons)

	counts = Counter(codons)
	print("Counts:")
	pp(dict(counts.most_common()[:COUNT_LIMIT]))

	decoded = decode(codons)
	print("Decoded:\n", len(np.unique(decoded)), decoded)
	if PLOT:
		plt.interactive(True)
		plt.show()

		# counts = Counter(decoded)
		# [index, count] = list(np.array(sorted(list(counts.items()))).T)
		# plt.plot(index, count)

		# ft = np.abs(fft(decoded, n=64, norm='ortho'))
		# plot(ft)

		# TODO Make subcommand
		seq_length = 1
		plot_sequence_counts(decoded, seq_length)

	encoded = ''.join([alphabet64.get(c, ' ') for c in decoded])
	print("Encoded:\n", encoded)

	# Write out base64 encoded data
	if fn:
		with open(fn + '.ser64', 'w', encoding='UTF-8') as file:
			file.write(''.join(map_array(str, encoded)))

	counts = Counter(encoded)
	print("Counts:")
	pp(dict(counts.most_common()[:COUNT_LIMIT]))

	# Bigrams:
	ch_list = list(chunked(encoded, 2))
	ch = map_array(lambda c: ''.join(c), ch_list, dtype='U2')
	counts = Counter(ch)
	print("Bigrams:\n")
	pp(dict(counts.most_common()[:COUNT_LIMIT]))

	print("Bigrams not appearing:\n")
	bigrams = combos[[cmb not in ch for cmb in combos]]
	print(bigrams)

	# Print codon sequences
	print("Codon sequences:\n")
	seq_length = 4
	[index, count] = count_sorted(codon_sequences(decoded, seq_length))
	twice_i = index[count == 2]
	codes = map_array(lambda a: pad_to_right(number_to_digits(a, 64), seq_length, 0), twice_i)
	b64_codes = map_array(lambda a: ''.join(list(map(alphabet64.get, number_to_digits(a)))), twice_i)
	print(b64_codes)

	catg = map_array(
		lambda a: ''.join(list(map(bases_inverse.get, pad_to_right(number_to_digits(a, 4), 3, 0)))),
		codes.flatten()
	)

	catg = catg.reshape(int(len(catg) / seq_length), seq_length)
	catg = map_array(lambda row: ''.join(row), catg)
	print(catg)

	show_image(decoded, width=64, fill=63, mode='RGB')

	return decoded


class Token(NamedTuple):
    type: str
    value: str
    line: int
    column: int


DATA_TOKENS = ['AMINO', 'BASE', 'DEGENERATE']


def tokenize(data, amino=False):
	"""Iterative FASTA sequence reader.
	See: https://docs.python.org/3/library/re.html#writing-a-tokenizer
	"""
	BASE = r'[ACGTU\n]'
	token_specification = [
		('DESCRIPTION',	 r'>[^\n]*'),
		('BASE', BASE + r'+'),
		('DEGENERATE', r'[WSMKRYBDHVNZ-]+?'),  # https://en.wikipedia.org/wiki/Nucleic_acid_sequence#Notation
		('NEWLINE',	 r'\n'),		   # Line endings
		('SKIP',	 r'[ \t]+'),	   # Skip over spaces and tabs
		('MISMATCH', r'.'),			   # Any other character
	]
	if amino:
		# https://en.wikipedia.org/wiki/Proteinogenic_amino_acid
		token_specification.insert(1, ('AMINO', r'[ARNDCQEGHILKMFPSTWYVUO]+'))

	tok_regex = '|'.join('(?P<%s>%s)' % pair for pair in token_specification)
	line_num = 1
	line_start = 0
	for mo in re.finditer(tok_regex, data):
		kind = mo.lastgroup
		value = mo.group()
		column = mo.start() - line_start
		if kind == 'NEWLINE':
			line_start = mo.end()
			line_num += 1
			continue
		elif kind == 'SKIP':
			continue
		elif kind == 'MISMATCH':
			raise RuntimeError(f'{value!r} unexpected on line {line_num}')
		yield Token(kind, value, line_num, column)


if __name__ == '__main__':
	args = sys.argv
	if len(args) < 2: print('Give a filename for DNA data.')
	fn = args[1]
	amino = '-a' in args
	writeout = '-o' in args

	data = []
	description_count = 0
	with open(fn, 'r', encoding='UTF-8') as file:
		while (line := file.readline().rstrip()) and description_count < 2:
			for token in tokenize(line, amino):
				# TODO Create a TUI or add CLI option to select sequences or
				# otherwise handle multiple sequences
				if token.type == 'DESCRIPTION':
					description_count += 1
				if description_count == 2:
					break
				if token.type in DATA_TOKENS:
					data.append(token.value)

	data = '\n'.join(data)

	outfile = fn if writeout else None
	decoded = main(data, fn=outfile)

	if description_count >= 2:
		print('Warning: File has more than one FASTA sequence!')

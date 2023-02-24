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


def plot_sequence_counts(decoded, n=4):
	numbers = codon_sequences(decoded, n)
	[index, count] = count_sorted(numbers)

	plt.plot(index, count)

	return [index, count]


def plot_sequence_repeats(decoded, n=4):
	numbers = codon_sequences(decoded, n)
	[index, count] = count_sorted(numbers)

	more_than_one_i = index[count > 1]
	more_than_one_c = count[count > 1]

	plt.plot(more_than_one_i, more_than_one_c)

	return [more_than_one_i, more_than_one_c]


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


def main(data):
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
		plot_sequence_counts(decoded, 4)

	encoded = ''.join([alphabet64.get(c, ' ') for c in decoded])
	print("Encoded:\n", encoded)

	counts = Counter(encoded)
	print("Counts:")
	pp(dict(counts.most_common()[:COUNT_LIMIT]))

	# Bigrams:
	ch_list = list(chunked(encoded, 2))
	ch = map_array(lambda c: ''.join(c), ch_list, dtype='U2')
	counts = Counter(ch)
	print("Bigrams:\m")
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

	show_image(decoded, width=108, fill=63, mode='RGB')


if __name__ == '__main__':
	args = sys.argv
	if len(args) < 2: print('Give a filename for DNA data.')
	fn = args[1]

	with open(fn, 'r', encoding='UTF-8') as file:
		lines = [line.rstrip() for line in file]
		# TODO Read and process files line by line
		# while (line := file.readline().rstrip()):
		# 	print(line)

	# TODO Use Regexp to extract data sequences OR
	# TODO Create iterative FASTA sequence reader OR
	# TODO Create a TUI to select a sequence or sequences
	data = '\n'.join([line for line in lines if line[0] != '>'])

	main(data)

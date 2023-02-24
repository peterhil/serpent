import matplotlib.pyplot as plt
import numpy as np
import re
import sys

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


bases = {'A': 0b00, 'C': 0b01, 'G': 0b10, 'T': 0b11, 'U': 0b11}
bases_inverse = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}

chars64= list(char_range('A', 'Z')) + \
	list(char_range('a', 'z')) + \
	list(char_range('0', '9')) + \
	(['+', '\n'] if LINE_FEED else ['+', '.'])

alphabet64 = dict({(i, char) for (i, char) in enumerate(chars64)})
combos = np.array(list(map(lambda cm: ''.join(cm), combinations(chars64, 2))))


def decode(dna):
	"""Returns dna’s codons encoded into numbers 0..63"""
	return np.array(list(map(decode_codon, dna)))


def decode_codon(codon):
	result = 0
	for num, char in enumerate(reversed(codon)):
		result += bases.get(char, 0) << num * 2 # Throw IndexError by using []?

	return result


def pad_to_left(data, divisor=3, fill='A'):
	"""Pad data to a length divisible by three by inserting fill character"""
	rem = len(data) % divisor
	if rem != 0:
		pad_length = divisor - rem
		data = data + pad_length * fill

	return data


def pad_to_right(data, divisor=3, fill='A'):
	"""Pad data to a length divisible by three by inserting fill character"""
	rem = len(data) % divisor
	if rem != 0:
		pad_length = divisor - rem
		data = pad_length * fill + data

	return data


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
	codons = np.array(list(map(lambda c: ''.join(c), codons_list)), dtype='U3')

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


if __name__ == '__main__':
	args = sys.argv
	if len(args) < 2: print('Give a filename for DNA data.')

	fn = args[1]
	with open(fn) as f:
		data = f.read().strip()

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
		counts = Counter(decoded)
		[index, count] = list(np.array(sorted(list(counts.items()))).T)
		plt.interactive(True)
		plt.plot(np.array(index + 1), count)
		plt.show()
		ft = np.abs(fft(decoded, n=64, norm='ortho'))
		plot(ft)

	encoded = ''.join([alphabet64.get(c, ' ') for c in decoded])
	print("Encoded:\n", encoded)

	counts = Counter(encoded)
	print("Counts:")
	pp(dict(counts.most_common()[:COUNT_LIMIT]))

	# Bigrams:
	ch_list = list(chunked(encoded, 2))
	ch = np.array(list(map(lambda c: ''.join(c), ch_list)), dtype='U2')
	counts = Counter(ch)
	print("Bigrams:\m")
	pp(dict(counts.most_common()[:COUNT_LIMIT]))

	print("Bigrams not appearing:\n")
	bigrams = combos[[cmb not in ch for cmb in combos]]
	print(bigrams)

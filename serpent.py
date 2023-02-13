import matplotlib.pyplot as plt
import numpy as np
import sys

from collections import Counter
from more_itertools import chunked
from pprint import pp
from scipy.fft import fft


LINE_FEED = False
PLOT = False


def char_range(c1, c2):
	"""Generates the characters from `c1` to `c2`, inclusive."""
	for c in range(ord(c1), ord(c2)+1):
		yield chr(c)


bases = {'A': 0b00, 'C': 0b01, 'G': 0b10, 'T': 0b11, 'U': 0b11}

chars64= list(char_range('A', 'Z')) + \
	list(char_range('a', 'z')) + \
	list(char_range('0', '9')) + \
	(['+', '\n'] if LINE_FEED else ['+', '.'])

alphabet64 = dict({(i, char) for (i, char) in enumerate(chars64)})


def decode(dna):
	"""Returns dna’s codons encoded into numbers 0..63"""
	return np.array(list(map(decode_codon, dna)))


def decode_codon(codon):
	result = 0
	for num, char in enumerate(reversed(codon)):
		result += bases.get(char, 0) << num * 2 # Throw IndexError by using []?

	return result


def pad_to_multiple(data, divisor=3, fill='A'):
	"""Pad data to a length divisible by three by inserting fill character"""
	pad_length = divisor - len(data) % divisor
	data = data + pad_length * fill

	return data


def get_codons(data):
	"""Get codons from data as Numpy array"""
	codons_list = list(chunked(data, 3, strict=True))
	codons = np.array(list(map(lambda c: ''.join(c), codons_list)), dtype='U3')

	return codons


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

	data = pad_to_multiple(data, 3, 'A')

	codons = get_codons(data)
	print("Codons:\n", codons)

	counts = Counter(codons)
	print("Counts:")
	pp(dict(counts.most_common()))

	decoded = decode(codons)
	print("Decoded:\n", decoded)
	if PLOT:
		ft = np.abs(fft(decoded, n=64, norm='ortho'))
		plot(ft)

	encoded = ''.join([alphabet64.get(c, ' ') for c in decoded])
	print("Encoded:\n", encoded)

	counts = Counter(encoded)
	print("Counts:")
	pp(dict(counts.most_common()))

#!/usr/bin/env python

import sys

import matplotlib.pyplot as plt
import numpy as np

from collections import Counter
from pprint import pp

from PIL import Image
from more_itertools import chunked
# from numpy.fft import fft

from serpent import dna
from serpent.encoding import alphabet64, combos
from serpent.fasta import read
from serpent.fun import map_array
from serpent.digit import number_to_digits
from serpent.mathematics import normalise
from serpent.padding import pad_to_left, pad_to_right
from serpent.stats import count_sorted


COUNT_LIMIT = 32
PLOT = True


def plot_sequence_counts(decoded, n=4, *args, **kwargs):
	numbers = dna.codon_sequences(decoded, n)
	[index, count] = count_sorted(numbers)

	size = 64**n
	data = np.zeros(size, dtype=np.uint64)
	data[index] = count

	plt.plot(np.arange(size), data, *args, **kwargs)

	return [index, count]


def show_image(decoded, width=64, fill=0, mode="RGB"):
	padded = np.array(pad_to_left(decoded, 3 * width, fill))
	norm = normalise(padded)
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


def analyse(data, fn=None):
	data = dna.clean_non_dna(data)
	data = pad_to_left(data, 3, "A")

	codons = dna.get_codons(data)
	print("Codons:\n", codons)

	counts = Counter(codons)
	print("Counts:")
	pp(dict(counts.most_common()[:COUNT_LIMIT]))

	decoded = dna.decode(codons)
	print("Decoded:\n", len(np.unique(decoded)), decoded)
	if PLOT:
		plt.interactive(True)
		plt.show()

		# TODO Make subcommand
		# ft = np.abs(fft(decoded, n=64, norm='ortho'))
		# plt.plot(ft)

		# TODO Make subcommand
		seq_length = 1
		plot_sequence_counts(decoded, seq_length)

	encoded = "".join([alphabet64.get(c, " ") for c in decoded])
	print("Encoded:\n", encoded)

	# Write out base64 encoded data
	if fn:
		with open(fn + ".ser64", "w", encoding="UTF-8") as file:
			file.write("".join(map_array(str, encoded)))

	counts = Counter(encoded)
	print("Counts:")
	pp(dict(counts.most_common()[:COUNT_LIMIT]))

	# Bigrams:
	ch_list = list(chunked(encoded, 2))
	ch = map_array(lambda c: "".join(c), ch_list, dtype="U2")
	counts = Counter(ch)
	print("Bigrams:\n")
	pp(dict(counts.most_common()[:COUNT_LIMIT]))

	print("Bigrams not appearing:\n")
	bigrams = combos[[cmb not in ch for cmb in combos]]
	print(bigrams)

	# Print codon sequences
	print("Codon sequences:\n")
	seq_length = 4
	[index, count] = count_sorted(dna.codon_sequences(decoded, seq_length))
	twice_i = index[count == 2]
	codes = map_array(
		lambda a: pad_to_right(number_to_digits(a, 64), seq_length, 0), twice_i
	)
	b64_codes = map_array(
		lambda a: "".join(list(map(alphabet64.get, number_to_digits(a)))), twice_i
	)
	print(b64_codes)

	catg = map_array(
		lambda a: "".join(
			list(map(dna.bases_inverse.get, pad_to_right(number_to_digits(a, 4), 3, 0)))
		),
		codes.flatten(),
	)

	catg = catg.reshape(int(len(catg) / seq_length), seq_length)
	catg = map_array(lambda row: "".join(row), catg)
	print(catg)

	show_image(decoded, width=64, fill=63, mode="RGB")

	return decoded


def main(args=None):
	# Ensure same behavior while testing and using the CLI
	args = args or sys.argv[1:]
	if len(args) == 0:
		print("Give a filename for DNA data.")
		sys.exit(1)
	fn = args[0]
	amino = "-a" in args
	writeout = "-o" in args

	data = read(fn, amino)

	outfile = fn if writeout else None
	decoded = analyse(data, fn=outfile)

	return decoded


if __name__ == "__main__":
	decoded = main()

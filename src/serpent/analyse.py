#!/usr/bin/env python

import sys

import matplotlib.pyplot as plt
import numpy as np

from collections import Counter
from pprint import pp

from more_itertools import chunked

from serpent import dna
from serpent.encoding import alphabet64, combos
from serpent.fasta import read
from serpent.fun import map_array
from serpent.digit import number_to_digits
from serpent.mathematics import logn
from serpent.padding import pad_to_left, pad_to_right
from serpent.stats import count_sorted
from serpent.visual import interactive, plot_fft, plot_histogram_sized, plot_sequence_counts, show_image


COUNT_LIMIT = 20
PLOT = True


def analyse(data, fn=None):
	data = dna.clean_non_dna(data)

	# Codons
	codons = dna.get_codons(data)
	counts = Counter(codons)
	print("Codons used:\n", np.unique(codons))
	print("Codons total:", len(codons), "unique:", len(np.unique(codons)))
	print("Counts:")
	pp(dict(counts.most_common()[:COUNT_LIMIT]))

	# Decode
	decoded = dna.decode(codons)
	print("Decoded:\n", len(np.unique(decoded)), decoded)

	# TODO Make subcommands
	if PLOT:
		interactive()
		# plot_fft(decoded, n=64)

		seq_length = 5
		seqs = dna.codon_sequences(decoded, seq_length)

		plot_histogram_sized(
			seqs,
			size='base',
			multi=max(16, 2 ** seq_length),  # cap to 16 * base = 1024
			density=False,
			cumulative=False,
		)
		# plot_sequence_counts(decoded, n=seq_length)

		# show_image(decoded, width=64, fill=63, mode="RGB")
	else:
		# Encode
		encoded = "".join([alphabet64.get(c, " ") for c in decoded])
		print("Encoded:\n", encoded)
		counts = Counter(encoded)
		print("Counts:")
		pp(dict(counts.most_common()[:COUNT_LIMIT]))

		# Write out base64 encoded data
		if fn:
			with open(fn + ".ser64", "w", encoding="UTF-8") as file:
				file.write("".join(map_array(str, encoded)))

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
			lambda a: pad_to_right(number_to_digits(a, 64), seq_length, 0),
			twice_i
		)
		b64_codes = map_array(
			lambda a: "".join(list(map(
				alphabet64.get,
				number_to_digits(a)))),
			twice_i
		)
		print(b64_codes)
		catg = map_array(
			lambda a: "".join(list(map(
				dna.bases_inverse.get,
				pad_to_right(number_to_digits(a, 4), 3, 0)))),
			codes.flatten(),
		)
		catg = catg.reshape(int(len(catg) / seq_length), seq_length)
		catg = map_array(lambda row: "".join(row), catg)
		print(catg)

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

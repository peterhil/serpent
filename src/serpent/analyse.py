#!/usr/bin/env python
"""Serpent analysis."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from pprint import pp

import argh
import numpy as np
from argh.decorators import arg
from more_itertools import chunked

from serpent import dna
from serpent.digit import number_to_digits
from serpent.encoding import alphabet64, combos
from serpent.fasta import read
from serpent.fun import map_array, str_join
from serpent.padding import pad_to_right
from serpent.stats import count_sorted
from serpent.visual import (
	dna_image,
	interactive,
	plot_fft,
	plot_histogram_sized,
)

COUNT_LIMIT = 20


def format_lines(data, width=80, sep=' '):
	"""Format lines for printing."""
	def fmt(chunk):
		return str_join(chunk, sep)
	lines = [
		f"{i * width}:\t{fmt(chunk)}"
		for i, chunk
		in enumerate(chunked(data, width))
	]

	return lines


def decode(data, verbose=False):
	data = dna.clean_non_dna(data)
	codons = dna.get_codons(data)
	decoded = dna.decode(codons)

	if verbose:
		uniq = len(np.unique(decoded))
		print(f"Decoded ({len(decoded)}, unique: {uniq}):")
		strings = map(str, iter(decoded))
		lines = format_lines(strings, 32)
		{print(line) for line in lines}

	return decoded


def analyse(decoded, plot=False, filename=None):
	"""Analyse data."""
	# TODO Make subcommands
	if plot:
		interactive()
		plot_fft(decoded, n=238)

		seq_length = 3
		seqs = dna.codon_sequences(decoded, seq_length)

		plot_histogram_sized(
			seqs,
			size='base',
			multi=max(16, 2 ** seq_length),  # cap to 16 * base = 1024
			density=False,
			cumulative=False,
		)
		# plot_sequence_counts(decoded, n=seq_length)

		# width = 64
		width = 119
		# width = 64 ** 2 - 1
		img = dna_image(decoded, width=width, fill=63, mode="RGB")
		img.show()

		if filename:
			img.save(filename + f".w{width}.png")
	else:
		# Encode
		encoded = str_join([alphabet64.get(c, " ") for c in decoded])
		print("Encoded:\n", encoded)
		counts = Counter(encoded)
		print("Counts:")
		pp(dict(counts.most_common()[:COUNT_LIMIT]))

		# Write out base64 encoded data
		if filename:
			with Path(filename + ".ser64").open("w", encoding="UTF-8") as file:
				file.write(str_join(map_array(str, encoded)))

		# Bigrams:
		ngram_list = list(chunked(encoded, 2))
		ngrams = map_array(str_join, ngram_list, dtype="U2")
		counts = Counter(ngrams)
		print("Bigrams:\n")
		pp(dict(counts.most_common()[:COUNT_LIMIT]))
		print("Bigrams not appearing:\n")
		bigrams = combos[[combo not in ngrams for combo in combos]]
		print(bigrams)

		# Print codon sequences
		print("Codon sequences:\n")
		seq_length = 4
		occurences = 2
		[index, count] = count_sorted(dna.codon_sequences(decoded, seq_length))
		twice_i = index[count == occurences]
		codes = map_array(
			lambda a: pad_to_right(number_to_digits(a, 64), fill=0, n=seq_length),
			twice_i
		)
		b64_codes = map_array(
			lambda a: str_join(list(map(
				alphabet64.get,
				number_to_digits(a)))),
			twice_i
		)
		print(b64_codes)
		catg = map_array(
			lambda a: str_join(list(map(
				dna.bases_inverse.get,
				pad_to_right(number_to_digits(a, 4), fill=0, n=3)))),
			codes.flatten(),
		)
		catg = catg.reshape(int(len(catg) / seq_length), seq_length)
		catg = map_array(str_join, catg)
		print(catg)

	return decoded


@arg('--stats',  '-s', help='Show statistics')
@arg('--width',  '-w', help='Codons per line')
def codons(filename, width=20, stats=False):
	"""Print codons and statistics."""
	data = read(filename)
	data = dna.clean_non_dna(data)

	codons = dna.get_codons(data)

	if stats:
		counts = Counter(codons)
		print("Codons used:\n", np.unique(codons))
		print("Codons total:", len(codons), "unique:", len(np.unique(codons)))
		print("Counts:")
		pp(dict(counts.most_common()[:COUNT_LIMIT]))
	else:
		print("Codons:")
		lines = format_lines(codons, width)
		{print(line) for line in lines}

	# return codons


@arg('--amino',    '-a', help='Read input as amino acids')
@arg('--plot',     '-p', help='Use plotting')
@arg('--verbose',  '-v', help='Verbose')
@arg('--writeout', '-w', help='Write base 64 encoded data out')
def serpent(filename, amino=False, plot=False, verbose=False, writeout=False):
	"""Explore DNA data with Serpent."""
	data = read(filename, amino)
	outfile = filename if writeout else None
	decoded = decode(data, verbose)

	analyse(decoded, plot, filename=outfile)

	return decoded


def main():
	parser = argh.ArghParser()
	parser.add_commands([
		codons,
		serpent,
	])
	# parser.set_default_command(serpent)
	return parser.dispatch()


if __name__ == '__main__':
	decoded = main()

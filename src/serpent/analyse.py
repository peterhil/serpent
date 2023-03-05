#!/usr/bin/env python
"""Serpent analysis."""

from __future__ import annotations

import itertools as itr
from collections import Counter
from pathlib import Path
from pprint import pp

import argh
import numpy as np
from argh.decorators import arg
from more_itertools import chunked

from serpent import dna
from serpent.digit import number_to_digits
from serpent.encoding import alphabet64, base64
from serpent.fasta import read
from serpent.fun import map_array, str_join
from serpent.mathematics import autowidth, phi
from serpent.padding import pad_to_right
from serpent.printing import format_decoded, format_lines
from serpent.stats import count_sorted
from serpent.visual import (
	dna_image,
	interactive,
	plot_fft,
	plot_histogram_sized,
	plot_sequence_counts,
)

COUNT_LIMIT = 20
DEFAULT_COLOR = '#70e'  # TODO Add config file


def analyse(decoded):
	"""Analyse data."""
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


@arg('--amino', '-a', help='Read input as amino acids')
@arg('--count', '-c', help='Print out counts')
@arg('--out',   '-o', help='Write image to file')
@arg('--verbose',  '-v', help='Verbose')
def encode(filename, amino=False, count=False, out=False, verbose=False):
	"""Encode data into various formats.

	Base 64: Uses characters A-Za-z0-9+. for numbers 0..63.
	"""
	data = read(filename, amino)
	decoded = dna.decode(data, amino)

	encoded = str_join([alphabet64.get(c, " ") for c in decoded])
	if verbose:
		print("Encoded:")
	print(encoded)

	if count:
		counts = Counter(encoded)
		if verbose:
			print("Counts:")
		pp(dict(counts.most_common()[:COUNT_LIMIT]))

	# Write out encoded data
	if out:
		with Path(filename + ".ser64").open("w", encoding="UTF-8") as file:
			file.write(str_join(map_array(str, encoded)))


@arg('--amino', '-a', help='Read input as amino acids')
@arg('--mode',  '-m', help='Image mode', choices=('RGB', 'L'))
@arg('--out',   '-o', help='Write image to file')
@arg('--width', '-w', help='Image width', type=int)
def image(filename, amino=False, width=None, mode="RGB", out=False):
	data = read(filename, amino)
	decoded = dna.decode(data, amino)

	if not width:
		width = autowidth(len(decoded) / len(mode), aspect=phi-1, base=64)

	img = dna_image(decoded, width=width, fill=63, mode=mode)
	img.show()

	if out:
		img.save(filename + f".w{width}.png")


@arg('--amino', '-a', help='Read input as amino acids')
@arg('--length', '-l', help='Fourier transform length', type=int)
def fft(filename, amino=False, length=64):
	"""Plot Fourier transform of the DNA data."""
	data = read(filename, amino)
	decoded = dna.decode(data, amino)

	interactive()
	plot_fft(decoded, n=length, color=DEFAULT_COLOR)


@arg('--amino', '-a', help='Read input as amino acids')
@arg('--cumulative', '-c', help='Cumulative distribution')
@arg('--density', '-d', help='Normalise histogram into density')
@arg('--length', '-l', help='Sequence length', type=int)
def hist(
	filename,
	amino=False,
	length=1,
	density=False,
	cumulative=False,
):
	"""Plot DNA data histograms."""
	data = read(filename, amino)
	decoded = dna.decode(data, amino)
	seqs = dna.codon_sequences(decoded, length)

	interactive()
	plot_histogram_sized(
		seqs,
		size='base',
		multi=max(16, 2 ** length),  # cap to 16 * base = 1024
		color=DEFAULT_COLOR,
		cumulative=cumulative,
		density=density,
	)


@arg('--amino', '-a', help='Read input as amino acids')
@arg('--length', '-l', help='Sequence length', type=int)
def seq(filename, amino=False, length=1):
	"""Plot DNA sequence count statistics."""
	data = read(filename, amino)
	decoded = dna.decode(data, amino)

	interactive()
	plot_sequence_counts(decoded, n=length, color=DEFAULT_COLOR)


@arg('--amino',   '-a', help='Read input as amino acids')
@arg('--length',  '-l', help='Peptide length', type=int)
@arg('--missing', '-m', help='Missing peptides')
def pep(filename, amino=False, length=2, missing=False):
	"""Peptide statistics."""
	data = read(filename, amino)
	decoded = dna.decode(data, amino)
	dtype = f'U{length}'

	# TODO Allow using amino acid codes?
	encoded = str_join([alphabet64.get(c, " ") for c in decoded])

	# TODO Use iterators only if Counter accepts them
	peptide_list = list(chunked(encoded, length))
	peptides = map_array(str_join, peptide_list, dtype=dtype)

	if not missing:
		print("Peptides:")
		lines = format_lines(peptides, 32)
		{print(line) for line in lines}

		print()
		print("Counts:")
		counts = Counter(peptides)
		# groups = itr.groupby(counts.most_common(), lambda x: x[1])
		# {print(f"{count}:\t{list(item[0])}") for [count, items] in groups}

		# TODO There must be batter way to use these iterators!
		grouped = [
			(c, [k for k, c in group])
			for c, group in itr.groupby(counts.most_common(), lambda x: x[1])
			if c > 1
		]
		for count, values in grouped:
			lines = format_lines(sorted(values), 32)
			print(f"-- {count} times --")
			{print(line) for line in lines}

	if missing:
		print("Peptides not appearing:\n")
		combos = np.fromiter(map(str_join, itr.product(base64, repeat=length)), dtype=dtype)
		absent = combos[[combo not in peptides for combo in combos]]

		lines = format_lines(absent, 64)
		{print(line) for line in lines}


@arg('--amino',    '-a', help='Read input as amino acids')
@arg('--verbose',  '-v', help='Verbose')
def serpent(filename, amino=False, verbose=False):
	"""Explore DNA data with Serpent."""
	data = read(filename, amino)
	decoded = dna.decode(data, amino)

	if verbose:
		lines = format_decoded(decoded)
		{print(line) for line in lines}

	analyse(decoded)

	return decoded


def main():
	parser = argh.ArghParser()
	parser.add_commands([
		codons,
		encode,
		fft,
		hist,
		image,
		pep,
		seq,
		serpent,
	])
	# parser.set_default_command(serpent)
	return parser.dispatch()


if __name__ == '__main__':
	decoded = main()

#!/usr/bin/env python
# flake8: noqa: PLR0913
"""Serpent analysis."""

from __future__ import annotations

import fileinput
import itertools as itr
from collections import Counter
from pathlib import Path
from pprint import pp

import argh
import matplotlib.pyplot as plt
import numpy as np
from argh.decorators import arg
from more_itertools import chunked

from serpent import dna
from serpent.config import COUNT_LIMIT, DEFAULT_COLOR
from serpent.convert.base64 import base64_to_num, num_to_base64
from serpent.digit import num_to_digits
from serpent.encoding import BASE64
from serpent.fasta import read
from serpent.fun import map_array, sort_values, str_join
from serpent.io import (
	check_inputs,
	echo,
	find_fasta_files,
	find_fasta_sequences,
	openhook,
)
from serpent.mathematics import autowidth, phi, phi_small
from serpent.padding import pad_start
from serpent.printing import format_decoded, format_lines
from serpent.stats import ac_peaks, autocorrelogram, count_sorted
from serpent.visual import (
	dna_image,
	interactive,
	plot_fft,
	plot_histogram_sized,
	plot_sequence_counts,
)


@arg('--amino',  '-a', help='Read input as amino acids')
@arg('--limit',  '-l', help='Peak limit', type=float)
@arg('--linear', '-n', help='Linear output (do not sort results)')
@arg('--seq',    '-s', help='Sequence length', type=int)
@arg('--width',  '-w', help='Autocorrelogram width', type=int)
def ac(filename, amino=False, limit=0.05, linear=False, width=256, seq=1):
	"""Plot autocorrelation of sequences."""
	data = read(filename, amino)
	data = dna.clean_non_dna(data, amino)

	decoded = dna.decode(data, amino)
	if seq > 1:
		decoded = dna.codon_sequences(decoded, seq)

	ac = autocorrelogram(decoded, width)

	interactive()
	plt.plot(ac, color=DEFAULT_COLOR)

	peaks = ac_peaks(ac, limit)
	peaks.pop(0)  # First item is always 1

	peaks = peaks.items() if linear else sort_values(peaks, reverse=True)
	return ("%s	%1.3f" % (peak, value) for peak, value in peaks)


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


def cat(*inputs):
	"""Concatenate and print FASTA sequences from files."""
	# TODO Implement read_sequences with more_itertools.split_with and chain
	# TODO Handle compressed files and archives
	inputs = check_inputs(inputs)
	with fileinput.input(inputs, mode='r', openhook=openhook) as fi:
		for line in fi:
			yield line.rstrip()


@arg('--seq', '-s', help='Print sequences')
def find(*inputs, seq=False):
	"""Find FASTA sequences from files and directories."""
	# TODO Handle compressed files and archives
	debug = False
	inputs = check_inputs(inputs, recurse=True)

	with fileinput.input(inputs, mode='r', openhook=openhook) as fi:
		if seq:
			yield from find_fasta_sequences(fi, debug)
		else:
			yield from find_fasta_files(fi, debug)


@arg('--amino', '-a', help='Read input as amino acids')
def decode(filename, amino=False):
	"""Explore DNA data with Serpent."""
	data = read(filename, amino)
	decoded = dna.decode(data, amino)

	lines = format_decoded(decoded)
	{print(line) for line in lines}


@arg('--amino', '-a', help='Read input as amino acids')
@arg('--count', '-c', help='Print counts')
@arg('--out',   '-o', help='Write out to file')
@arg('--width', '-w', help='Line width', type=int)
def encode(filename, amino=False, count=False, out=False, width=64):
	"""Encode data into various formats."""
	# TODO Read and decode data iteratively
	data = read(filename, amino)
	decoded = dna.decode(data, amino)

	encoded = (num_to_base64.get(c, " ") for c in decoded)

	if count:
		counts = Counter(encoded)
		return (f"{count}\t{codon}" for codon, count in counts.most_common())

	lines = (str_join(line) for line in chunked(encoded, width))
	if out:
		outfile = filename + ".ser64"
		with Path(outfile).open("w", encoding="UTF-8") as file:
			file.write(str_join(lines, "\n"))
			print(f"Wrote: {outfile}")
		return None
	else:
		return lines


@arg('--amino', '-a', help='Read input as amino acids')
@arg('--mode',  '-m', help='Image mode', choices=('RGB', 'L'))
@arg('--out',   '-o', help='Write image to file')
@arg('--width', '-w', help='Image width', type=int)
def image(filename, amino=False, width=None, mode="RGB", out=False):
	"""Visualise FASTA data as images."""
	data = read(filename, amino)
	decoded = dna.decode(data, amino)

	if not width:
		width = autowidth(len(decoded) / len(mode), aspect=phi-1, base=64)

	img = dna_image(decoded, width=width, fill=63, mode=mode)
	img.show()

	if out:
		img.save(filename + f".w{width}.{dna.BASE_ORDER}.png")


@arg('--amino', '-a', help='Read input as amino acids')
@arg('--length', '-l', help='Fourier transform length', type=int)
def fft(filename, amino=False, length=64):
	"""Plot Fourier transform of the DNA data."""
	data = read(filename, amino)
	decoded = dna.decode(data, amino)

	interactive()
	plot_fft(decoded, n=length, color=DEFAULT_COLOR)


@arg('--amino', '-a', help='Read input as amino acids')
@arg('--bins', '-b', help='Histogram bin sizing', choices=[
	'base', 'auto', 'fd', 'doane', 'scott', 'stone', 'rice', 'sturges', 'sqrt',
])
@arg('--num', '-n', help='Number of bins (overrides bins option)', type=int)
@arg('--cumulative', '-c', help='Cumulative distribution')
@arg('--density', '-d', help='Normalise histogram into density')
@arg('--length', '-l', help='Sequence length', type=int)
def hist(
	filename,
	amino=False,
	bins='base',
	num=None,
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
		size=num or bins,
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
	encoded = str_join([num_to_base64.get(c, " ") for c in decoded])

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
		combos = np.fromiter(map(str_join, itr.product(BASE64, repeat=length)), dtype=dtype)
		absent = combos[[combo not in peptides for combo in combos]]

		lines = format_lines(absent, 64)
		{print(line) for line in lines}


def analyse_repeats(decoded, length=4, limit=2, encode=False):
	"""Analyse codon sequence repeats."""
	[index, count] = count_sorted(dna.codon_sequences(decoded, length))
	repeats = index[count >= limit]

	# TODO Move conversions to helper functions
	echo("Repeated codon sequences:")
	if encode:
		b64_codes = map_array(
			lambda a: str_join(map(num_to_base64.get, num_to_digits(a))),
			repeats
		)

		lines = format_lines(b64_codes, 32)
		{print(line) for line in lines}
	else:
		codes = map_array(
			lambda a: pad_start(num_to_digits(a, 64), fill=0, n=length),
			repeats
		)
		catg = map_array(
			lambda a: str_join(map(
				dna.num_to_base.get,
				pad_start(num_to_digits(a, 4), fill=0, n=3))),
			codes.flatten(),
		)
		catg = catg.reshape(int(len(catg) / length), length)
		catg = map_array(str_join, catg)

		lines = format_lines(catg, 8)
		{print(line) for line in lines}


@arg('--amino',   '-a', help='Read input as amino acids')
@arg('--encode',  '-e', help='Encode output as base 64')
@arg('--limit',   '-l', help='Limit to at least this many repeats', type=int)
@arg('--seq',     '-s', help='Sequence length', type=int)
def repeats(filename, amino=False, seq=2, limit=2, encode=False):
	"""Find repeated sequences."""
	data = read(filename, amino)
	decoded = dna.decode(data, amino)

	analyse_repeats(decoded, length=seq, limit=limit, encode=encode)


def main():
	parser = argh.ArghParser()
	parser.add_commands([
		ac,
		cat,
		codons,
		decode,
		encode,
		fft,
		find,
		hist,
		image,
		pep,
		repeats,
		seq,
	])
	# parser.set_default_command(serpent)
	return parser.dispatch()


if __name__ == '__main__':
	decoded = main()

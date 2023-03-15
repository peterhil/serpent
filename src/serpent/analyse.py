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
from argh.decorators import arg, aliases, wrap_errors
from more_itertools import chunked

from serpent import dna
from serpent.convert.amino import aa_tables
from serpent.convert.base64 import base64_to_num, num_to_base64
from serpent.convert.digits import num_to_digits
from serpent.convert.nucleotide import num_to_nt
from serpent.encoding import BASE64
from serpent.fasta import (
	ParseError,
	auto_select_amino,
	find_fasta_files,
	find_fasta_sequences,
	read,
)
from serpent.fun import map_array, sort_values, str_join
from serpent.io import (
	check_paths,
	echo,
	openhook,
)
from serpent.mathematics import autowidth, phi, phi_small
from serpent.padding import pad_start
from serpent.printing import format_decoded, format_lines
from serpent.settings import COUNT_LIMIT, DEFAULT_COLOR
from serpent.stats import ac_peaks, autocorrelogram, count_sorted
from serpent.visual import (
	bin_choices,
	dna_image,
	interactive,
	plot_fft,
	plot_histogram_sized,
	plot_sequence_counts,
)
from serpent.zigzag import zigzag_blocks, zigzag_text


wrapped_errors = [AssertionError, ParseError]


@arg('--amino',  '-a', help='Amino acid input')
@arg('--table',  '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen',  '-g', help='Degenerate data')
@arg('--limit',  '-l', help='Peak limit', type=float)
@arg('--linear', '-n', help='Linear output (do not sort results)')
@arg('--seq',    '-s', help='Sequence length', type=int)
@arg('--width',  '-w', help='Autocorrelogram width', type=int)
@wrap_errors(wrapped_errors)
def ac(
	filename,
	limit=0.05, linear=False, width=256, seq=1,
	amino=False, degen=False, table=1,
):
	"""Plot autocorrelation of sequences."""
	amino = auto_select_amino(filename, amino)
	data = read(filename, amino)
	decoded = dna.decode(data, amino, table, degen)

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
	inputs = check_paths(inputs)

	with fileinput.input(inputs, mode='r', openhook=openhook) as fi:
		for line in fi:
			yield line.rstrip()


@arg('--seq', '-s', help='Print sequences')
def find(*inputs, seq=False):
	"""Find FASTA sequences from files and directories."""
	# TODO Handle compressed files and archives
	debug = False
	inputs = check_paths(inputs, recurse=True)

	with fileinput.input(inputs, mode='r', openhook=openhook) as fi:
		if seq:
			yield from find_fasta_sequences(fi, debug)
		else:
			yield from find_fasta_files(fi, debug)


@arg('--text', '-x', help='Page as text')
@aliases('zz')
@wrap_errors(wrapped_errors)
def zigzag(*inputs, text=False):
	"""Browse DNA data paged into variable line widths."""
	if text:
		zigzag_text(inputs)
	else:
		zigzag_blocks(inputs)


@arg('--amino', '-a', help='Amino acid input')
@arg('--table', '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen', '-g', help='Degenerate data')
@wrap_errors(wrapped_errors)
def decode(filename, amino=False, degen=False, table=1):
	"""Explore DNA data with Serpent."""
	amino = auto_select_amino(filename, amino)
	data = read(filename, amino)
	decoded = dna.decode(data, amino, table, degen)

	lines = format_decoded(decoded)
	{print(line) for line in lines}


@arg('--amino', '-a', help='Amino acid input')
@arg('--table', '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen', '-g', help='Degenerate data')
@arg('--count', '-c', help='Print counts')
@arg('--out',   '-o', help='Write out to file')
@arg('--width', '-w', help='Line width', type=int)
@wrap_errors(wrapped_errors)
def encode(
	filename,
	count=False, out=False, width=64,
	amino=False, degen=False, table=1,
):
	"""Encode data into various formats."""
	# TODO Read and decode data iteratively
	amino = auto_select_amino(filename, amino)
	data = read(filename, amino)
	decoded = dna.decode(data, amino, table, degen)

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


@arg('--amino', '-a', help='Amino acid input')
@arg('--table', '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen', '-g', help='Degenerate data')
@arg('--mode',  '-m', help='Image mode', choices=('RGB', 'L'))
@arg('--out',   '-o', help='Write image to file')
@arg('--width', '-w', help='Image width', type=int)
@wrap_errors(wrapped_errors)
def image(
	filename,
	width=None, mode="RGB", out=False,
	amino=False, degen=False, table=1,
):
	"""Visualise FASTA data as images."""
	amino = auto_select_amino(filename, amino)
	data = read(filename, amino)
	decoded = dna.decode(data, amino, table, degen)

	if not width:
		width = autowidth(len(decoded) / len(mode), aspect=phi-1, base=64)

	img = dna_image(decoded, width=width, fill=0, mode=mode)
	img.show()

	if out:
		if amino and table != 1:
			img.save(filename + f".w{width}.{dna.BASE_ORDER}.t{table}.png")
		else:
			img.save(filename + f".w{width}.{dna.BASE_ORDER}.png")


@arg('--amino',  '-a', help='Amino acid input')
@arg('--table',  '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen',  '-g', help='Degenerate data')
@arg('--length', '-l', help='Fourier transform length', type=int)
@wrap_errors(wrapped_errors)
def fft(
	filename,
	length=64,
	amino=False, degen=False, table=1,
):
	"""Plot Fourier transform of the DNA data."""
	amino = auto_select_amino(filename, amino)
	data = read(filename, amino)
	decoded = dna.decode(data, amino, table, degen)

	interactive()
	plot_fft(decoded, n=length, color=DEFAULT_COLOR)


@arg('--amino', '-a', help='Amino acid input')
@arg('--table', '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen', '-g', help='Degenerate data')
@arg('--bins',  '-b', help='Histogram bin sizing', choices=bin_choices)
@arg('--num',   '-n', help='Number of bins (overrides bins option)', type=int)
@arg('--cumulative', '-c', help='Cumulative distribution')
@arg('--density', '-d', help='Normalise histogram into density')
@arg('--length', '-l', help='Sequence length', type=int)
@wrap_errors(wrapped_errors)
def hist(
	filename,
	bins='base',
	num=None,
	length=1,
	density=False,
	cumulative=False,
	amino=False, degen=False, table=1,
):
	"""Plot DNA data histograms."""
	amino = auto_select_amino(filename, amino)
	data = read(filename, amino)
	decoded = dna.decode(data, amino, table, degen)
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


@arg('--amino',  '-a', help='Amino acid input')
@arg('--table',  '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen',  '-g', help='Degenerate data')
@arg('--length', '-l', help='Sequence length', type=int)
@wrap_errors(wrapped_errors)
def seq(filename, length=1, amino=False, degen=False, table=1):
	"""Plot DNA sequence count statistics."""
	amino = auto_select_amino(filename, amino)
	data = read(filename, amino)
	decoded = dna.decode(data, amino, table, degen)

	interactive()
	plot_sequence_counts(decoded, n=length, color=DEFAULT_COLOR)


@arg('--amino',   '-a', help='Amino acid input')
@arg('--table',   '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen',   '-g', help='Degenerate data')
@arg('--length',  '-l', help='Peptide length', type=int)
@arg('--missing', '-m', help='Missing peptides')
@wrap_errors(wrapped_errors)
def pep(
	filename,
	length=2, missing=False,
	amino=False, table=1, degen=False,
):
	"""Peptide statistics."""
	amino = auto_select_amino(filename, amino)
	data = read(filename, amino)
	decoded = dna.decode(data, amino, table, degen)
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
				num_to_nt.get,
				pad_start(num_to_digits(a, 4), fill=0, n=3))),
			codes.flatten(),
		)
		catg = catg.reshape(int(len(catg) / length), length)
		catg = map_array(str_join, catg)

		lines = format_lines(catg, 8)
		{print(line) for line in lines}


@arg('--amino',  '-a', help='Amino acid input')
@arg('--table',  '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen',  '-g', help='Degenerate data')
@arg('--encode', '-e', help='Encode output as base 64')
@arg('--limit',  '-l', help='Limit to at least this many repeats', type=int)
@arg('--seq',    '-s', help='Sequence length', type=int)
@wrap_errors(wrapped_errors)
def repeats(
	filename,
	seq=2, limit=2, encode=False,
	amino=False, degen=False, table=1,
):
	"""Find repeated sequences."""
	amino = auto_select_amino(filename, amino)
	data = read(filename, amino)
	decoded = dna.decode(data, amino, table, degen)

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
		zigzag,
	])
	# parser.set_default_command(serpent)
	return parser.dispatch()


if __name__ == '__main__':
	decoded = main()

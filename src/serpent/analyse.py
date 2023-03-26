#!/usr/bin/env python
# flake8: noqa: PLR0913
"""Serpent analysis."""

from __future__ import annotations

import fileinput
import itertools as itr
import os
import re
from PIL import Image
from collections import Counter
from pathlib import Path
from pprint import pp

import argh
import matplotlib.pyplot as plt
import numpy as np
from argh.decorators import arg, aliases, wrap_errors
from more_itertools import chunked, grouper, partition, split_before, take

from serpent import dna
from serpent.convert.amino import aa_tables
from serpent.convert.base64 import base64_to_num, num_to_base64
from serpent.convert.codon import num_to_codon
from serpent.convert.digits import num_to_digits
from serpent.convert.nucleotide import num_to_nt
from serpent.encoding import BASE64
from serpent.fasta import (
	ParseError,
	RE_DESCRIPTION,
	auto_select_amino,
	data_and_descriptions,
	find_fasta_files,
	find_fasta_sequences,
	read,
	read_sequences,
	read_tokens,
)
from serpent.fun import map_array, sort_values, str_join
from serpent.io import (
	check_paths,
	echo,
	file_extension_for,
	openhook,
)
from serpent.mathematics import autowidth, phi, phi_small
from serpent.padding import pad_start
from serpent.printing import format_counts, format_decoded, format_lines
from serpent.settings import COUNT_LIMIT, DEFAULT_COLOR
from serpent.stats import ac_peaks, autocorrelogram, count_sorted
from serpent.visual import (
	bin_choices,
	dna_image,
	dna_image_seq,
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
@arg('--limit',  '-l', help='Limit codon stats to N most common')
@arg('--width',  '-w', help='Codons per line')
def codons(filename, width=20, stats=False, limit=COUNT_LIMIT):
	"""Print codons and statistics."""
	seqs = read_sequences(filename)

	for seq in seqs:
		[tokens, descriptions] = data_and_descriptions(seq)

		yield from (token.value for token in descriptions)
		data = str_join(token.data for token in tokens)
		codons = dna.get_codons(data)

		if stats:
			unique = np.unique(codons)
			print(f"Unique {len(unique)} codons used: (total: {len(codons)})")
			yield from format_lines(unique, width)
			print("Counts:")
			counts = Counter(codons)
			yield from format_counts(counts, limit)
		else:
			lines = format_lines(codons, width)
			yield from lines


def cat(*inputs):
	"""Concatenate and print FASTA sequences from files."""
	# TODO Implement read_sequences with more_itertools.split_with and chain
	# TODO Handle compressed files and archives
	paths = check_paths(inputs)

	with fileinput.input(paths, mode='r', openhook=openhook) as fi:
		for line in fi:
			yield line.rstrip()


@arg('--num', '-n', help='Print line numbers')
@arg('--seq', '-s', help='Print sequences')
def find(*inputs, num=False, seq=False):
	"""Find FASTA sequences from files and directories."""
	# TODO Handle compressed files and archives
	debug = False
	paths = check_paths(inputs, recurse=True)

	with fileinput.input(paths, mode='r', openhook=openhook) as fi:
		if seq:
			yield from find_fasta_sequences(fi, num, debug)
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
	"""Decode DNA data into numbers."""
	amino = auto_select_amino(filename, amino)
	data = read(filename, amino)
	decoded = dna.decode(data, amino, table, degen)

	lines = format_decoded(decoded)
	{print(line) for line in lines}


@arg('--amino', '-a', help='Amino acid input')
@arg('--table', '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen', '-g', help='Degenerate data')
@arg('--count', '-c', help='Print counts')
@arg('--fmt',   '-f', help='Output format', choices=['b', 'base64', 'c', 'codon'])
@arg('--out',   '-o', help='Write out to file')
@arg('--width', '-w', help='Line width', type=int)
@wrap_errors(wrapped_errors)
def encode(
	filename,
	count=False, fmt='base64', out=False, width=64,
	amino=False, degen=False, table=1,
):
	"""Encode data into various formats."""
	# TODO Read and decode data iteratively
	amino = auto_select_amino(filename, amino)
	data = read(filename, amino)
	decoded = dna.decode_iter(data, amino, table, degen)

	encoded = dna.encode(decoded, fmt)

	if count:
		counts = Counter(encoded)
		return (f"{count}\t{codon}" for codon, count in counts.most_common())

	lines = (str_join(line) for line in chunked(encoded, width))
	if out:
		file_ext = file_extension_for(fmt)
		outfile = f'{filename}.{file_ext}'
		with Path(outfile).open("w", encoding="UTF-8") as file:
			file.write(str_join(lines, "\n"))
			file.write("\n")  # Newline at the end
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
	seqs = list(read_sequences(filename, amino))

	if not width:
		# Set width from the total data size divided by sequence count
		file_size = os.path.getsize(filename)
		item_size = 1 if amino else 3
		size = file_size / item_size
		width = autowidth(size / len(mode), aspect=phi-1, base=64)
		echo(f'Automatically set image width: {width} px')

	rgb = np.vstack([
		dna_image_seq(
			seq, width, mode, fill=0,
			amino=amino, degen=degen, table=table)
		for seq in seqs
	])

	img = Image.fromarray(rgb, mode=mode)

	if amino and table != 1:
		outfile = filename + f".w{width}.{dna.BASE_ORDER}.t{table}.png"
	else:
		outfile = filename + f".w{width}.{dna.BASE_ORDER}.png"

	img.show(title=outfile)

	if out:
		img.save(outfile)


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
	decoded = dna.decode_iter(data, amino, table, degen)
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
	decoded = dna.decode_iter(data, amino, table, degen)

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
	decoded = dna.decode_iter(data, amino, table, degen)
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


def analyse_repeats(decoded, length=4, limit=2, fmt='codon'):
	"""Analyse codon sequence repeats."""
	[index, count] = count_sorted(dna.codon_sequences(decoded, length))
	repeats = index[count >= limit]

	echo("Repeated codon sequences:")
	# TODO Move conversion to a helper function
	b64 = (num_to_digits(num, base=64) for num in repeats)
	encoded = (str_join(dna.encode(codon, fmt)) for codon in b64)

	width = 32 if fmt in ['b', 'base64'] else 8
	yield from format_lines(encoded, width)


@arg('--amino',  '-a', help='Amino acid input')
@arg('--table',  '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen',  '-g', help='Degenerate data')
@arg('--fmt',    '-f', help='Output format', choices=['b', 'base64', 'c', 'codon'])
@arg('--limit',  '-l', help='Limit to at least this many repeats', type=int)
@arg('--seq',    '-s', help='Sequence length', type=int)
@wrap_errors(wrapped_errors)
def repeats(
	filename,
	seq=2, limit=2, fmt='codon',
	amino=False, degen=False, table=1,
):
	"""Find repeated codon sequences."""
	amino = auto_select_amino(filename, amino)
	data = read(filename, amino)
	decoded = dna.decode_iter(data, amino, table, degen)

	yield from analyse_repeats(decoded, length=seq, limit=limit, fmt=fmt)


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

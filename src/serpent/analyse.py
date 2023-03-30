#!/usr/bin/env python
# ruff: noqa: PLR0913
"""Serpent analysis."""

from __future__ import annotations

import fileinput
import itertools as itr
import os
from collections import Counter
from pathlib import Path

import argh
import matplotlib.pyplot as plt
import more_itertools as mit
import numpy as np
from argh.decorators import aliases, arg, wrap_errors
from PIL import Image

from serpent import ansi, dna
from serpent.bitmap import num_to_pixel
from serpent.block_elements import pixels_to_blocks
from serpent.convert.amino import aa_tables
from serpent.convert.digits import num_to_digits
from serpent.dsp import fft_spectra
from serpent.encoding import BASE64
from serpent.fasta import (
	ParseError,
	auto_select_amino,
	data_and_descriptions,
	find_fasta_files,
	find_fasta_sequences,
	read,
	read_sequences,
)
from serpent.fun import sort_values, str_join
from serpent.io import (
	check_paths,
	echo,
	file_extension_for,
	openhook,
	wait_user,
)
from serpent.mathematics import autowidth, phi
from serpent.palette import set_palette, spectrum_palette
from serpent.printing import format_counts, format_decoded, format_lines
from serpent.settings import (
	BASE_ORDER,
	COUNT_LIMIT,
	DEFAULT_COLOR,
	FLOW_DESCRIPTION_COLOR,
)
from serpent.stats import ac_peaks, autocorrelogram, count_sorted
from serpent.visual import (
	bin_choices,
	dna_image_seq,
	interactive,
	plot_histogram_sized,
	plot_sequence_counts,
)
from serpent.zigzag import zigzag_blocks, zigzag_text

fmt_choices = ['a', 'amino', 'b', 'base64', 'c', 'codon']
wrapped_errors = [AssertionError, ParseError]


@arg('--amino',  '-a', help='Amino acid input')
@arg('--table',  '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen',  '-g', help='Degenerate data')
@arg('--limit',  '-l', help='Peak limit', type=float)
@arg('--linear', '-n', help='Linear output (do not sort results)')
@arg('--plot',   '-p', help='Plot data')
@arg('--seq',    '-s', help='Sequence length', type=int)
@arg('--width',  '-w', help='Autocorrelogram width', type=int)
@wrap_errors(wrapped_errors)
def ac(
	filename,
	limit=0.05, linear=False, width=256, seq=1,
	plot=False,
	amino=False, degen=False, table=1,
):
	"""Plot autocorrelation of sequences."""
	amino = auto_select_amino(filename, amino)
	data = read(filename, amino)
	decoded = dna.decode(data, amino, table, degen)

	if seq > 1:
		decoded = dna.codon_sequences(decoded, seq)

	ac = autocorrelogram(decoded, width)

	if plot:
		plt.plot(ac, color=DEFAULT_COLOR)
		interactive()
		wait_user()
	else:
		peaks = ac_peaks(ac, limit)
		peaks.pop(0)  # First item is always 1
		peaks = peaks.items() if linear else sort_values(peaks, reverse=True)

		yield from (f"{peak}	{value:1.3f}" for peak, value in peaks)


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
			yield from format_lines(codons, width)


def cat(*inputs):
	"""Concatenate and print FASTA sequences from files."""
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

	yield from format_decoded(decoded)


@arg('--amino', '-a', help='Amino acid input')
@arg('--table', '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen', '-g', help='Degenerate data')
@arg('--count', '-c', help='Print counts')
@arg('--fmt',   '-f', help='Output format', choices=fmt_choices)
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

	if fmt in ['a', 'amino']:
		encoded = dna.to_amino(data, amino, table, degen)
	else:
		decoded = dna.decode_iter(data, amino, table, degen)
		encoded = dna.encode(decoded, fmt)

	if count:
		counts = Counter(encoded)
		return (f"{count}\t{codon}" for codon, count in counts.most_common())

	lines = (str_join(line) for line in mit.chunked(encoded, width))
	if out:
		file_ext = file_extension_for(fmt)
		outfile = f'{filename}.{file_ext}'
		with Path(outfile).open("w", encoding="UTF-8") as file:
			file.write(str_join(lines, "\n"))
			file.write("\n")  # Newline at the end
			print(f"Wrote: {outfile}")
	else:
		yield from lines


@arg('--amino', '-a', help='Amino acid input')
@arg('--table', '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen', '-g', help='Degenerate data')
@arg('--desc',  '-c', help='Output descriptions')
@arg('--width', '-w', help='Line width', type=int)
@arg('--mode',  '-m', help='Image mode', choices=('RGB', 'L'))
@wrap_errors(wrapped_errors)
def flow(
	*inputs,
	desc=False, width=64, mode='RGB',
	amino=False, degen=False, table=1,
):
	"""Encode data into Unicode block graphics."""
	paths = check_paths(inputs)
	amino_opt = amino

	for filename in paths:
		if len(inputs) > 1:
			print('file:', filename)

		amino = amino_opt
		amino = auto_select_amino(filename, amino)
		seqs = read_sequences(filename, amino)

		for seq in seqs:
			[decoded, description] = dna.decode_seq(seq, amino, table, degen)

			if desc:
				# TODO Add ansi.rgb_text and fix argument handling (needs back currently)
				colour = ansi.rgb(front=FLOW_DESCRIPTION_COLOR, back=(0, 0, 0))
				yield colour + description + ansi.RESET

			pixels = num_to_pixel(decoded, amino, degen)
			yield from pixels_to_blocks(pixels, width, mode=mode)


@arg('--amino', '-a', help='Amino acid input')
@arg('--table', '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen', '-g', help='Degenerate data')
@arg('--mode',  '-m', help='Image mode', choices=('RGB', 'L', 'P'))
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
	seqs = read_sequences(filename, amino)

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

	if mode == 'P':
		num_colours = 22 if amino else 64
		set_palette(img, spectrum_palette(num_colours, 0.75))

	if amino and table != 1:
		outfile = filename + f".w{width}.{BASE_ORDER}.t{table}.png"
	else:
		outfile = filename + f".w{width}.{BASE_ORDER}.png"

	img.show(title=outfile)

	if out:
		img.save(outfile)


@arg('--amino',  '-a', help='Amino acid input')
@arg('--table',  '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen',  '-g', help='Degenerate data')
@arg('--length', '-l', help='Fourier transform length', type=int)
@arg('--seq',    '-s', help='Sequence length', type=int)
@wrap_errors(wrapped_errors)
def fft(
	filename,
	length=None, seq=1,
	amino=False, degen=False, table=1,
):
	"""Plot Fourier transform of the DNA data."""
	amino = auto_select_amino(filename, amino)
	data = read(filename, amino)
	decoded = dna.decode(data, amino, table, degen)

	if seq > 1:
		decoded = np.fromiter(dna.codon_sequences(decoded, n=seq), dtype=np.uint64)

	# TODO Make a spectrogram by windowing
	[freqs, spectra] = fft_spectra(decoded, n=length)

	plt.plot(freqs, spectra, color=DEFAULT_COLOR)
	interactive()
	wait_user()


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
	seqs = np.fromiter(dna.codon_sequences(decoded, length), dtype=np.uint64)

	plot_histogram_sized(
		seqs,
		size=num or bins,
		multi=max(16, 2 ** length),  # cap to 16 * base = 1024
		color=DEFAULT_COLOR,
		cumulative=cumulative,
		density=density,
	)
	interactive()
	wait_user()


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

	plot_sequence_counts(decoded, n=length, color=DEFAULT_COLOR)
	interactive()
	wait_user()


@arg('--amino',   '-a', help='Amino acid input')
@arg('--table',   '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen',   '-g', help='Degenerate data')
@arg('--fmt',     '-f', help='Output format', choices=['b', 'base64', 'c', 'codon'])
@arg('--count',   '-c', help='Print counts')
@arg('--length',  '-l', help='Peptide length', type=int)
@arg('--missing', '-m', help='Missing peptides')
@wrap_errors(wrapped_errors)
def pep(
	filename,
	length=2, count=False, missing=False, fmt='base64',
	amino=False, table=1, degen=False,
):
	"""Peptide statistics."""
	amino = auto_select_amino(filename, amino)
	data = read(filename, amino)
	decoded = dna.decode_iter(data, amino, table, degen)
	dtype = f'U{length}'

	# TODO Allow using amino acid codes?
	encoded = dna.encode(decoded, fmt=fmt)
	peptides = (str_join(pep) for pep in mit.chunked(encoded, length))

	if count:
		counts = Counter(peptides)

		# TODO There must be a better way to use these iterators!
		groups = itr.groupby(counts.most_common(), lambda x: x[1])
		grouped = [
			(c, [k for k, c in group])
			for c, group in groups
			if c > 1
		]
		for count, values in grouped:
			yield f"-- {count} times --"
			yield from format_lines(sorted(values), 32)
		return None

	if not missing:
		yield from format_lines(peptides, 32)
	else:
		print("Peptides not appearing:\n")
		combos = np.fromiter(map(str_join, itr.product(BASE64, repeat=length)), dtype=dtype)
		absent = combos[[combo not in peptides for combo in combos]]

		yield from format_lines(absent, 32)


def analyse_repeats(decoded, length=4, limit=2, fmt='codon'):
	"""Analyse codon sequence repeats."""
	[index, count] = count_sorted(dna.codon_sequences(decoded, length))
	repeats = index[count >= limit]

	echo("Repeated codon sequences:")
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
		flow,
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

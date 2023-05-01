#!/usr/bin/env python
# ruff: noqa: PLR0913
"""Serpent analysis."""

from __future__ import annotations

import fileinput
import itertools as itr
import os
import sys
from collections import Counter, OrderedDict

import argh
import matplotlib.pyplot as plt
import more_itertools as mit
import numpy as np
from argh.decorators import aliases, arg, wrap_errors
from PIL import Image

from serpent import ansi, dna
from serpent.cli.encode import encode_data
from serpent.cli.flow import flow_blocks, verbose_flow_blocks
from serpent.cli.image import dna_image
from serpent.cli.repeats import analyse_repeats
from serpent.cli.zigzag import zigzag_blocks, zigzag_text
from serpent.convert.amino import aa_tables, aminos_for_table
from serpent.convert.degenerate import is_degenerate
from serpent.convert.split import split_aminos, split_encoded
from serpent.dsp import fft_spectra
from serpent.encoding import BASE64
from serpent.fasta import (
	ParseError,
	auto_select_amino,
	descriptions_and_data,
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
	image_name_for,
	openhook,
	wait_user,
	write_iterable,
)
from serpent.mathematics import autowidth_for
from serpent.palette import rgb_to_hex, spectrum_layer_colours
from serpent.printing import (
	format_counter,
	format_decoded,
	format_lines,
	format_quasar,
	format_quasar_pulses,
	format_split,
	reflow,
)
from serpent.settings import (
	COUNT_LIMIT,
	DEFAULT_COLOR,
)
from serpent.spatial import amino_path_3d
from serpent.stats import ac_peaks, autocorrelogram, count_sorted, quasar_pulses
from serpent.visual import (
	bin_choices,
	dna_quasar_seq,
	interactive,
	plot_amino_labels,
	plot_directions,
	plot_histogram_sized,
	plot_sequence_counts,
)

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


@arg('--degen-only', '-G', help='Only degenerate')
@arg('--stats',  '-s', help='Show statistics')
@arg('--limit',  '-l', help='Limit codon stats to N most common')
@arg('--width',  '-w', help='Codons per line')
def codons(filename, degen_only=False, width=20, stats=False, limit=COUNT_LIMIT):
	"""Print codons and statistics."""
	seqs = read_sequences(filename)

	for seq in seqs:
		[descriptions, data] = descriptions_and_data(seq)

		yield from descriptions
		codons = dna.get_codons_iter(data)
		if degen_only:
			# TODO Print descriptions only when there is data by using mit.peekable
			codons = filter(is_degenerate, codons)

		if stats:
			counts = Counter(codons)
			print("Counts:")
			yield from format_counter(counts, limit=limit)

			unique = counts.keys()  # np.unique(codons)
			yield from format_lines(unique, width)

			total = np.sum([*counts.values()])
			print(f"Unique codons: {len(counts)} (total: {total})")
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
	paths = check_paths(inputs, recurse=True)
	with fileinput.input(paths, mode='r', openhook=openhook) as fi:
		if seq:
			yield from find_fasta_sequences(fi, num)
		else:
			yield from find_fasta_files(fi)


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

	yield from format_decoded(decoded, degen)


@arg('--amino', '-a', help='Amino acid input')
@arg('--table', '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen', '-g', help='Degenerate data')
@arg('--fmt',   '-f', help='Output format', choices=fmt_choices)
@arg('--out',   '-o', help='Write out to file')
@arg('--width', '-w', help='Line width', type=int)
@wrap_errors(wrapped_errors)
def encode(
	*inputs,
	fmt='codon', out=False, width=64,
	amino=False, degen=False, table=1,
):
	"""Encode data into various formats."""
	amino_opt = amino

	for filename in check_paths(inputs):
		if len(inputs) > 1:
			yield f';file:{filename}'
		amino = auto_select_amino(filename, amino_opt)
		seqs = read_sequences(filename, amino)

		for seq in seqs:
			[descriptions, data] = descriptions_and_data(seq)
			yield from descriptions

			encoded = encode_data(data, fmt, amino, table, degen)
			lines = reflow(encoded, width)

			if out:
				ext = file_extension_for(fmt)
				outfile = f'{filename}.enc.{ext}'
				write_iterable(lines, outfile)
			else:
				yield from lines


@arg('--amino', '-a', help='Amino acid input')
@arg('--table', '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen', '-g', help='Degenerate data')
@arg('--fmt',   '-f', help='Output format', choices=fmt_choices)
@arg('--summary', '-s', help='Summarise counts per file')
@wrap_errors(wrapped_errors)
def count(
	*inputs,
	fmt=None, summary=False,
	amino=False, degen=False, table=1,
):
	"""Count data."""
	amino_opt = amino

	for filename in check_paths(inputs):
		if len(inputs) > 1:
			yield f';file:{filename}'
		amino = auto_select_amino(filename, amino_opt)
		show_gc = fmt is None and not amino

		seqs = read_sequences(filename, amino)
		total = Counter()

		for seq in seqs:
			[descriptions, data] = descriptions_and_data(seq)

			if fmt is not None:
				data = encode_data(data, fmt, amino, table, degen)

			counts = Counter(data)
			total.update(counts)

			if not summary:
				yield from descriptions
				yield from format_counter(counts, show_gc)

		if summary:
			print(f'Summary: {filename}')
			yield from format_counter(total, show_gc)


@arg('--amino', '-a', help='Amino acid input')
@arg('--table', '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen', '-g', help='Degenerate data')
@arg('--fmt',   '-f', help='Output format', choices=['a', 'amino', 'c', 'codon'])
@arg('--words', '-r', help='Split by either stop or start codons')
@arg('--width', '-w', help='Line width', type=int)
@wrap_errors(wrapped_errors)
def split(
	*inputs,
	fmt='codon', words=False, width=64,
	amino=False, degen=False, table=1,
):
	"""Split data in various ways."""
	amino_opt = amino
	split = 'r' if words else 'n'

	for filename in check_paths(inputs):
		if len(inputs) > 1:
			yield f';file:{filename}'
		amino = auto_select_amino(filename, amino_opt)
		seqs = read_sequences(filename, amino)

		for seq in seqs:
			[descriptions, data] = descriptions_and_data(seq)
			yield from descriptions

			encoded = encode_data(data, fmt, amino, table, degen)
			regions = split_encoded(encoded, fmt, split=split)
			yield from format_split(regions, width, split=split)


@arg('--amino',   '-a', help='Amino acid input')
@arg('--table',   '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen',   '-g', help='Degenerate data')
@arg('--desc',    '-c', help='Output descriptions')
@arg('--fmt',     '-f', help='Output format', choices=fmt_choices)
@arg('--mode',    '-m', help='Image mode', choices=('RGB', 'L', 'P'))
@arg('--verbose', '-v', help='Verbose mode')
@arg('--width',   '-w', help='Line width', type=int)
@wrap_errors(wrapped_errors)
def flow(
	*inputs,
	width=64, mode='RGB', fmt=None, desc=False, verbose=False,
	amino=False, degen=False, table=1,
):
	"""Encode data into Unicode block graphics."""
	amino_opt = amino

	for filename in check_paths(inputs):
		if len(inputs) > 1:
			yield f';file:{filename}'
		amino = auto_select_amino(filename, amino_opt)
		seqs = read_sequences(filename, amino)

		for seq in seqs:
			[descriptions, data] = descriptions_and_data(seq)

			if verbose:
				yield from descriptions
				yield from verbose_flow_blocks(
					data, width, mode, fmt,
					amino=amino, degen=degen, table=table)
			else:
				if desc:
					yield from (ansi.dim_text(d) for d in descriptions)
				yield from flow_blocks(
					data, width, mode,
					amino=amino, degen=degen, table=table)


@arg('--amino',  '-a', help='Amino acid input')
@arg('--table',  '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen',  '-g', help='Degenerate data')
@arg('--mode',   '-m', help='Image mode', choices=('RGB', 'L', 'P', 'Q'))
@arg('--out',    '-o', help='Write image to file')
@arg('--length', '-l', help='Peptide length (for Q mode)', type=int)
@arg('--width',  '-w', help='Image width', type=int)
@wrap_errors(wrapped_errors)
def image(
	filename,
	length=1, width=None, mode="RGB", out=False,
	amino=False, degen=False, table=1,
):
	"""Visualise FASTA data as images."""
	amino = auto_select_amino(filename, amino)
	seqs = read_sequences(filename, amino)

	if not width:
		channels = 3 if mode == 'Q' else len(mode)
		file_size = os.path.getsize(filename)
		width = autowidth_for(file_size, amino, channels)
		echo(f'Automatically set image width: {width} px')

	img = dna_image(
		seqs, width, mode,
		amino=amino, degen=degen, table=table, length=length)

	outfile = image_name_for(
		filename, width, mode,
		amino=amino, degen=degen, table=table, length=length)

	img.show(title=outfile)
	if out:
		img.save(outfile)


@arg('--amino', '-a', help='Amino acid input')
@arg('--table', '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen', '-g', help='Degenerate data')
@arg('--plot',  '-p', help='Plot data')
@arg('--count', '-c', help='Print counts')
@arg('--cumulative', '-m', help='Cumulative')
@wrap_errors(wrapped_errors)
def pulse(
	filename,
	count=False, cumulative=False, plot=False,
	amino=False, degen=False, table=1,
):
	"""Pulse repetition intervals for symbol repeat lengths."""
	amino = auto_select_amino(filename, amino)
	seqs = read_sequences(filename, amino)
	key = aminos_for_table(table)

	if plot:
		colour_map = spectrum_layer_colours(amino)
		colours = OrderedDict([
			[aa, rgb_to_hex(rgb)]
			for aa, rgb in zip(key, colour_map)
		])
		size = 10
		title = 'SRI counts' if count else 'Symbol repetition intervals'

		ax = plt.axes()
		ax.set_title(f'{title} of\n{filename}')

		maxheight = maxscale = 0

		for seq in seqs:
			[aminos, descriptions] = dna.decode_seq(seq, amino, table, degen, dna.to_amino)
			pulses, height, scale = quasar_pulses(aminos, cumulative=cumulative, key=key)
			maxheight = max(height, maxheight)
			maxscale = max(scale, maxscale)

			for i, item in enumerate(reversed(pulses.items())):
				[aa, pulse] = item

				colour = colours[aa]
				mx = len(pulses) - 1
				pos = mx - i

				if count:
					base = 10
					data = count_sorted(pulse)
					ax.plot(*data, colour)

					ax.set_xscale('log', base=base)
					ax.set_yscale('log', base=base)

					ax.set_xlabel('repeat length')
					ax.set_ylabel('count')
				else:
					data = pulse
					ax.plot(data + pos * 100, colour)

					# ax.set_xscale('log', base=10)
					# ax.set_yscale('log', base=10)

					ax.text(
						-15, pos * 100, aa, color=colour,
						ha='right', size=size, fontweight='semibold'
					)
					ax.set_xlabel('index of repetition')
					ax.set_ylabel('interval length (space between symbols)')

		if count:
			mx = len(key) - 1
			for i, aa in enumerate(reversed(key)):
				colour = colours[aa]
				ax.text(
					1 + maxscale * (i / mx),
					2 + maxheight * (1 - i / mx),
					aa, color=colour,
					size=size, fontweight='semibold'
				)

		interactive()
		wait_user()
	else:
		for seq in seqs:
			[aminos, descriptions] = dna.decode_seq(seq, amino, table, degen, dna.to_amino)
			pulses, height, scale = quasar_pulses(aminos, cumulative=cumulative, key=key)

		yield from descriptions
		yield from format_quasar(pulses.keys())  # Print symbols
		yield from format_quasar_pulses(pulses, height)
		yield f'scale: {scale}'


@arg('--amino', '-a', help='Amino acid input')
@arg('--table', '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen', '-g', help='Degenerate data')
@arg('--mod',   '-m', help='Modulo of pattern data', type=int)
@arg('--log',   '-l', help='Logarithmic result')
@arg('--cumulative', '-c', help='Cumulative')
@arg('--test',  help='Show test images (replace data with linear numbers)')
@wrap_errors(wrapped_errors)
def quasar(
	filename,
	cumulative=False, log=False, mod=0, test=False,
	amino=False, degen=False, table=1,
):
	"""Visualise symbol repeats as image."""
	# TODO Show data as block flow
	# TODO Show animated image of multiple sequences?
	amino = auto_select_amino(filename, amino)
	seqs = read_sequences(filename, amino)
	key = aminos_for_table(table)

	rgb = np.vstack([
		dna_quasar_seq(seq, amino, degen, table, cumulative, log, mod, test, key=key)
		for seq in seqs
	])

	img = Image.fromarray(rgb, mode='L')
	img.show()


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


@arg('--amino',  '-a', help='Amino acid input')
@arg('--table',  '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen',  '-g', help='Degenerate data')
@arg('--split',  '-s', help='Split by stop aNd/oR start codons', choices=('n', 'r'))
@aliases('vec')
def vectors(filename, split='', amino=False, table=1, degen=False):
	"""Visualise amino acids in 3D vector space."""
	amino = auto_select_amino(filename, amino)
	seqs = read_sequences(filename, amino)

	scale = np.ones(3)
	title = f'{filename} (peptides)' if split else filename
	ax = plt.axes(projection='3d')
	ax.set_title(title)

	for index, seq in enumerate(seqs):
		[aminos, _] = dna.decode_seq(seq, amino, table, degen, dna.to_amino)

		if split:
			for peptide in split_aminos(aminos, split=split):
				# TODO Does it occur often that start and stop codons are next
				# to each other, and what does it mean? (Stop not interpreted
				# as stop, but amino acid?)
				if len(peptide) == 0:
					continue
				dirs = amino_path_3d(peptide)
				maximum = np.amax(np.abs(dirs), axis=0)
				scale = np.max([scale, maximum], axis=0)
				plot_directions(ax, dirs)
		else:
			# TODO Map descriptions and regions to different colours from text?
			color = DEFAULT_COLOR if index == 0 and not split else None
			dirs = amino_path_3d(aminos)
			maximum = np.amax(np.abs(dirs), axis=0)
			scale = np.max([scale, maximum], axis=0)
			plot_directions(ax, dirs, color=color)

	plot_amino_labels(ax, scale)
	plt.axis('equal')
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
		count,
		decode,
		encode,
		fft,
		find,
		flow,
		hist,
		image,
		pep,
		pulse,
		quasar,
		repeats,
		seq,
		split,
		vectors,
		zigzag,
	])
	# parser.set_default_command(serpent)
	try:
		return parser.dispatch()
	except (BrokenPipeError, OSError):
		sys.exit(32)


if __name__ == '__main__':
	decoded = main()

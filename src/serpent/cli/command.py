#!/usr/bin/env python
# ruff: noqa: PLR0913
"""Serpent analysis."""

from __future__ import annotations

import fileinput
import itertools as itr
import sys
from collections import Counter
from pathlib import Path

import argh
import matplotlib.pyplot as plt
import numpy as np
from argh.decorators import aliases, arg, wrap_errors
from PIL import Image

from serpent import dna
from serpent.cli.encode import encode_data, encode_sequences
from serpent.cli.entropy import format_entropy, plot_entropy
from serpent.cli.flow import flow_blocks, verbose_flow_blocks
from serpent.cli.image import dna_image
from serpent.cli.pep import all_symbols_for, peptides_of_length
from serpent.cli.pulse import (
	pulse_plot_sequences,
	pulse_text,
)
from serpent.cli.quasar import dna_quasar_seq
from serpent.cli.tide import plot_tides, tide_sequence, tide_total
from serpent.cli.walk import walk_sequence
from serpent.cli.zigzag import zigzag_blocks, zigzag_text
from serpent.convert.amino import aa_tables, aminos_for_table
from serpent.convert.dnt import is_degenerate
from serpent.convert.split import split_aminos, split_encoded
from serpent.fun import second, sort_values, str_join
from serpent.io.fasta import (
	ParseError,
	auto_select_amino,
	descriptions_and_data,
	find_fasta_files,
	find_fasta_sequences,
	is_fasta,
	read,
	read_sequences,
	regex_no_match,
)
from serpent.io.files import (
	check_paths,
	file_extension_for,
	image_name_for,
	openhook,
	readlines,
	write_iterable,
)
from serpent.io.printing import (
	auto_line_width_for,
	echo,
	format_counter,
	format_decoded,
	format_lines,
	format_split,
	info,
	wait_user,
)
from serpent.math.basic import autowidth_for
from serpent.math.dsp import fft_spectra
from serpent.math.statistic import ac_peaks, autocorrelogram, ewma, quasar_pulses
from serpent.settings import (
	BASE_ORDER,
	COUNT_LIMIT,
	DEFAULT_COLOR,
)
from serpent.spatial.path import amino_path_3d
from serpent.visual import ansi
from serpent.visual.plot import (
	bin_choices,
	interactive,
	plot_amino_labels,
	plot_directions,
	plot_histogram_sized,
	plot_sequence_counts,
)

fmt_choices = ['a', 'amino', 'c', 'codon']
wrapped_errors = [AssertionError, ParseError]


@arg('--amino',  '-a', help='Amino acid input')
@arg('--table',  '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen',  '-g', help='Degenerate data')
@arg('--limit',  '-l', help='Peak limit', type=float)
@arg('--linear', '-n', help='Linear output (do not sort results)')
@arg('--plot',   '-p', help='Plot data')
@arg('--seql',   '-q', help='Sequence length', type=int)
@arg('--width',  '-w', help='Autocorrelogram width', type=int)
@wrap_errors(wrapped_errors)
def ac(
	filename,
	limit=0.05, linear=False, width=256, seql=1,
	plot=False,
	amino=False, degen=False, table=1,
):
	"""Autocorrelation of sequences (with plot option)."""
	amino = auto_select_amino(filename, amino)
	data = read(filename, amino)
	decoded = dna.decode(data, amino, table, degen)

	if seql > 1:
		# FIXME Change autocorrelogram to accept iterators
		decoded = [*dna.codon_sequences(decoded, seql)]

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

	for sequence in seqs:
		[descriptions, data] = descriptions_and_data(sequence)

		yield from descriptions
		codons = dna.get_codons_iter(data)
		if degen_only:
			# TODO Print descriptions only when there is data by using mit.peekable
			codons = filter(is_degenerate, codons)

		if stats:
			counts = Counter(codons)
			yield from format_counter(counts, top=limit)

			unique = counts.keys()  # np.unique(codons)
			yield f"unique\t{len(unique)}"
			yield from format_lines(unique, width)
		else:
			yield from format_lines(codons, width)


@arg('--amino', '-a', help='Amino acid input')
@arg('--base',  '-b', help='Base information unit', type=float)
@arg('--plot',  '-p', help='Plot data')
@arg('--reg',   '-r', help='Filter sequences by regexp', type=str)
@arg('--seql',  '-q', help='Sequence length', type=int)
@arg('--step',  '-s', help='Step size', type=int)
@aliases('ent')
def entropy(*inputs, base=2.0, amino=False, plot=False, reg=None, seql=None, step=None):
	"""Entropy and other information theory statistics."""
	amino_opt = amino

	for filename in check_paths(inputs):
		if len(inputs) > 1:
			info(f'file: {filename}')

		if is_fasta(filename):
			amino = auto_select_amino(filename, amino_opt)
			seqs = read_sequences(filename, amino)

			for index, sequence in enumerate(seqs):
				[descriptions, data] = descriptions_and_data(sequence)
				if regex_no_match(reg, descriptions):
					continue

				if plot:
					color = DEFAULT_COLOR if index == 0 else None
					plot_entropy(data, base, seql, step, color=color)
				else:
					yield from descriptions
					yield from format_entropy(data, base, seql, step)
		else:
			lines = readlines(filename)
			counts = Counter()
			for line in lines:
				if regex_no_match(reg, line):
					continue
				counts.update(line)
			yield from format_entropy(counts, base, seql, step)

	if plot:
		interactive()
		wait_user()


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
@arg('--reg',   '-r', help='Filter sequences by regexp on descriptions', type=str)
@arg('--width', '-w', help='Line width', type=int)
@wrap_errors(wrapped_errors)
def encode(
	*inputs,
	fmt='codon', out=False, width=64, reg=None,
	amino=False, degen=False, table=1,
):
	"""Encode data into various formats."""
	amino_opt = amino

	for filename in check_paths(inputs):
		if len(inputs) > 1:
			info(f'file: {filename}')
		amino = auto_select_amino(filename, amino_opt)
		seqs = read_sequences(filename, amino)

		lines = encode_sequences(
			seqs, width, fmt, reg=reg,
			amino=amino, table=table, degen=degen)

		if out:
			ext = file_extension_for(fmt)
			outfile = f'{filename}.enc.{ext}'
			write_iterable(lines, outfile)
		else:
			yield from lines


@arg('--amino',   '-a', help='Amino acid input')
@arg('--table',   '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen',   '-g', help='Degenerate data')
@arg('--fmt',     '-f', help='Output format', choices=fmt_choices)
@arg('--seql',    '-q', help='Sequence length', type=int)
@arg('--limit',   '-l', help='Limit to at least this many repeats', type=int)
@arg('--summary', '-s', help='Summarise counts per file')
@wrap_errors(wrapped_errors)
def count(
	*inputs,
	fmt=None, summary=False, seql=1, limit=1,
	amino=False, degen=False, table=1,
):
	"""Count data."""
	amino_opt = amino

	for filename in check_paths(inputs):
		if len(inputs) > 1:
			info(f'file: {filename}')
		amino = auto_select_amino(filename, amino_opt)
		show_gc = fmt is None and not amino

		seqs = read_sequences(filename, amino)
		total = Counter()

		for sequence in seqs:
			[descriptions, data] = descriptions_and_data(sequence)

			if fmt is not None:
				data = encode_data(data, fmt, amino, table, degen)

			if seql > 1:
				data = peptides_of_length(data, seql)

			counts = Counter(data)
			total.update(counts)

			if not summary:
				yield from descriptions
				yield from format_counter(counts, show_gc, limit=limit)

		if summary:
			print(f'Summary: {filename}')
			yield from format_counter(total, show_gc, limit=limit)


@arg('--amino', '-a', help='Amino acid input')
@arg('--table', '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen', '-g', help='Degenerate data')
@arg('--fmt',   '-f', help='Output format', choices=['a', 'amino', 'c', 'codon'])
@arg('--split', '-s', help='Splitter type', choices=['f', 'n', 'r'])
@arg('--width', '-w', help='Line width', type=int)
@wrap_errors(wrapped_errors)
def split(
	*inputs,
	fmt='codon', split='f', width=64,
	amino=False, degen=False, table=1,
):
	"""Split data in various ways."""
	amino_opt = amino

	for filename in check_paths(inputs):
		if len(inputs) > 1:
			info(f'file: {filename}')
		amino = auto_select_amino(filename, amino_opt)
		seqs = read_sequences(filename, amino)

		for sequence in seqs:
			[descriptions, data] = descriptions_and_data(sequence)
			yield from descriptions

			encoded = encode_data(data, fmt, amino, table, degen)
			regions = split_encoded(encoded, fmt, table, split=split)
			yield from format_split(regions, width, split=split)


@arg('--amino',   '-a', help='Amino acid input')
@arg('--table',   '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen',   '-g', help='Degenerate data')
@arg('--desc',    '-c', help='Output descriptions')
@arg('--fmt',     '-f', help='Output format', choices=fmt_choices)
@arg('--mode',    '-m', help='Image mode', choices=('RGB', 'L', 'P'))
@arg('--reg',     '-r', help='Filter sequences by regexp on descriptions', type=str)
@arg('--verbose', '-v', help='Verbose mode')
@arg('--width',   '-w', help='Line width', type=int)
@wrap_errors(wrapped_errors)
def flow(
	*inputs,
	width=64, mode='RGB', fmt=None, desc=False, reg=None, verbose=False,
	amino=False, degen=False, table=1,
):
	"""Visualise and explore data with Unicode block graphics."""
	amino_opt = amino

	for filename in check_paths(inputs):
		if len(inputs) > 1:
			yield f';file:{filename}'
		amino = auto_select_amino(filename, amino_opt)
		seqs = read_sequences(filename, amino)

		for sequence in seqs:
			[descriptions, data] = descriptions_and_data(sequence)
			if regex_no_match(reg, descriptions):
				continue

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
@arg('--seql',   '-q', help='Sequence length (for Q mode)', type=int)
@arg('--width',  '-w', help='Image width', type=int)
@wrap_errors(wrapped_errors)
def image(
	filename,
	seql=1, width=None, mode="RGB", out=False,
	amino=False, degen=False, table=1,
):
	"""Visualise FASTA data as images."""
	amino = auto_select_amino(filename, amino)
	seqs = read_sequences(filename, amino)

	if not width:
		channels = 3 if mode == 'Q' else len(mode)
		file_size = Path(filename).stat().st_size
		width = autowidth_for(file_size, amino, channels)
		echo(f'Automatically set image width: {width} px')

	img = dna_image(
		seqs, width, mode,
		amino=amino, degen=degen, table=table, length=seql)

	outfile = image_name_for(
		filename, width, mode,
		amino=amino, degen=degen, table=table, length=seql)

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
	"""Pulse repetition intervals for symbol repeat lengths (with plot option)."""
	amino = auto_select_amino(filename, amino)
	seqs = read_sequences(filename, amino)
	key = aminos_for_table(table)

	if plot:
		title = 'SRI counts' if count else 'Symbol repetition intervals'
		ax = plt.axes()
		ax.set_title(f'{title} of\n{filename}')

		pulse_plot_sequences(
			ax, seqs, key,
			amino=amino, table=table, degen=degen,
			count=count, cumulative=cumulative,
		)
		interactive()
		wait_user()
	else:
		for sequence in seqs:
			[aminos, descriptions] = dna.decode_seq(
				sequence,
				amino, table, degen,
				dna.to_amino
			)
			pulses, height, scale = quasar_pulses(
				aminos,
				cumulative=cumulative,
				key=key
			)

			yield from descriptions
			yield from pulse_text(pulses, height, scale)


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
	"""Visualise symbol repeats as images."""
	# TODO Show data as block flow
	# TODO Show animated image of multiple sequences?
	amino = auto_select_amino(filename, amino)
	seqs = read_sequences(filename, amino)
	key = aminos_for_table(table)

	rgb = np.vstack([
		dna_quasar_seq(
			sequence,
			cumulative=cumulative,
			log=log,
			mod=mod,
			test=test,
			key=key,
			amino=amino, table=table, degen=degen)
		for sequence in seqs
	])

	img = Image.fromarray(rgb, mode='L')
	img.show()


@arg('--amino',  '-a', help='Amino acid input')
@arg('--table',  '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen',  '-g', help='Degenerate data')
@arg('--length', '-l', help='Fourier transform length', type=int)
@arg('--seql',   '-q', help='Sequence length', type=int)
@wrap_errors(wrapped_errors)
def fft(
	filename,
	length=None, seql=1,
	amino=False, degen=False, table=1,
):
	"""Plot Fourier transform of the DNA data."""
	amino = auto_select_amino(filename, amino)
	data = read(filename, amino)
	decoded = dna.decode(data, amino, table, degen)

	if seql > 1:
		decoded = np.fromiter(dna.codon_sequences(decoded, n=seql), dtype=np.uint64)

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
@arg('--seql',  '-q', help='Sequence length', type=int)
@wrap_errors(wrapped_errors)
def hist(
	filename,
	bins='base',
	num=None,
	seql=1,
	density=False,
	cumulative=False,
	amino=False, degen=False, table=1,
):
	"""Plot DNA data histograms."""
	amino = auto_select_amino(filename, amino)
	data = read(filename, amino)
	decoded = dna.decode_iter(data, amino, table, degen)
	seqs = np.fromiter(dna.codon_sequences(decoded, n=seql), dtype=np.uint64)

	plot_histogram_sized(
		seqs,
		size=num or bins,
		multi=max(16, 2 ** seql),  # cap to 16 * base = 1024
		color=DEFAULT_COLOR,
		cumulative=cumulative,
		density=density,
	)
	interactive()
	wait_user()


@arg('--amino',  '-a', help='Amino acid input')
@arg('--table',  '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen',  '-g', help='Degenerate data')
@arg('--seql',   '-q', help='Sequence length', type=int)
@wrap_errors(wrapped_errors)
def seq(filename, seql=1, amino=False, degen=False, table=1):
	"""Plot DNA sequence count statistics."""
	amino = auto_select_amino(filename, amino)
	data = read(filename, amino)
	decoded = dna.decode_iter(data, amino, table, degen)

	plot_sequence_counts(decoded, n=seql, color=DEFAULT_COLOR)
	interactive()
	wait_user()


@arg('--cumulative', '-m', help='Cumulative probabilities')
@arg('--norm',   '-n', help='Normalise series to unity')
@arg('--slopes', '-l', help='Cumulative series')
@arg('--ref',    '-r', help='Reference sequence length', type=int)
@arg('--refstep', '-t', help='Reference step size (smoothing by overlap)', type=int)
@arg('--seql',   '-q', help='Sequence length', type=int)
@arg('--step',   '-s', help='Step size (smoothing by overlap)', type=int)
@wrap_errors(wrapped_errors)
def tide(
	*inputs, seql=1024, step=None, ref=64, refstep=1,
	cumulative=False, norm=False, slopes=False,
):
	"""Visualise averaged relative symbol frequencies."""
	symbols = BASE_ORDER  # + DEGENERATE

	for filename in check_paths(inputs):
		if len(inputs) > 1:
			info(f'file: {filename}')
		seqs = read_sequences(filename, amino=False)

		for sequence in seqs:
			[descriptions, data] = descriptions_and_data(sequence)
			yield from descriptions

			# Reference
			alpha = '33'
			tides = tide_sequence(
				data, symbols, ref, step=refstep,
				cumulative=cumulative,
				slopes=slopes,
			)
			yield f'ref:\t{tides.shape}'
			plot_tides(tides, symbols, alpha)
			tide_total(tides, norm=norm, color='#000000' + alpha)

			# Actual
			alpha = 'ff'
			tides = tide_sequence(
				data, symbols, seql, step=step,
				cumulative=cumulative,
				slopes=slopes,
			)
			yield f'tides:\t{tides.shape}'
			plot_tides(tides, symbols, alpha)
			tide_total(tides, norm=norm, color='#000000' + alpha)

			# EWMA
			alpha = '77'
			tides = tide_sequence(
				data, symbols, ref, step=1,
				cumulative=cumulative,
				slopes=slopes,
			)
			ewa = np.array([
				ewma(tide, alpha=0.1, window_size=ref)
				for tide in tides.T
			]).T
			yield f'ewa:\t{ewa.shape}'
			plot_tides(ewa, symbols, alpha)
			tide_total(ewa, norm=norm, color='#000000' + alpha)

	interactive()
	wait_user()


@arg('--amino',  '-a', help='Amino acid input')
@arg('--table',  '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen',  '-g', help='Degenerate data')
@arg('--split',  '-s', help='Split by stop aNd/oR start codons', choices=('n', 'r'))
@aliases('vec')
def vectors(filename, split='', amino=False, table=1, degen=False):
	"""Spatial visualisation of amino acids in 3D vector space."""
	amino = auto_select_amino(filename, amino)
	seqs = read_sequences(filename, amino)

	scale = np.ones(3)
	title = f'{filename} (peptides)' if split else filename
	ax = plt.axes(projection='3d')
	ax.set_title(title)

	for index, sequence in enumerate(seqs):
		[aminos, _] = dna.decode_seq(sequence, amino, table, degen, dna.to_amino)

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


# @arg('--amino',   '-a', help='Amino acid input')
# @arg('--table',   '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen',   '-g', help='Degenerate data')
@arg('--norm',    '-n', help='Normalise paths to unit square')
@arg('--length',  '-l', help='Exact sequence length as --seql/-q')
@arg('--seql',    '-q', help='Path segments (sequence length with --length)', type=int)
@arg('--step',    '-s', help='Step size (smoothing by overlap)', type=int)
@arg('--test',    '-t', help='Test path quality (compare to blueprint)')
@arg('--unit',    '-u', help='Stay on unit square')
@arg('--reg',     '-r', help='Filter sequences by regexp on descriptions', type=str)
@wrap_errors(wrapped_errors)
def walk(
	filename,
	reg=None,
	length=False, seql=1024, step=None, test=False, unit=False, norm=False,
	# amino=False, table=1,
	degen=False,
):
	"""Spatial walk visualisation of nucleotides in 2D space."""
	amino = auto_select_amino(filename, False)
	if amino:
		err_msg = 'Walk command only works for nucleotide data.'
		raise NotImplementedError(err_msg)
	seqs = read_sequences(filename, amino)

	for index, sequence in enumerate(seqs):
		[descriptions, data] = descriptions_and_data(sequence)
		if regex_no_match(reg, descriptions):
			continue
		yield from descriptions

		color = DEFAULT_COLOR if index == 0 else None
		seql = seql if length else max([1, len(data) // seql])

		if test:
			path = walk_sequence(
				data, seql=1, step=1,
				degen=degen, norm=norm, unit=unit
			)
			plt.plot(*path, color='#0ff7')

		path = walk_sequence(data, seql, step, degen=degen, norm=norm, unit=unit)
		plt.plot(*path, color=color)

	plt.axis('equal')
	interactive()
	wait_user()


@arg('--amino',   '-a', help='Amino acid input')
@arg('--table',   '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen',   '-g', help='Degenerate data')
@arg('--fmt',     '-f', help='Output format', choices=fmt_choices)
@arg('--seql',    '-q', help='Sequence length', type=int)
@arg('--missing', '-m', help='Missing peptides')
@arg('--unique',  '-u', help='Unique peptides')
@wrap_errors(wrapped_errors)
def pep(
	filename,
	fmt=None, seql=2, missing=False, unique=False,
	amino=False, table=1, degen=False,
):
	"""Peptides found or missing in the data."""
	amino = auto_select_amino(filename, amino)
	fmt = fmt or ('amino' if amino else 'codon')
	width = auto_line_width_for(fmt, seql, base=8, indent=8)

	data = read(filename, amino)
	encoded = encode_data(data, fmt, amino, table, degen)

	peptides = peptides_of_length(encoded, seql)

	if not missing:
		if unique:
			peptides = map(str_join, sorted(set(peptides)))

		yield from format_lines(peptides, width)
	else:
		print("Peptides not appearing:\n")
		combos = all_symbols_for(fmt, seql, table)
		absent = combos[[combo not in peptides for combo in combos]]

		yield from format_lines(absent, width)


@arg('--amino',   '-a', help='Amino acid input')
@arg('--table',   '-t', help='Amino acid translation table', choices=aa_tables)
@arg('--degen',   '-g', help='Degenerate data')
@arg('--fmt',     '-f', help='Output format', choices=fmt_choices)
@arg('--limit',   '-l', help='Limit to at least this many repeats', type=int)
@arg('--seql',    '-q', help='Sequence length', type=int)
@aliases('repeats')
@wrap_errors(wrapped_errors)
def pepcount(
	filename,
	fmt=None, limit=2, seql=2,
	amino=False, table=1, degen=False,
):
	"""Repeated peptide counts."""
	amino = auto_select_amino(filename, amino)
	fmt = fmt or ('amino' if amino else 'codon')
	width = auto_line_width_for(fmt, seql, base=8, indent=8)

	data = read(filename, amino)
	encoded = encode_data(data, fmt, amino, table, degen)

	peptides = peptides_of_length(encoded, seql)
	counts = Counter(peptides)

	groups = itr.groupby(counts.most_common(), second)
	grouped = (
		(count, [peptide for peptide, _ in group])
		for count, group in groups
		if count >= limit
	)
	for count, values in grouped:
		yield f"-- {count} times --"
		yield from format_lines(sorted(values), width)


def main():
	parser = argh.ArghParser()
	parser.add_commands([
		ac,
		cat,
		codons,
		count,
		decode,
		encode,
		entropy,
		fft,
		find,
		flow,
		hist,
		image,
		pep,
		pepcount,
		pulse,
		quasar,
		seq,
		split,
		tide,
		vectors,
		walk,
		zigzag,
	])
	# parser.set_default_command(serpent)
	try:
		return parser.dispatch()
	except (BrokenPipeError, OSError):
		sys.exit(32)


if __name__ == '__main__':
	decoded = main()

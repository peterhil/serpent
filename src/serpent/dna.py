#!/usr/bin/env python

import re

import numpy as np

from more_itertools import chunked

from serpent.digit import digits_to_number
from serpent.fun import map_array
from serpent.padding import pad_to_left


bases = {"A": 0b00, "C": 0b01, "G": 0b10, "T": 0b11, "U": 0b11}
bases_inverse = {0: "A", 1: "C", 2: "G", 3: "T"}


def decode(dna):
	"""Returns dna’s codons encoded into numbers 0..63"""
	return map_array(decode_codon, dna)


def decode_codon(codon):
	result = 0
	for num, char in enumerate(reversed(codon)):
		result += bases.get(char, 0) << num * 2  # Throw IndexError by using []?

	return result


def clean_non_dna(data):
	"""Clean up non DNA or RNA data. Warns if there are residual characters."""
	cleaned = "".join(re.sub(r"[^ACGTU]{6,}", "", data).split("\n"))
	residual = "".join(re.findall(r"[^ACGTU\n]", data))

	if len(residual) > 0:
		# TODO Use logger.warn with warnings.warn?
		print("Residual characters:", residual)

	return cleaned


def get_codons(data):
	"""Get codons from data as Numpy array"""
	codons_list = list(chunked(data, 3, strict=True))
	codons = map_array(lambda c: "".join(c), codons_list, dtype="U3")

	return codons


def codon_sequences(decoded, n=4, fill=0):
	"""Chunk data into length N sequences of codons.
	Count the occurences of different kmers as numbers between 0..64**n.

	Return index and counts.
	"""
	padded = pad_to_left(list(decoded), n, fill)
	sequences = list(chunked(padded, n, strict=True))
	numbers = np.apply_along_axis(digits_to_number, 1, sequences)

	return numbers

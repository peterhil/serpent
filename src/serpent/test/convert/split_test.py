from __future__ import annotations

from serpent.convert.split import split_nucleotides


def test_split_nucleotides_f():
	encoded = ['CCC', 'TAA', 'CCC', 'ATG', 'CCC', 'ATG', 'TAG', 'ATG', 'CCC', 'TAG']
	expected = [
		['CCC', 'TAA', 'CCC'],
		['ATG', 'CCC', 'ATG'],
		['TAG'],
		['ATG', 'CCC'],
		['TAG']
	]
	split = 'f'
	actual = list(split_nucleotides(encoded, table=11, split=split))

	assert actual == expected


def test_split_nucleotides_n():
	encoded = ['CCC', 'TAA', 'CCC', 'ATG', 'CCC', 'ATG', 'TAG', 'ATG', 'CCC', 'TAG']
	expected = [
		['CCC', 'TAA', 'CCC', 'ATG', 'CCC', 'ATG', 'TAG'],
		['ATG', 'CCC', 'TAG']
	]
	split = 'n'
	actual = list(split_nucleotides(encoded, table=11, split=split))

	assert actual == expected


def test_split_nucleotides_r():
	encoded = ['CCC', 'TAA', 'CCC', 'ATG', 'CCC', 'ATG', 'TAG', 'ATG', 'CCC', 'TAG']
	expected = [
		['CCC', 'TAA'],
		['CCC'],
		['ATG', 'CCC'],
		['ATG', 'TAG'],
		['ATG', 'CCC', 'TAG']
	]
	split = 'r'
	actual = list(split_nucleotides(encoded, table=11, split=split))

	assert actual == expected

from __future__ import annotations

from serpent.io import image_name_for


def test_image_name_for_isnt_filename():
	filename = 'sequence.fasta'
	assert image_name_for(filename) != filename


def test_image_name_for():
	expected = 'sequence.faa.w64.GACT.png'
	actual = image_name_for('sequence.faa', 64)
	assert actual == expected


def test_image_name_for_palette():
	expected = 'sequence.faa.Pa.w64.GACT.t11.png'
	actual = image_name_for('sequence.faa', 64, mode='P', amino=True, table=11)
	assert actual == expected

from __future__ import annotations

from serpent.mathematics import autowidth, phi_small


def test_autowidth_square():
	n = 256
	assert autowidth(n ** 2, aspect=1) == n


def test_autowidth_n_less_than_base_squared():
	base = 64
	n = 1365.333  # n < 4096 == base ** 2
	aspect = phi_small  # < 1
	assert autowidth(n, base, aspect) == base

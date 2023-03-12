from __future__ import annotations

from hypothesis import given
from hypothesis import strategies as st

from serpent.mathematics import autowidth, phi_small


def test_autowidth_square():
	n = 256
	assert autowidth(n ** 2, aspect=1) == n


def test_autowidth_n_less_than_base_squared():
	base = 64
	n = 1365.333  # n < 4096 == base ** 2
	aspect = phi_small  # < 1
	assert autowidth(n, base, aspect) == base


@given(st.data())
def test_hypo_autowidth_n_less_than_base_squared(data):
	base = data.draw(st.integers(min_value=4, max_value=256))
	n = data.draw(st.integers(min_value=0, max_value=base**2))
	aspect = 1
	assert autowidth(n, base, aspect) == base


@given(st.integers(min_value=0, max_value=2 ** 32 - 1))
def test_hypo_autowidth_is_multiple_of_base(n):
	base = 64
	aspect=1
	result = autowidth(n, base, aspect)
	assert (result % base) == 0

from __future__ import annotations

import pytest
from hypothesis import given
from hypothesis import strategies as st
from numpy.testing import assert_array_equal

from serpent.math.combinatorics import spread, unspread


@pytest.mark.xfail()
@given(st.data())
def test_spread_unspread(data):
	seq = data.draw(st.lists(
		st.one_of(
			st.integers(min_value=0),
			# st.characters(),
		),
		min_size=1,
		max_size=16,
		unique=True,
	))

	limit = len(seq)
	n = data.draw(st.integers(min_value=1, max_value=len(seq)))
	offset = data.draw(st.integers(min_value=-limit, max_value=limit))

	assert_array_equal(
		unspread(spread(seq, n, offset), n, offset),
		seq,
	)

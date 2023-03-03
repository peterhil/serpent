from __future__ import annotations

from collections.abc import Iterable

import numpy as np


def map_array(fn, arr, dtype=None):
	return np.array(list(map(fn, arr)), dtype=dtype)


def str_join(seq: Iterable) -> str:
	return "".join(seq)

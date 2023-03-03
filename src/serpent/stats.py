"""Statistical function."""
from __future__ import annotations

from collections import Counter

import numpy as np


def count_sorted(items):
	"""Count items and return as sorted array."""
	counts = Counter(items)
	return np.array(sorted(counts.items())).T

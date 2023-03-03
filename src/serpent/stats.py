from __future__ import annotations

from collections import Counter

import numpy as np


def count_sorted(items):
	counts = Counter(items)
	return np.array(sorted(counts.items())).T

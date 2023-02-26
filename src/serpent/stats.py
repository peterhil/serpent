#!/usr/bin/env python

import numpy as np

from collections import Counter


def count_sorted(items):
	counts = Counter(items)
	return np.array(list(sorted(counts.items()))).T

#!/usr/bin/env python

import numpy as np


def map_array(fn, arr, dtype=None):
	return np.array(list(map(fn, arr)), dtype=dtype)
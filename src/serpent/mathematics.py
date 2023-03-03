from __future__ import annotations

import numpy as np


def logn(x, base=np.e):
    """
    Logarithm of x on some base.
    """
    return np.log2(x) / np.log2(base)


def magnitude(data, base=64):
	return int(np.ceil(logn(np.max(data), base)))


def normalise(seq):
	return np.asanyarray(seq) / np.amax(seq)

#!/usr/bin/env python

import numpy as np


def normalise(seq):
	return np.asanyarray(seq) / np.amax(seq)

"""Lab for interactive work."""
# flake8: noqa F401, F403

from __future__ import annotations

import itertools as itr
import re
from collections import Counter, OrderedDict, defaultdict
from functools import partial

import matplotlib.pyplot as plt
import more_itertools as mit
import numpy as np
from PIL import Image

from serpent import *
from serpent.convert.amino import *
from serpent.convert.codon import *
from serpent.convert.degenerate import *
from serpent.convert.digits import *
# from serpent.convert.dnt import *
from serpent.convert.format import *
from serpent.convert.genetic_code import *
from serpent.convert.nucleotide import *
from serpent.convert.quad import *
from serpent.convert.split import *
from serpent.fun import *
from serpent.io.fasta import *
from serpent.io.files import *
from serpent.io.printing import *
from serpent.math.basic import *
from serpent.math.combinatorics import *
from serpent.math.dsp import *
from serpent.math.statistic import *
from serpent.padding import *
from serpent.settings import *
from serpent.spatial.path import *
from serpent.visual import ansi
from serpent.visual.bitmap import *
from serpent.visual.block_elements import *
from serpent.visual.palette import *
from serpent.visual.plot import *

interactive()
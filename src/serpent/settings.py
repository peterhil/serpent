"""Serpent settings."""
from __future__ import annotations

# TODO Read settings from a `.serpentrc` file

# BASE_ORDER notes:
# -----------------
#
# Palette colour impressions:
# (based on Buchnera aphidicola genome with 25.348% GC content)
#
# * Reddish (very interesting & beautiful, clear patterns, esp. on AGCT / TCGA)
# ACGT, AGCT, TCGA, TGCA
# * Orange-bluish-grey (quite beautiful, patterns shown):
# ACTG, AGTC, TCAG, TGAC
# * Yellowish (ugly red-green, patterns disapper into noise):
# ATCG, ATGC, TACG, TAGC
# * Green-purple (ugly colours, but strong patterns)
# CAGT, CTGA, GACT, GTCA
# * Turquoise (very beautiful, patterns visible & look like water waves):
# CATG, CTAG, GATC, GTAC
# * Violet (very beautiful, but patterns not so visible)
# CGAT, CGTA, GCAT, GCTA
#
# Greyscale seems most interesting when TA and GC are close together.

BASE_ORDER = 'GACT'
COUNT_LIMIT = 20
DEBUG = False
DEFAULT_COLOR = '#70e'
DEFAULT_TERM_COLOR = 'magenta'
FLOW_DESCRIPTION_COLOR = (85, 85, 85)

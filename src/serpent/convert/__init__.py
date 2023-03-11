"""Convert data between numbers and encodings."""
# flake8: noqa F401

from __future__ import annotations

from .amino import amino_to_num, num_to_amino
from .base64 import base64_to_num, num_to_base64
from .digits import digits_to_num, num_to_digits
from .nucleotide import nt_to_num, num_to_nt

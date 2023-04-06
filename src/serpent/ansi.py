"""ANSI escape codes module.

See:
https://en.wikipedia.org/wiki/ANSI_escape_code
"""
from __future__ import annotations

from serpent.fun import str_join
from serpent.settings import FLOW_DESCRIPTION_COLOR

CSI = '\x1b['
RESET = f'{CSI}0m'


def sgr(codes):
	return f'{CSI}{codes}m'


def reset():
	return RESET


def colour_code(*rgb, bg=False):
	code = 38 if not bg else 48
	RGB_LEN = 3

	if len(rgb) == 1:
		n = rgb[0]
		return f'{code};5;{n}'
	elif len(rgb) == RGB_LEN:
		(r, g, b) = rgb
		return f'{code};2;{r};{g};{b}'
	else:
		err_msg = f'Invalid rgb values: {rgb}'
		raise ValueError(err_msg)


def rgb(front=None, back=None):
	fg = colour_code(*front) if front else ''
	bg = colour_code(*back, bg=True) if back else ''

	return sgr(str_join([fg, bg], ';'))


def grey(front=None, back=None):
	fg = colour_code(*(front, front, front)) if front else ''
	bg = colour_code(*(back, back, back), bg=True) if back else ''

	return sgr(str_join([fg, bg], ';'))


def dim_text(text):
	return rgb_text(text, front=FLOW_DESCRIPTION_COLOR)


def rgb_text(text, front=(255, 255, 255), back=(0, 0, 0)):
	colour = rgb(front, back)

	return colour + text + RESET

from __future__ import annotations

from collections import defaultdict

import numpy as np

from serpent import dna
from serpent.math.statistic import count_sorted, quasar_pulses
from serpent.printing import format_quasar, format_quasar_pulses
from serpent.settings import PLOT_FONT_SIZE
from serpent.visual.palette import spectrum_layer_colours_for


def pulse_plot(ax, pulse, colour, *, symbol, y_offset=0):
	data = pulse
	ax.plot(data + y_offset, colour)

	ax.text(
		-15, y_offset, symbol, color=colour,
		ha='right', size=PLOT_FONT_SIZE, fontweight='semibold'
	)
	ax.set_xlabel('index of repetition')
	ax.set_ylabel('interval length (space between symbols)')


def pulse_plot_counts(ax, pulse, colour, base: float=10):
	data = count_sorted(pulse)
	ax.plot(*data, colour)

	ax.set_xscale('log', base=base)
	ax.set_yscale('log', base=base)

	ax.set_xlabel('repeat length')
	ax.set_ylabel('count')

	return data


def pulse_plot_symbols_legend(ax, key, colours, width, height):
	mx = len(key) - 1
	for i, aa in enumerate(reversed(key)):
		colour = colours[aa]
		ax.text(
			1 + width * (i / mx),
			2 + height * (1 - i / mx),
			aa,
			color=colour,
			size=PLOT_FONT_SIZE,
			fontweight='semibold'
		)


def pulse_plot_sequences(
	ax, seqs, key,
	*,
	amino=False, table=1, degen=False,
	count=False, cumulative=False,
):
	colours = spectrum_layer_colours_for(key, amino)
	maxheight = maxscale = 0
	sum_pulses = defaultdict(int)

	for seq in seqs:
		[aminos, descriptions] = dna.decode_seq(seq, amino, table, degen, dna.to_amino)
		pulses, height, scale = quasar_pulses(aminos, cumulative=cumulative, key=key)
		maxheight = max(height, maxheight)
		maxscale = max(scale, maxscale)

		for i, item in enumerate(reversed(pulses.items())):
			[aa, pulse] = item

			colour = colours[aa]
			mx = len(pulses) - 1
			pos = mx - i

			if count:
				sorted_pulse = pulse_plot_counts(ax, pulse, colour, base=10)
				for repeat, count in sorted_pulse.T:
					sum_pulses[repeat] += count
			else:
				pulse_plot(ax, pulse, colour, symbol=aa, y_offset=pos * 100)

	if count:
		# sum pulses
		sum_pulses = np.array([*sorted(sum_pulses.items())]).T
		ax.plot(*sum_pulses, color='#2227')

		# symbols legend
		pulse_plot_symbols_legend(ax, key, colours, width=maxscale, height=maxheight)


def pulse_text(pulses, height, scale):
	yield from format_quasar(pulses.keys())  # Print symbols
	yield from format_quasar_pulses(pulses, height)
	yield f'scale: {scale}'

from __future__ import annotations

from serpent.printing import format_quasar, format_quasar_pulses
from serpent.stats import count_sorted


def pulse_plot(ax, pulse, colour, *, symbol, y_offset=0, size=10):
	data = pulse
	ax.plot(data + y_offset, colour)

	# ax.set_xscale('log', base=10)
	# ax.set_yscale('log', base=10)

	ax.text(
		-15, y_offset, symbol, color=colour,
		ha='right', size=size, fontweight='semibold'
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


def pulse_plot_symbols_legend(ax, key, colours, width, height, *, size=10):
	mx = len(key) - 1
	for i, aa in enumerate(reversed(key)):
		colour = colours[aa]
		ax.text(
			1 + width * (i / mx),
			2 + height * (1 - i / mx),
			aa, color=colour,
			size=size, fontweight='semibold'
		)


def pulse_text(pulses, height, scale):
	yield from format_quasar(pulses.keys())  # Print symbols
	yield from format_quasar_pulses(pulses, height)
	yield f'scale: {scale}'

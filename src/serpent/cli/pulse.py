from __future__ import annotations

from serpent.printing import format_quasar, format_quasar_pulses


def pulse_count_symbols_legend(ax, key, colours, width, height, *, size=10):
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

from __future__ import annotations

from serpent.printing import format_quasar, format_quasar_pulses


def pulse_text(pulses, height, scale):
	yield from format_quasar(pulses.keys())  # Print symbols
	yield from format_quasar_pulses(pulses, height)
	yield f'scale: {scale}'

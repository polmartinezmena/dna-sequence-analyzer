"""Render sequence coordinates and translations without a Streamlit dependency."""

from html import escape

from Bio.Seq import Seq

from .translate import translate_frames
from .ui_theme import BASE_COLORS, format_es


def _base_markup(base, complement=False):
    color = BASE_COLORS.get(base, BASE_COLORS["other"])
    classes = "base complement" if complement else "base"
    return f'<span class="{classes}" style="--base:{color}">{escape(base)}</span>'


def render_sequence_html(sequence, start, end, frames=None, table=1,
                         show_complement=True, rna=False):
    """Render a 1-based inclusive interval in 60-base blocks."""
    sequence = sequence.upper()
    if not 1 <= start <= end <= len(sequence):
        raise ValueError("La región debe estar dentro de la secuencia.")
    if frames is None:
        frames = translate_frames(sequence, table)
    display = sequence.replace("T", "U") if rna else sequence
    fragment = display[start - 1:end]
    complement = str(Seq(fragment).complement_rna() if rna else Seq(fragment).complement())
    blocks = []
    for block_start in range(start, end + 1, 60):
        block_end = min(end, block_start + 59)
        positions = range(block_start, block_end + 1)
        ticks = "".join(
            f'<span class="tick" data-position="{pos}">{format_es(pos)}</span>'
            if pos == block_start or pos % 10 == 0 else "<span></span>"
            for pos in positions
        )
        rows = [f'<div class="sequence-row"><span class="track">Posición</span>{ticks}</div>']
        forward = "".join(_base_markup(display[pos - 1]) for pos in positions)
        rows.append(f'<div class="sequence-row"><span class="track">5′→3′</span>{forward}</div>')
        if show_complement:
            reverse = "".join(_base_markup(complement[pos - start], True) for pos in positions)
            rows.append(f'<div class="sequence-row"><span class="track">3′→5′</span>{reverse}</div>')
        for frame in frames:
            offset = abs(int(frame.frame)) - 1
            reverse = frame.frame.startswith("-")
            cells = []
            for pos in positions:
                distance = len(sequence) - pos if reverse else pos - 1
                index, phase = divmod(distance - offset, 3)
                if phase != 1 or not 0 <= index < len(frame.protein):
                    cells.append("<span></span>")
                    continue
                residue = frame.protein[index]
                tone = {"M": "met", "*": "stop"}.get(residue, "residue")
                cells.append(f'<span class="{tone}" data-frame="{escape(frame.frame)}" '
                             f'data-codon-center="{pos}">{escape(residue)}</span>')
            rows.append(f'<div class="sequence-row"><span class="track">{escape(frame.frame)}</span>'
                        f'{"".join(cells)}</div>')
        blocks.append(f'<div class="sequence-scroll">{"".join(rows)}</div>')
    return "".join(blocks)

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


def render_base_tiles_html(sequence, start, end, record_id, rna=False):
    """Render aligned ruler marks and base tiles for an inclusive interval."""
    if not 1 <= start <= end <= len(sequence):
        raise ValueError("La región debe estar dentro de la secuencia.")
    fragment = sequence[start - 1:end]
    cells = []
    for position, symbol in enumerate(fragment, start):
        base = symbol.upper()
        if rna and base == "T":
            base = "U"
        tick = "tile-tick mark" if position == start or position % 10 == 0 else "tile-tick"
        color = BASE_COLORS.get(base, BASE_COLORS["other"])
        tone = " tile-ambiguous" if base not in BASE_COLORS else ""
        char = escape(base)
        cells.append(
            f'<div class="tile-cell"><span class="{tick}" data-position="{position}">'
            f'{format_es(position)}</span>'
            f'<span class="tile-base{tone}" data-base="{char}" data-position="{position}" '
            f'style="--base:{color}">{char}</span></div>'
        )
    return (
        '<div class="base-tiles"><p class="tile-ruler-label">Posición (pb)</p>'
        f'<p class="tile-record">{escape(record_id)}</p>'
        f'<div class="tile-grid">{"".join(cells)}</div></div>'
    )


def render_orf_track_html(hits, start, end):
    """Render the visible portion of each overlapping ORF on a strand-labelled track."""
    length = end - start + 1
    tracks = []
    for hit in hits:
        left = (max(hit.start, start) - start) / length * 100
        width = (min(hit.end, end) - max(hit.start, start) + 1) / length * 100
        if width <= 0:
            continue
        frame = escape(hit.frame.replace("-", "−"))
        label = f"Marco {frame}, {format_es(hit.start)}–{format_es(hit.end)} pb"
        tracks.append(
            f'<div class="orf-track"><span>{frame}</span><div class="orf-track-rail">'
            f'<span aria-label="{label}" title="{label}" role="img" class="orf-track-bar" '
            f'style="left:{left:.6f}%;width:{width:.6f}%"></span></div></div>'
        )
    return f'<div class="orf-tracks">{"".join(tracks)}</div>'

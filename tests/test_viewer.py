from dna_sequence_analyzer.translate import translate_frames
from dna_sequence_analyzer.viewer import render_sequence_html


def test_rule_uses_absolute_positions_and_wraps():
    html = render_sequence_html("ACGT" * 40, 7, 80)
    assert 'data-position="10"' in html
    assert 'data-position="70"' in html
    assert html.count('class="sequence-scroll"') == 2


def test_codon_centers_for_all_frames():
    seq = "ATGAAATAGCCT"
    html = render_sequence_html(seq, 1, len(seq), translate_frames(seq))
    expected = {"+1": (2, "M"), "+2": (3, "*"), "+3": (4, "E"),
                "-1": (11, "R"), "-2": (10, "G"), "-3": (9, "A")}
    for frame, (position, residue) in expected.items():
        assert f'data-frame="{frame}" data-codon-center="{position}">{residue}' in html


def test_markup_is_escaped_and_rna_uses_uracil():
    html = render_sequence_html("AT<x", 1, 4, frames=[], show_complement=True, rna=True)
    assert "&lt;" in html and "<x" not in html
    assert ">U</span>" in html and ">T</span>" not in html


def test_tile_ruler_matches_bases():
    import re

    from dna_sequence_analyzer.ui_theme import BASE_COLORS
    from dna_sequence_analyzer.viewer import render_base_tiles_html

    sequence = "ATGC" * 8
    html = render_base_tiles_html(sequence, 5, 21, "demo")
    assert re.findall(r'class="tile-tick mark" data-position="(\d+)"', html) == ["5", "10", "20"]
    cells = re.findall(r'data-base="(.)" data-position="(\d+)" style="--base:(#[A-Fa-f0-9]+)"', html)
    assert cells == [(base, str(i), BASE_COLORS[base]) for i, base in enumerate(sequence[4:21], 5)]


def test_tile_html_escapes_record_and_converts_rna():
    from dna_sequence_analyzer.ui_theme import BASE_COLORS
    from dna_sequence_analyzer.viewer import render_base_tiles_html

    html = render_base_tiles_html("AT<", 1, 3, '<script title="x">', rna=True)
    assert "<script" not in html
    assert "&lt;script title=&quot;x&quot;&gt;" in html
    assert 'data-base="U" data-position="2"' in html
    assert f'data-base="&lt;" data-position="3" style="--base:{BASE_COLORS["other"]}"' in html

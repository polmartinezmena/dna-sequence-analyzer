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

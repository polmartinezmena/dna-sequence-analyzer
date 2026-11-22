from dna_sequence_analyzer.ui_theme import BASE_COLORS, format_es, stylesheet
from dna_sequence_analyzer.plots import lab_template


def test_spanish_numbers():
    assert format_es(1387, unit="pb") == "1.387 pb"
    assert format_es(63.875, 2, "%") == "63,88 %"
    assert format_es(-1234.5, 1) == "-1.234,5"
    assert format_es(None, 2, "%") == "—"


def test_theme_pairs():
    css = stylesheet()
    assert BASE_COLORS["U"] == BASE_COLORS["T"]
    assert "#FFFFFF" in css and "#0F1318" in css


def test_plotly_template_uses_both_palettes():
    from plotly import io

    light = io.templates[lab_template("light")].layout
    dark = io.templates[lab_template("dark")].layout
    assert light.font.color != dark.font.color
    assert light.height == dark.height == 320
    assert light.separators == dark.separators == ",."

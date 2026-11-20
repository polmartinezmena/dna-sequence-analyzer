from pathlib import Path

from streamlit.testing.v1 import AppTest

from dna_sequence_analyzer.lab import quality, report, run_metadata
from dna_sequence_analyzer.web_input import describe_records, parse_sequence_text


def test_quality():
    records = parse_sequence_text(">empty\n>bad\nGCN!\n")
    checks = quality(records, describe_records(records))
    assert [level for level, _ in checks] == ["Error", "Error", "Aviso", "Aviso"]
    clean = parse_sequence_text("ATGC")
    assert quality(clean, describe_records(clean))[0][0] == "OK"


def test_run_metadata():
    records = parse_sequence_text("ATGC")
    first = run_metadata(records, {"Paso": 1})
    assert run_metadata(records, {"Paso": 1}, first) == first
    assert run_metadata(records, {"Paso": 2}, first)["ID"] != first["ID"]
    assert run_metadata(parse_sequence_text("AAAA"), {"Paso": 1}, first)["ID"] != first["ID"]


def test_report():
    html = report({"ID": "run-1"}, [("Aviso", "<script>")], {"GC": [{"GC %": 42}]})
    assert all(value in html for value in ["run-1", "Aviso", "&lt;script&gt;", "42"])


def test_app_example():
    app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app/streamlit_app.py"))
    app.run(timeout=30)
    assert not app.exception
    app.button[0].click().run(timeout=30)
    assert not app.exception
    assert len(app.tabs) == 6
    assert app.dataframe

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
    assert app.session_state["records"][0].id == "NM_000207.3"
    app.button(key="nav-Traducción").click().run()
    assert not app.exception
    assert any('<pre class="lab-protein">' in entry.proto.body for entry in app.get("html"))
    app.button(key="nav-ORFs").click().run()
    next(field for field in app.number_input if field.label == "Longitud mínima (aa)").set_value(1)
    app.button(key="buscar-orfs").click().run()
    assert not app.exception
    assert app.session_state["orf_result"][2]


def test_global_gc_is_weighted_and_excludes_ambiguous_bases():
    from dna_sequence_analyzer.lab import global_gc_percent
    from dna_sequence_analyzer.web_input import describe_records, parse_sequence_text

    stats = describe_records(parse_sequence_text(">corta\nGCNN\n>larga\nAAAAAAAA"))
    assert global_gc_percent(stats) == 20
    assert global_gc_percent([]) is None

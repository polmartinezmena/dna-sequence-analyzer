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
    assert app.selectbox(key="source").value == "Subir un FASTA"
    assert "Visor" not in [tab.label for tab in app.tabs]
    app.selectbox(key="source").set_value("Colección INS y GFP").run(timeout=30)
    next(field for field in app.selectbox if field.label == "Registro en el visor").select("L29345.1").run(timeout=30)
    assert any("<h2>L29345.1</h2>" in item.proto.body for item in app.get("html"))
    assert not app.exception
    assert "ORFs" in [tab.label for tab in app.tabs]
    next(field for field in app.number_input if field.label == "Longitud mínima (aa)").set_value(1)
    next(button for button in app.button if button.label == "Buscar ORFs").click().run(timeout=30)
    assert not app.exception
    assert len(app.get("download_button")) >= 4
    assert app.dataframe

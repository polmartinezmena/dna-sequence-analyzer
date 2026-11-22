from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest  # noqa: E402

from dna_sequence_analyzer import blast  # noqa: E402

LAUNCHER = Path(__file__).resolve().parents[1] / "app" / "streamlit_app.py"


@pytest.fixture
def app(monkeypatch):
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    monkeypatch.setattr(blast, "_request", lambda *args: pytest.fail("Petición real a NCBI"))
    return AppTest.from_file(str(LAUNCHER), default_timeout=60).run()


def open_view(app, view):
    app.button(key=f"nav-{view}").click().run()
    assert not app.exception


def test_navigation_with_example(app):
    assert not app.exception
    assert app.session_state["vista"] == "Visor"
    assert app.session_state["nombre_archivo"] == "examples.fasta"
    for view in ("Resumen", "Visor", "Traducción", "ORFs", "BLAST", "Oligo Tm", "Aprende"):
        open_view(app, view)
        assert any(view in heading.value for heading in app.title)


def test_default_orfs_for_insulin(app):
    open_view(app, "ORFs")
    app.button(key="buscar-orfs").click().run()
    assert not app.exception
    assert len(app.dataframe[0].value) == 2
    assert app.session_state["metadata"]["Inicios alternativos"] is False


def test_blast_disabled_without_email(app):
    open_view(app, "BLAST")
    assert app.button(key="lanzar-blast").disabled


def test_import_pasted_fasta_with_preamble(app):
    app.button(key="pegar-muestra").click().run()
    app.text_area(key="texto-fasta").set_value(";comentario\n>uno\nACGT\n>dos\nGC")
    app.button(key="cargar-muestra").click().run()
    assert not app.exception
    assert app.session_state["nombre_archivo"] == "pegado.fasta"
    assert [r.id for r in app.session_state["records"]] == ["uno", "dos"]
    assert any("se ignoraron 1" in element.proto.body for element in app.get("html"))


def html_elements(app):
    return [element.proto.body for element in app.get("html")]


def apply_region(app, start, end):
    next(field for field in app.number_input if field.label == "Inicio").set_value(start)
    next(field for field in app.number_input if field.label == "Fin").set_value(end)
    next(button for button in app.button if button.label == "Aplicar").click().run()
    assert not app.exception


def test_visor_region_and_validation(app):
    apply_region(app, 5, 40)
    html = next(body for body in html_elements(app) if 'class="base-tiles"' in body)
    assert html.split('data-position="', 1)[1].startswith('5"')
    assert 'data-position="40"' in html and 'data-position="41"' not in html
    apply_region(app, 5, 9999)
    assert any("longitud" in error.value for error in app.error)
    assert app.session_state["region"][0] == (5, 40)
    app.button(key="registro-1").click().run()
    apply_region(app, 1, 503)
    assert any("500 pb" in error.value for error in app.error)


def test_visor_switch_record_preserves_each_region(app):
    apply_region(app, 5, 40)
    app.button(key="registro-1").click().run()
    assert not app.exception
    assert app.session_state["registro_activo"] == 1
    assert any('class="sample-filename">examples.fasta' in body for body in html_elements(app))
    assert any("L29345.1" in body and "Registro seleccionado" in body for body in html_elements(app))
    app.button(key="registro-0").click().run()
    assert next(field for field in app.number_input if field.label == "Inicio").value == 5


def test_table_selection_changes_active_record(app):
    app.session_state["tabla-registros"] = {"selection": {"rows": [1], "columns": [], "cells": []}}
    app.run()
    assert not app.exception
    assert app.session_state["registro_activo"] == 1
    assert any("L29345.1" in body and "Registro seleccionado" in body for body in html_elements(app))
    app.button(key="registro-0").click().run()
    assert app.session_state["registro_activo"] == 0
    app.run()
    assert app.session_state["registro_activo"] == 0


def test_blast_limits_block_forged_clicks(app, monkeypatch):
    from time import monotonic

    monkeypatch.setenv("NCBI_EMAIL", "test@example.org")
    open_view(app, "BLAST")
    next(box for box in app.checkbox if "autorizo" in box.label).check().run()
    click = app.button(key="lanzar-blast")
    assert not click.disabled
    app.session_state["last_blast_search"] = monotonic()
    click.click().run()
    assert not app.exception
    assert app.button(key="lanzar-blast").disabled
    assert any("próxima búsqueda BLAST" in info.value for info in app.info)
    assert "blast_result" not in app.session_state


def test_oversized_blast_query(app, monkeypatch):
    monkeypatch.setenv("NCBI_EMAIL", "test@example.org")
    app.button(key="pegar-muestra").click().run()
    app.text_area(key="texto-fasta").set_value(">larga\n" + "ACGT" * 2501)
    app.button(key="cargar-muestra").click().run()
    open_view(app, "BLAST")
    assert app.button(key="lanzar-blast").disabled
    assert any("10.000 nt" in error.value for error in app.error)


def test_blast_result_and_cooldown_use_mocked_search(app, monkeypatch):
    from dna_sequence_analyzer.app.views import blast as blast_view

    calls = []

    def fake_search(query, options, **kwargs):
        calls.append(options)
        return [blast.BlastHit("TEST.1", "Control", 97.5, 80, 1e-20, 200)]

    monkeypatch.setenv("NCBI_EMAIL", "test@example.org")
    monkeypatch.setattr(blast_view, "search", fake_search)
    open_view(app, "BLAST")
    next(box for box in app.checkbox if "autorizo" in box.label).check().run()
    app.button(key="lanzar-blast").click().run()
    assert not app.exception
    assert len(calls) == 1 and calls[0].megablast
    assert app.button(key="lanzar-blast").disabled
    assert app.dataframe[0].value.iloc[0]["accession"] == "TEST.1"
    assert app.session_state["metadata"]["MegaBLAST"] is True


def download_labels(app):
    return {button.proto.label for button in app.get("download_button")}


def test_visor_segments_annotations_and_exports(app):
    assert {"Informe HTML", "CSV de la vista actual", "Metadatos JSON"} <= download_labels(app)
    app.session_state["pestana-visor"] = "Vista por tramos"
    app.run()
    assert not app.exception
    assert any('data-frame="+3"' in body for body in html_elements(app))
    app.session_state["pestana-visor"] = "Anotaciones"
    app.run()
    assert not app.exception
    assert any('class="orf-track-bar"' in body for body in html_elements(app))
    assert {"ORFs FASTA", "ORFs GFF3"} <= download_labels(app)
    apply_region(app, 1, 20)
    app.session_state["pestana-visor"] = "Anotaciones"
    app.run()
    assert any("No hay ORFs" in message.value for message in app.info)


def test_protein_download_is_in_sample_exports(app):
    open_view(app, "Traducción")
    assert "Proteína FASTA" in download_labels(app)


def test_viewer_gc_uses_centered_windows_and_full_record_reference(app):
    import json

    app.button(key="pegar-muestra").click().run()
    app.text_area(key="texto-fasta").set_value(">demo\nAGNC")
    app.button(key="cargar-muestra").click().run()
    app.number_input(key="ventana-visor").set_value(2).run()
    assert not app.exception
    spec = json.loads(app.get("plotly_chart")[0].proto.spec)
    assert spec["data"][0]["x"] == [1.5, 2.5, 3.5]
    assert spec["data"][0]["y"] == [50.0, 100.0, 100.0]
    assert spec["layout"]["shapes"][0]["y0"] == pytest.approx(200 / 3)

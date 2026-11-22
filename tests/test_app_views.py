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

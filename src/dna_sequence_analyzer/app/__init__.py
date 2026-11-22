"""Streamlit workspace and view routing."""

import streamlit as st

from dna_sequence_analyzer.lab import quality, run_metadata
from dna_sequence_analyzer.ui_theme import stylesheet
from dna_sequence_analyzer.web_input import describe_records

from . import layout
from .views import aprende, blast, oligo, orfs, resumen, traduccion, visor

VIEWS = {
    "Resumen": resumen, "Visor": visor, "Traducción": traduccion,
    "ORFs": orfs, "BLAST": blast, "Oligo Tm": oligo, "Aprende": aprende,
}


def main():
    st.set_page_config(page_title="Analizador de secuencias", layout="wide")
    state = st.session_state
    view = state.setdefault("vista", "Visor")
    st.html(stylesheet(active_view=view))
    if "records" not in state:
        try:
            layout.load_example()
        except (OSError, ValueError) as exc:
            st.error(f"No se pudo cargar el ejemplo: {exc}")
            if st.button("Importar FASTA", icon=":material/upload:"):
                layout.import_dialog("Subir archivo")
            return
    records = state.records
    if view == "Visor":
        visor.follow_table_selection(state, len(records))
    active = state.registro_activo
    record = records[active]
    stats = describe_records(records)
    checks = quality(records, stats, comments=state.fasta_comments)
    footer = layout.sidebar(records, active, view)
    layout.topbar(view)
    params = {"Origen": state.source, "Archivo": state.nombre_archivo, "Registro activo": record.id, "Vista": view}
    previous = state.get("metadata")
    popover = layout.sample_header(stats)

    st.title(layout.NAVIGATION[view][0])
    results, exports = {}, {}
    with st.container(key="espacio-vista"):
        body, side = st.columns([3, 1], gap="medium")
        with side:
            settings = st.container(border=True, key="panel-lateral" if view == "Visor" else "parametros")
            if view != "Visor":
                settings.html('<h2 class="section-label">Parámetros</h2>')
        with body:
            if view not in ("Resumen", "Visor"):
                body = st.container(border=True, key="tarjeta-contenido")
            with body:
                VIEWS[view].render(record, stats[active], stats, checks, results, params, settings, exports)
    if view == "Visor":
        visor.records_table(stats, results)
    metadata = run_metadata(records, params, previous)
    state.metadata = metadata
    saved_orfs = state.get("orf_result")
    if saved_orfs and saved_orfs[2] and saved_orfs[0] == orfs.record_key(record):
        layout.orf_downloads(exports, saved_orfs[2])
    layout.downloads(popover, metadata, checks, results, exports)
    layout.run_footer(footer, metadata)
    if mode := state.get("import_mode"):
        layout.import_dialog(mode)

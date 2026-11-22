"""Shared layout, sample import and downloads."""

import json
from datetime import datetime, timezone
from html import escape
from pathlib import Path

import streamlit as st

from dna_sequence_analyzer.lab import global_gc_percent, report, result_csv
from dna_sequence_analyzer.ui_theme import format_es, stylesheet
from dna_sequence_analyzer.web_input import (
    bundled_example_fasta, looks_like_fasta, parse_sequence_text, read_fasta_text, text_from_upload,
)

NAVIGATION = {
    "Resumen": ("Resumen", "dashboard"),
    "Visor": ("Visor de secuencia", "view_week"),
    "Traducción": ("Traducción", "translate"),
    "ORFs": ("ORFs", "linear_scale"),
    "BLAST": ("BLAST", "travel_explore"),
    "Oligo Tm": ("Oligo Tm", "thermostat"),
    "Aprende": ("Aprende", "school"),
}
IMPORT_MODES = ("Subir archivo", "Pegar texto", "Ejemplo")


def load_sample(text, name, source):
    if looks_like_fasta(text) or source != "Pegar texto":
        records, comments = read_fasta_text(text)
    else:
        records, comments = parse_sequence_text(text), 0
    if not records:
        raise ValueError("El archivo no contiene registros FASTA.")
    st.session_state.update(
        records=records, nombre_archivo=name, source=source, fasta_comments=comments,
        sample_date=datetime.now(timezone.utc).strftime("%d/%m/%Y"),
        registro_activo=0, region={}, table_selection=tuple(
            st.session_state.get("tabla-registros", {}).get("selection", {}).get("rows", ())
        ),
        table_sample=st.session_state.get("table_sample", 0) + 1,
    )
    for key in ("orf_result", "blast_result", "blast_error"):
        st.session_state.pop(key, None)


def load_example():
    path = bundled_example_fasta()
    load_sample(path.read_text(encoding="utf-8-sig"), path.name, "Ejemplo")


def close_import():
    st.session_state.pop("import_mode", None)


@st.dialog("Importar FASTA", on_dismiss=close_import)
def import_dialog(mode):
    source = st.segmented_control("Origen del FASTA", IMPORT_MODES, default=mode, key=f"origen-{mode}") or mode
    upload = text = None
    if source == "Subir archivo":
        upload = st.file_uploader("Archivo FASTA", type=["fasta", "fa", "fna", "txt"])
    elif source == "Pegar texto":
        text = st.text_area("FASTA o secuencia", height=200, placeholder=">mi_secuencia\nATGC",
                            key="texto-fasta")
    else:
        st.caption("Carga los registros de insulina humana y proteína fluorescente del ejemplo.")
    if st.button("Cargar muestra", type="primary", key="cargar-muestra"):
        try:
            if source == "Ejemplo":
                load_example()
            elif source == "Pegar texto":
                load_sample(text, "pegado.fasta", source)
            elif upload is not None:
                load_sample(text_from_upload(upload.getvalue()), Path(upload.name).name, source)
            else:
                raise ValueError("Selecciona un archivo FASTA.")
        except (OSError, ValueError) as exc:
            st.error(str(exc))
        else:
            close_import()
            st.rerun()


def sidebar(records, active, view):
    with st.sidebar:
        with st.container(key="muestra"):
            title, folder, add = st.columns([3, 1, 1], vertical_alignment="center", gap="small")
            title.html('<h2 class="section-label">Muestra</h2>')
            if folder.button("", key="abrir-muestra", icon=":material/folder_open:", help="Importar FASTA"):
                st.session_state.import_mode = "Subir archivo"
            if add.button("", key="pegar-muestra", icon=":material/add:", help="Pegar texto"):
                st.session_state.import_mode = "Pegar texto"
            st.html(
                '<p class="sample-tree"><span class="file-icon" aria-hidden="true">description</span>'
                f'<span><strong>{escape(st.session_state.nombre_archivo)}</strong> ({format_es(len(records))})</span></p>'
            )
            for index, record in enumerate(records):
                label = record.id.replace("`", "′")
                st.button(f"`{label}`", key=f"registro-{index}", width="stretch",
                          type="primary" if index == active else "tertiary",
                          on_click=select_record, args=(index,))
        with st.container(key="navegacion"):
            for index, (name, (_, icon)) in enumerate(NAVIGATION.items()):
                if index == 5:
                    st.divider()
                st.button(name, key=f"nav-{name}", type="tertiary", icon=f":material/{icon}:",
                          width="stretch", on_click=set_view, args=(name,),
                          help="Vista activa" if name == view else None)
        return st.container(key="pie-lateral")


def set_view(view):
    st.session_state.vista = view


def select_record(index):
    st.session_state.registro_activo = index


def topbar(view):
    title, _ = NAVIGATION[view]
    st.html(f'<div class="topbar"><strong>DNA Sequence Analyzer</strong>'
            f'<span class="topbar-separator">|</span><span>{escape(title)}</span></div>')


def sample_header(stats):
    name = st.session_state.nombre_archivo
    initials = (Path(name).stem[:2] or "FA").upper()
    date = st.session_state.sample_date
    with st.container(border=True, key="cabecera-muestra"):
        identity, count, bases, gc, actions = st.columns(
            [2.8, 1, 1.15, 1, 2.45], gap="medium", vertical_alignment="center"
        )
        identity.html(
            '<div class="sample-identity">'
            f'<span aria-hidden="true" class="sample-avatar">{escape(initials)}</span>'
            f'<div><p class="sample-date">{date} UTC</p>'
            f'<p class="sample-filename">{escape(name)}</p></div></div>'
        )
        for column, label, value in (
            (count, "Registros", format_es(len(stats))),
            (bases, "Bases totales", format_es(sum(item.length for item in stats), unit="pb")),
            (gc, "GC medio", format_es(global_gc_percent(stats), 1, "%")),
        ):
            column.html(f'<div class="sample-stat"><strong>{value}</strong><small>{label}</small></div>')
        with actions:
            upload, export = st.columns([1.2, 1], gap="small")
            if upload.button("Importar FASTA", key="importar-cabecera", icon=":material/upload:"):
                st.session_state.import_mode = "Subir archivo"
            return export.popover("Exportar", type="primary", icon=":material/download:")


def downloads(popover, metadata, checks, results, exports):
    run = metadata["ID"][:12]
    document = report(metadata, checks, results)
    document = document.replace("<title>", stylesheet() + "<title>", 1)
    document = document.replace("<h1>", '<main class="lab-report"><h1>', 1).replace("</html>", "</main></html>")
    with popover:
        st.download_button("Informe HTML", document, f"analisis-{run}.html", "text/html")
        available = [title for title, rows in results.items() if rows]
        if available:
            selected = st.selectbox("Tabla de resultados", available, key="exportar-tabla")
            st.download_button("CSV de la vista actual", result_csv(results[selected]),
                               f"{selected}-{run}.csv", "text/csv")
        for label, (data, filename, mime) in exports.items():
            st.download_button(label, data, filename, mime)
        st.download_button("Metadatos JSON", json.dumps(metadata, ensure_ascii=False, indent=2),
                           f"metadatos-{run}.json", "application/json")


def run_footer(placeholder, metadata):
    placeholder.html(f'<div class="lab-run">Run {escape(metadata["ID"][:12])} · '
                     f'SHA-256 {escape(metadata["SHA-256 de la muestra"][:12])}</div>')

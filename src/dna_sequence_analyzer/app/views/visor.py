"""Regional sequence viewer, annotations and sample quality."""

from html import escape
from inspect import signature

import plotly.graph_objects as go
import streamlit as st

from dna_sequence_analyzer.lab import global_gc_percent
from dna_sequence_analyzer.orfs import find_orfs
from dna_sequence_analyzer.plots import lab_template
from dna_sequence_analyzer.ui_theme import CHART_LINE, INDIGO_SOFT, REFERENCE, format_es, spanish_table
from dna_sequence_analyzer.viewer import render_base_tiles_html, render_orf_track_html, render_sequence_html
from dna_sequence_analyzer.windows import gc_windows_for_records

_PROGRAMMATIC_SELECTION = "selection_default" in signature(st.dataframe).parameters
_TABLE_KEY = "tabla-registros"
_STATUS_ICONS = {"OK": "check_circle", "Aviso": "warning", "Error": "error"}
MAX_REGION = 500


def render(record, item, stats, checks, results, params, settings, exports):
    with settings:
        _quality_panel(item, stats, checks)
        if item.length:
            start, end = _region(item.length)
        else:
            st.info("El registro está vacío.")
            return
    params.update({"Inicio (pb)": start, "Fin (pb)": end})
    sequence = str(record.seq)
    with st.container(key="panel-oscuro"):
        tiles, segments, annotations = st.tabs(
            ["Secuencia", "Vista por tramos", "Anotaciones"], key="pestana-visor", on_change="rerun"
        )
        if tiles.open:
            with tiles:
                st.html(render_base_tiles_html(sequence, start, end, record.id))
                _gc_plot(record, item, start, end, results, params)
        if segments.open:
            with segments:
                try:
                    html = render_sequence_html(sequence, start, end)
                except ValueError as exc:
                    st.warning(str(exc))
                    html = render_sequence_html(sequence, start, end, frames=[])
                st.html(html)
        if annotations.open:
            with annotations:
                _annotations(sequence, start, end, results)


def records_table(stats, results):
    rows = [
        {"Identificador": item.id, "Longitud (pb)": item.length,
         "GC (%)": item.gc_percent, "Bases ambiguas": item.ambiguous_bases}
        for item in stats
    ]
    results["Registros"] = rows
    active = st.session_state.registro_activo
    with st.container(border=True, key="tarjeta-registros"):
        st.html('<h2 class="section-label">Registros</h2>')
        styler = spanish_table(rows, {"GC (%)": 1})
        selection_options = {}
        if _PROGRAMMATIC_SELECTION:
            selection = {"selection": {"rows": [active]}}
            selection_options["selection_default"] = selection
            if _table_rows(st.session_state) != (active,):
                st.session_state[_TABLE_KEY] = selection
                st.session_state.table_selection = (active,)
        else:
            styler.apply(
                lambda row: [f"background-color: {INDIGO_SOFT}" if row.name == active else ""] * len(row),
                axis=1,
            )
        st.dataframe(styler, on_select="rerun", selection_mode="single-row", hide_index=True,
                     key=_TABLE_KEY, width="stretch", **selection_options)


def _table_rows(state):
    return tuple(state.get(_TABLE_KEY, {}).get("selection", {}).get("rows", ()))


def follow_table_selection(state, count):
    """Apply a new table selection without overriding a subsequent sidebar click."""
    selected = _table_rows(state)
    if selected != state.get("table_selection", ()):
        state.table_selection = selected
        if selected and 0 <= selected[0] < count:
            state.registro_activo = selected[0]


def _quality_panel(item, stats, checks):
    summary = (
        ("Archivo", st.session_state.nombre_archivo),
        ("Registros", format_es(len(stats))),
        ("GC", format_es(global_gc_percent(stats), 1, "%")),
        ("Bases ambiguas", format_es(sum(entry.ambiguous_bases for entry in stats))),
    )
    details = (
        ("Identificador", item.id),
        ("Longitud", format_es(item.length, unit="pb")),
        ("GC", format_es(item.gc_percent, 1, "%")),
        ("Bases ambiguas", format_es(item.ambiguous_bases)),
    )
    lines = "".join(_quality_line("OK", f"{name}: {value}") for name, value in summary)
    lines += "".join(_quality_line(level, f"{level}: {message}") for level, message in checks)
    entries = "".join(
        f'<div><dt>{name}</dt><dd class="{"record-id" if name == "Identificador" else ""}">'
        f'{escape(value)}</dd></div>' for name, value in details
    )
    st.html(f'<h2 class="section-label">Control de calidad</h2><ul class="qc-list">{lines}</ul>'
            f'<h2 class="section-label record-heading">Registro seleccionado</h2>'
            f'<dl class="record-details">{entries}</dl>')


def _quality_line(level, message):
    return (f'<li class="qc-{level.lower()}"><span aria-hidden="true" class="qc-icon">'
            f'{_STATUS_ICONS[level]}</span><span>{escape(message)}</span></li>')


def _region(length):
    state = st.session_state
    index = state.registro_activo
    start, end = state.region.get(index, (1, min(100, length)))
    start, end = min(start, length), min(end, length)
    state.region[index] = start, end
    st.html('<h2 class="section-label record-heading">Región (1-based)</h2>')
    with st.form(key="form-region", border=False):
        left, dash, right = st.columns([5, 1, 5], gap=None, vertical_alignment="center")
        key = f"{state.table_sample}-{index}"
        requested_start = left.number_input("Inicio", value=start, step=1, key=f"inicio-{key}")
        dash.html('<div class="region-dash" aria-hidden="true">–</div>')
        requested_end = right.number_input("Fin", value=end, step=1, key=f"fin-{key}")
        if st.form_submit_button("Aplicar", type="primary", width="stretch"):
            if not 1 <= requested_start <= requested_end <= length:
                st.error(f"Usa 1 ≤ Inicio ≤ Fin ≤ longitud del registro ({format_es(length, unit='pb')}).")
            elif requested_end - requested_start + 1 > MAX_REGION:
                st.error(f"La región puede ocupar hasta {format_es(MAX_REGION, unit='pb')}.")
            else:
                start, end = requested_start, requested_end
                state.region[index] = start, end
    return start, end


def _gc_plot(record, item, start, end, results, params):
    st.markdown("**GC por posición (%)**")
    size = end - start + 1
    window = st.number_input("Ventana (pb)", min_value=1, max_value=size, value=min(10, size),
                             key=f"ventana-visor-{size}", width=148)
    params["Ventana del visor (pb)"] = window
    rows = gc_windows_for_records([record[start - 1:end]], window, 1)
    positions = [start - 1 + (row.start + row.end) / 2 for row in rows]
    values = [row.gc_percent for row in rows]
    results["GC por posición"] = [
        {"Posición (pb)": position, "GC (%)": value} for position, value in zip(positions, values)
    ]
    figure = go.Figure(go.Scatter(
        x=positions, y=values, mode="lines" if len(rows) > 1 else "markers",
        line={"color": CHART_LINE}, name="GC de la ventana",
        hovertemplate="%{x:,.1f} pb<br>%{y:.1f} %<extra></extra>",
    ))
    if item.gc_percent is not None:
        figure.add_hline(y=item.gc_percent, line_dash="dash", line_color=REFERENCE,
                         annotation_text=f"GC del registro · {format_es(item.gc_percent, 1, '%')}",
                         annotation_position="bottom right", annotation_font_color=REFERENCE)
    figure.update_layout(
        template=lab_template("dark"), height=260, margin={"t": 8, "l": 42, "r": 12, "b": 38},
        xaxis={"title": "Posición (pb)", "range": [start - 0.5, end + 0.5], "separatethousands": True},
        yaxis={"range": [0, 100], "ticksuffix": " %"},
    )
    st.plotly_chart(figure, theme=None, width="stretch", config={"displayModeBar": False})


def _annotations(sequence, start, end, results):
    try:
        hits = [hit for hit in find_orfs(sequence, min_aa=30) if hit.start <= end and hit.end >= start]
    except ValueError as exc:
        st.warning(str(exc))
        return
    if not hits:
        st.info("No hay ORFs ≥ 30 aa en esta región.")
        return
    st.html(render_orf_track_html(hits, start, end))
    results["Anotaciones"] = [
        {"Marco": hit.frame.replace("-", "−"), "Hebra": hit.strand.replace("-", "−"),
         "Inicio": hit.start, "Fin": hit.end, "Longitud (aa)": len(hit.protein)}
        for hit in hits
    ]
    st.dataframe(spanish_table(results["Anotaciones"]), hide_index=True, width="stretch")

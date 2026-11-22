"""Streamlit interface for the DNA sequence analyzer.

Run it with ``dna-analyzer web`` or ``streamlit run app/streamlit_app.py``.
Every number comes from the same functions as the command line.
"""

from __future__ import annotations

import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from html import escape
import re
from dataclasses import asdict
from dna_sequence_analyzer.lab import quality, run_metadata, report, result_csv

from dna_sequence_analyzer.translate import TABLES, transcribe, translate_frames
from dna_sequence_analyzer.orfs import export_orfs, find_orfs
from dna_sequence_analyzer.remote_blast import submit, retrieve
from dna_sequence_analyzer.gc import GC_HIGH_PERCENT, GC_LOW_PERCENT
from dna_sequence_analyzer.tm import OLIGO_GUIDANCE_LENGTH, MeltingConditions, melting_temperature
from dna_sequence_analyzer.web_input import (
    bundled_example_fasta,
    describe_records,
    explain_tm_error,
    parse_sequence_text,
    record_labels,
    stats_csv,
    text_from_upload,
)
from dna_sequence_analyzer.windows import (
    gc_windows_for_records,
    skew_extrema,
    skew_windows_for_records,
)

_BLUE = "#6654ce"
_AMBER = "#aa641e"
_INDIGO = "#a13776"
_STYLE = '''<style>
.block-container {max-width:none;padding:3.4rem 1rem 1.1rem}
[data-testid="stAppViewContainer"] > section {padding-top:0}
[data-testid="stHeader"] {height:2.4rem;background:#252a32;color:#fff}
[data-testid="stHeader"] button {color:#fff}
.app-bar {position:fixed;z-index:1000000;inset:0 0 auto;height:2.4rem;padding-inline:3.45rem 1rem;
  display:flex;align-items:center;gap:1.1rem;
  background:#252a32;color:#fff;pointer-events:none}
.app-bar strong {font-size:.91rem}
.app-bar span {border-left:1px solid #667180;padding-left:1rem;font-size:.8rem;color:#dde4ef}
[data-testid="stSidebar"] {width:235px;min-width:235px;max-width:235px;border-right:1px solid #cbd3e0}
[data-testid="stSidebarHeader"] {background:#252a32}
[data-testid="stSidebarHeader"] button {color:#fff}
[data-testid="stSidebar"] [data-testid="stSidebarUserContent"] {padding:.8rem .85rem 1rem}
[data-testid="stSidebar"] h2 {font-size:1.05rem;padding:0}
[data-testid="stSidebar"] [role="radiogroup"] {gap:.12rem}
[data-testid="stSidebar"] [role="radiogroup"] label {padding:.3rem .5rem;border-radius:4px}
[data-testid="stSidebar"] [role="radiogroup"] label:has(input:checked) {background:#b9b9f155}
[data-testid="stSidebar"] [role="radiogroup"] label:has(input:checked) p {font-weight:650}
.st-key-tool_panels [role="tablist"] {display:none}
.st-key-sample_toolbar {padding:.53rem .7rem;border:1px solid #cbd3e0;border-radius:6px;background:#848fac0c}
.st-key-sample_toolbar [data-testid="stHorizontalBlock"] {align-items:center}
.st-key-sample_toolbar [data-testid="stPopoverButton"] {min-height:2rem;padding:.18rem .65rem}
.st-key-sample_toolbar [data-testid="stColumn"]:last-child [data-testid="stPopoverButton"] {
  background:#3032a8;color:#fff;border-color:#3032a8}
.sample-mast {display:flex;flex-wrap:wrap;align-items:center;gap:.3rem 1rem}
.sample-mast h1 {font-size:1.25rem;padding:0;overflow-wrap:anywhere}
.sample-mast p {font-size:.78rem;margin:0;font-variant-numeric:tabular-nums}
.sample-mast p span {display:inline-block;border-left:1px solid #cbd3e0;padding:0 .65rem}
.st-key-sequence_instrument {background:#20252d;color:#f1f4fa;border:1px solid #343c48;
  border-radius:6px;padding:.7rem .85rem .2rem}
.st-key-sequence_plot {overflow-x:auto;min-width:0}
.st-key-sequence_plot [data-testid="stPlotlyChart"] {min-width:650px}
.st-key-sequence_instrument [role="tablist"] {display:flex;border-bottom:1px solid #626a78}
.st-key-sequence_instrument [role="tab"] {color:#e2e8f4;min-height:2.15rem}
.st-key-sequence_instrument [role="tab"][aria-selected="true"] {color:#fff;background:#3a4557}
.sequence-view h2 {color:#f3f6fc;font:600 1.25rem/1.2 Consolas,monospace;
  padding:0;overflow-wrap:anywhere}
.sequence-view p {color:#bec9d9;font-size:.82rem;margin:.35rem 0 0;overflow-wrap:anywhere}
.st-key-record_inspector {min-height:29.3rem}
.st-key-record_inspector hr {margin:.25rem 0}
.st-key-record_inspector [data-testid="stForm"] {padding:.45rem}
.st-key-record_inspector h3 {font-size:1rem;padding:.1rem 0 .35rem}
.st-key-record_inspector [data-testid="stMarkdownContainer"] p {font-size:.86rem}
.st-key-viewer_registry h3 {font-size:1rem;padding:.1rem 0 .35rem}
.st-key-composition_metrics [data-testid="stMetric"] {
  padding:.75rem;border-radius:9px;background:rgba(126,101,199,.10)}
[data-testid="stMetricValue"] {font-variant-numeric:tabular-nums;font-size:1.6rem}
[role="tablist"] {padding:.35rem;gap:.3rem;border:1px solid var(--st-border-color);
  border-radius:9px;background:var(--st-secondary-background-color)}
[role="tab"] {padding-inline:.85rem;min-height:42px;border-radius:6px;font-weight:600}
[role="tab"][aria-selected="true"] {background:var(--st-background-color);color:var(--st-primary-color)}
[data-testid="stTabs"] [data-baseweb="tab-highlight"] {display:none}
code.protein {overflow-wrap:anywhere}
@media(max-width:760px) {
  .block-container {padding-top:3.1rem}
  .app-bar {background:linear-gradient(90deg,transparent 3rem,#252a32 3rem)}
  .app-bar span {display:none}
  .st-key-record_inspector {min-height:auto}
  .st-key-composition_metrics [data-testid="stHorizontalBlock"] {
    display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:.6rem}
  .st-key-composition_metrics [data-testid="stColumn"] {width:100%;min-width:0}
}
</style>'''


_COLLECTION = "Colección INS y GFP"


_VIEWER_SPAN = 100
_BASE_COLORS = {"A": "#f6c253", "C": "#67d9cd", "G": "#6ec7f4",
                "T": "#f5a677", "U": "#f5a677"}


def _region_figure(record, start, end):
    """Plot bases and GC windows against the same inclusive coordinates."""
    sequence = str(record.seq[start - 1:end]).upper()
    positions = list(range(start, end + 1))
    base_size = min(16, max(7, 575 // len(positions)))
    window = min(20, len(sequence))
    rows = gc_windows_for_records([record[start - 1:end]], window, min(4, window))
    figure = make_subplots(rows=2, cols=1, shared_xaxes=True,
                           row_heights=[.27, .73], vertical_spacing=.12)
    figure.add_trace(go.Scatter(
        x=positions, y=[0] * len(positions), mode="markers+text",
        marker=dict(symbol="square", size=base_size,
                    color=[_BASE_COLORS.get(base, "#b6b9ca") for base in sequence]),
        text=list(sequence), textposition="middle center",
        textfont=dict(size=max(6, base_size - 2), color="#1c2129", family="Consolas, monospace"),
        hovertemplate="Posición %{x} · %{text}<extra></extra>",
    ), row=1, col=1)
    figure.add_trace(go.Scatter(
        x=[start - 1 + (row.start + row.end) / 2 for row in rows],
        y=[row.gc_percent for row in rows],
        customdata=[(start - 1 + row.start, start - 1 + row.end) for row in rows],
        mode="lines+markers", line=dict(color="#53dfeb", width=2), marker=dict(size=4),
        hovertemplate="%{customdata[0]}–%{customdata[1]} · GC %{y:.1f} %<extra></extra>",
    ), row=2, col=1)
    figure.update_xaxes(range=[start - .5, end + .5], gridcolor="#394350", zeroline=False)
    figure.update_xaxes(side="top", showticklabels=True, ticks="outside", ticklen=5,
                        showline=True, linecolor="#a0acba", row=1, col=1)
    figure.update_xaxes(title_text="Posición (base)", row=2, col=1)
    figure.update_yaxes(visible=False, fixedrange=True, row=1, col=1)
    figure.update_yaxes(title_text="GC (%)", range=[0, 100], dtick=25, gridcolor="#394350",
                        fixedrange=True, row=2, col=1)
    figure.update_layout(height=290, margin=dict(t=28, b=44, l=53, r=14),
                         paper_bgcolor="#20252d", plot_bgcolor="#20252d",
                         font=dict(color="#d9e2f0"), showlegend=False, dragmode="pan")
    return figure, rows


def _quality_panel(checks):
    issues = sum(level != "OK" for level, _ in checks)
    if not issues:
        st.badge("Sin incidencias", color="green")
        return
    label = "1 incidencia" if issues == 1 else f"{issues} incidencias"
    with st.expander(label, expanded=True):
        for level, message in checks:
            st.text(f"{level} · {message}")


def _select_viewer_row(labels):
    rows = st.session_state.viewer_rows.selection.rows
    if rows and rows[0] < len(labels):
        st.session_state.viewer_record = labels[rows[0]]


def _tab_viewer(record, item, stats, labels, checks, params, results):
    plot, inspector = st.columns([3.2, 1.25], gap="small")
    length = len(record.seq)
    region = (1, min(length, _VIEWER_SPAN))
    with inspector, st.container(key="record_inspector", border=True):
        st.markdown("### Control de calidad")
        _quality_panel(checks)
        st.markdown("---")
        st.markdown("### Registro seleccionado")
        gc = "—" if item.gc_percent is None else f"{item.gc_percent:.2f} %"
        st.text(f"Longitud: {length:,} bases\nGC: {gc}\nBases ambiguas: {item.ambiguous_bases}")
        if length:
            record_key = (st.session_state.viewer_record, length)
            if st.session_state.get("viewer_region_record") == record_key:
                region = st.session_state.viewer_region
            with st.form("viewer-region"):
                st.markdown("### Región (base 1)")
                first, last = st.columns(2)
                start = first.number_input("Desde", 1, length, region[0])
                end = last.number_input("Hasta", 1, length, region[1])
                if st.form_submit_button("Aplicar", type="primary", width="stretch"):
                    if end < start or end - start + 1 > _VIEWER_SPAN:
                        st.warning(f"Elige un tramo de 1 a {_VIEWER_SPAN} bases.")
                    else:
                        region = (start, end)
                        st.session_state.viewer_region_record = record_key
                        st.session_state.viewer_region = region
    with plot, st.container(key="sequence_instrument"):
        start, end = region
        detail = f"Bases {start}–{end} de {length:,}" if length else "Registro vacío"
        description = record.description.removeprefix(record.id).strip()
        if description:
            detail = f"{detail} · {description}"
        sequence_tab, windows_tab, annotation_tab = st.tabs(["Secuencia", "Vista por tramos", "Anotaciones"])
        with sequence_tab:
            st.html(f'<section class="sequence-view"><h2>{escape(record.id)}</h2>'
                    f'<p>{escape(detail)}</p></section>')
        if length:
            figure, rows = _region_figure(record, start, end)
            with sequence_tab, st.container(key="sequence_plot"):
                st.plotly_chart(figure, theme=None, config={"displaylogo": False}, width="stretch")
            params["Región del visor"] = f"{record.id}: {start}–{end}"
            results["GC de la región"] = [
                {**asdict(row), "start": start - 1 + row.start, "end": start - 1 + row.end}
                for row in rows
            ]
            with windows_tab:
                st.dataframe(results["GC de la región"], hide_index=True,
                             column_order=("start", "end", "gc_percent", "ambiguous_bases"),
                             column_config={"start": "Desde", "end": "Hasta",
                                            "gc_percent": "GC (%)", "ambiguous_bases": "Bases ambiguas"})
        with annotation_tab:
            annotation = escape(description or "El encabezado FASTA solo contiene el identificador.")
            st.html(f'<section class="sequence-view"><p>{annotation}</p></section>')
    with st.container(key="viewer_registry", border=True):
        st.markdown("### Registros")
        _records_table(
            _record_rows(stats), key="viewer_rows", selection_mode="single-row",
            on_select=lambda: _select_viewer_row(labels),
        )


def main() -> None:
    st.set_page_config(page_title="DNA Sequence Analyzer", layout="wide")
    st.html(_STYLE)
    st.html('<header class="app-bar"><strong>DNA Sequence Analyzer</strong>'
            '<span>Visor de secuencia</span></header>')
    with st.container(key="sample_toolbar"):
        sample, intake, export = st.columns([5, 1.35, 1.1], gap="small")
    with intake.popover("Importar FASTA", width="stretch"):
        records, error, filename = _read_records()
    records = records or []
    stats = describe_records(records)
    checks = quality(records, stats, error) if records or error else [
        ("Aviso", "Sin muestra cargada.")
    ]
    with sample:
        if records:
            defined = [item.gc_percent for item in stats if item.gc_percent is not None]
            mean_gc = f"{sum(defined) / len(defined):.1f} %" if defined else "—"
            total = sum(item.length for item in stats)
            st.html(f'<section class="sample-mast"><h1>{escape(filename)}</h1>'
                    f'<p><span>Registros <b>{len(records)}</b></span>'
                    f'<span>Bases totales <b>{total:,}</b></span>'
                    f'<span>GC medio <b>{mean_gc}</b></span></p></section>')
        else:
            st.title("Importar secuencias")
    params = {"Origen": st.session_state.source}
    results = {}
    if records:
        labels = record_labels(records)
        st.sidebar.header("Muestra")
        st.sidebar.text(filename)
        if st.session_state.get("viewer_record") not in labels:
            st.session_state.viewer_record = labels[0]
        if len(labels) <= 8:
            selected = st.sidebar.radio("Registro en el visor", labels, key="viewer_record", label_visibility="collapsed")
        else:
            selected = st.sidebar.selectbox("Registro en el visor", labels, key="viewer_record")
        index = labels.index(selected)
        tools = ["Resumen", "Visor", "Traducción", "ORFs", "BLAST", "GC por ventanas",
                 "Sesgo GC", "Tm", "Guía"]
        st.sidebar.divider()
        active_tool = st.sidebar.radio("Herramientas", tools, index=1, key="workspace_tool",
                                       label_visibility="collapsed")
        with st.container(key="tool_panels"):
            with st.container(key=f"panel_{tools.index(active_tool)}"):
                overview, viewer, translation, orfs, blast, windows, skew, tm, guide = st.tabs(
                    tools, default=active_tool,
                )
        with viewer:
            _tab_viewer(records[index], stats[index], stats, labels, checks, params, results)
        with overview:
            _tab_resumen(stats, results)
        with windows:
            _tab_ventanas(records, results, params)
        with skew:
            _tab_skew(records, results, params)
        with translation:
            _tab_traduccion(records, results, params)
        with orfs:
            _tab_orfs(records)
    else:
        if error:
            st.error(error)
        else:
            st.info("Abre Importar FASTA para cargar un archivo o pegar una secuencia.")
        with st.expander("Herramientas disponibles sin archivo"):
            tm, blast, guide = st.tabs(["Tm", "BLAST", "Guía"])
    with tm:
        _tab_tm(results, params)
    with blast:
        _tab_blast(records, results)
    with guide:
        _tab_aprende()
    metadata = run_metadata(records, params, st.session_state.get("metadata"))
    st.session_state.metadata = metadata
    with export.popover("Exportar", width="stretch"):
        st.download_button("Descargar informe HTML", report(metadata, checks, results),
                           f"analisis-{metadata['ID']}.html", "text/html")
        available = [title for title, rows in results.items() if rows]
        if available:
            selected = st.selectbox("Tabla de resultados", available)
            st.download_button("Descargar tabla CSV", result_csv(results[selected]),
                               f"{selected}-{metadata['ID']}.csv", "text/csv")
        with st.expander("Metadatos del análisis"):
            st.json(metadata, expanded=False)


def _tab_blast(records, results):
    st.subheader("Buscar secuencias similares")
    source = st.radio("Consulta BLAST", ["Accesión pública", "Secuencia de la muestra"], horizontal=True)
    if source == "Accesión pública":
        query = st.text_input("Accesión de nucleótidos", placeholder="NM_000207.3").strip()
        if query and not re.fullmatch(r"[A-Za-z]{1,6}_?[0-9]+(?:\.[0-9]+)?", query):
            st.info("Introduce un identificador de nucleótidos con su número de accesión.")
            return
        confirmed = True
    else:
        if not records:
            st.info("Carga una muestra para seleccionar un registro.")
            return
        record, _ = _choose_record(records, "blast")
        query = str(record.seq)
        confirmed = st.checkbox("Confirmo que esta secuencia es pública, no es material genómico humano y autorizo su envío a NCBI.")
    st.caption("BLASTn compara la consulta con core_nt en NCBI. Se solicitan hasta 20 coincidencias.")
    if not query:
        return
    searches = st.session_state.setdefault("blast_searches", {})
    try:
        if query not in searches:
            if not st.button("Enviar consulta a NCBI", disabled=not confirmed):
                return
            searches[query] = {"job": submit(query), "hits": None}
        search = searches[query]
        job = search["job"]
        st.caption(f"Identificador de búsqueda · {job.rid}")
        if search["hits"] is None:
            st.info("Consulta el resultado pasado un minuto. NCBI puede tardar varios minutos.")
            if st.button("Consultar resultado"):
                search["hits"] = retrieve(job)
        if search["hits"] is not None:
            if search["hits"]:
                rows = [asdict(hit) for hit in search["hits"]]
                st.dataframe(rows, hide_index=True)
                results["BLAST"] = rows
                st.caption("La identidad corresponde al segmento con menor E-value de cada coincidencia.")
            else:
                st.info("NCBI completó la búsqueda sin coincidencias.")
    except (OSError, ValueError) as exc:
        st.error(str(exc))


@st.fragment
def _tab_orfs(records) -> None:
    st.subheader("Marcos abiertos de lectura")
    if not records:
        st.info("Carga una muestra FASTA para buscar ORFs.")
        return
    record, _ = _choose_record(records, "orfs")
    with st.form("orf-settings"):
        table_col, minimum_col = st.columns(2)
        table = table_col.selectbox("Tabla genética NCBI", TABLES)
        minimum = minimum_col.number_input("Longitud mínima (aa)", min_value=1, value=30)
        alternative = st.checkbox("Incluir inicios alternativos")
        nested = st.checkbox("Incluir ORFs anidados", value=True)
        submitted = st.form_submit_button("Buscar ORFs", type="primary")
    if submitted:
        try:
            hits = find_orfs(str(record.seq), table=table, min_aa=minimum,
                             alternative_starts=alternative, include_nested=nested)
        except ValueError as exc:
            st.error(str(exc))
            return
        st.session_state.orf_result = (record.id, str(record.seq), hits, table, minimum, alternative, nested)
    stored = st.session_state.get("orf_result")
    if stored is None or stored[:2] != (record.id, str(record.seq)):
        return
    _, _, hits, table, minimum, alternative, nested = stored
    st.caption(f"Tabla {table} · Mínimo {minimum} aa · "
               f"Inicios {'de la tabla' if alternative else 'ATG'} · "
               f"Anidados {'incluidos' if nested else 'excluidos'}")
    if not hits:
        st.info("No hay ORFs completos para estos parámetros. Prueba una longitud mínima menor.")
        return
    st.markdown(f"**{len(hits):,} ORFs completos**")
    figure = go.Figure()
    for strand, color in (("+", _BLUE), ("-", _INDIGO)):
        x, y = [], []
        for hit in hits:
            if hit.frame.startswith(strand):
                x.extend((hit.start + 1, hit.end, None))
                y.extend((hit.frame, hit.frame, None))
        figure.add_scatter(x=x, y=y, mode="lines", name=f"Hebra {strand}",
                           line=dict(color=color, width=7), hoverinfo="skip")
    figure.update_layout(height=300, margin=dict(l=0, t=12, r=12, b=0),
                         xaxis=dict(title="Posición en la secuencia (nt)", range=[0, len(record)]),
                         yaxis=dict(title="Marco", type="category", categoryorder="array",
                                    categoryarray=["-3", "-2", "-1", "+3", "+2", "+1"]))
    st.plotly_chart(figure)
    table_rows = [{"Marco": hit.frame, "Inicio": hit.start + 1, "Fin": hit.end,
                   "Longitud (aa)": len(hit.protein), "Proteína": hit.protein} for hit in hits]
    st.dataframe(table_rows, hide_index=True)
    st.caption("Coordenadas 1-based inclusivas, con el codón de parada. "
               "Los ORFs son candidatos a regiones codificantes.")
    for column, kind, suffix in zip(st.columns(3), ("csv", "fasta", "gff"), ("csv", "faa", "gff3")):
        column.download_button(f"Descargar {suffix.upper()}",
                               export_orfs([(record.id, hit) for hit in hits], kind),
                               file_name=f"orfs.{suffix}", mime="text/plain")


def _tab_traduccion(records, results, params) -> None:
    if not records:
        st.info("Añade una secuencia para ver su traducción.")
        return
    with st.expander("Secuencia y código genético"):
        record, _ = _choose_record(records, "traduccion")
        table = st.selectbox("Tabla del código genético (NCBI)", TABLES)
    if record is None:
        return
    params.update({"Registro de traducción": st.session_state["pick-traduccion"], "Código genético": table})
    sequence = str(record.seq)
    try:
        frames = translate_frames(sequence, table)
    except ValueError as exc:
        results["Traducción"] = [{"Error": str(exc)}]
        st.error(str(exc))
        return
    results["Traducción"] = [asdict(frame) for frame in frames]
    results["Transcripción"] = [{"ARN": transcribe(sequence)}]
    with st.expander("ARN transcrito de la hebra de entrada"):
        st.code(transcribe(sequence), language=None)
    st.caption("El asterisco (*) indica una parada. X indica un aminoácido indeterminado.")
    for item in frames:
        st.markdown(f"**Marco {item.frame}**")
        if item.protein:
            protein = escape(item.protein).replace("*", "<mark>*</mark>")
            st.markdown(f'<code class="protein">{protein}</code>',
                        unsafe_allow_html=True)
        else:
            st.caption("Sin codones completos.")


def _read_records():
    """Read the selected input and return records or a parsing error."""
    source = st.selectbox(
        "Origen", ["Subir un FASTA", "Pegar texto", _COLLECTION], key="source",
    )
    try:
        if source == _COLLECTION:
            st.caption("Transcritos públicos de insulina humana y GFP.")
            path = bundled_example_fasta()
            return parse_sequence_text(path.read_text(encoding="utf-8")), None, _COLLECTION
        if source == "Pegar texto":
            pasted = st.text_area("FASTA o secuencia", height=130)
            return parse_sequence_text(pasted) if pasted.strip() else [], None, "Secuencia pegada"
        uploaded = st.file_uploader("Archivo FASTA", type=["fasta", "fa", "fna", "txt"])
        if uploaded is None:
            return [], None, ""
        return parse_sequence_text(text_from_upload(uploaded.getvalue())), None, uploaded.name
    except (ValueError, FileNotFoundError, OSError) as exc:
        return [], str(exc), ""


def _record_rows(stats):
    return [
        {
            "id": item.id,
            "longitud": item.length,
            "GC %": None if item.gc_percent is None else round(item.gc_percent, 2),
            "ambiguas": item.ambiguous_bases,
            "tipo": item.nucleic_acid,
            "aviso": item.gc_flag or "",
        }
        for item in stats
    ]


def _records_table(rows, **selection):
    st.dataframe(rows, row_height=34, width="stretch", hide_index=True, column_config={
        "id": "Identificador", "longitud": "Longitud (bases)", "ambiguas": "Bases ambiguas",
        "tipo": "Ácido nucleico", "aviso": "Aviso",
        "GC %": st.column_config.ProgressColumn("GC (%)", min_value=0, max_value=100, format="%.2f"),
    }, **selection)


def _tab_resumen(stats, results) -> None:
    st.header("Composición de la muestra")
    method = (
        "Cada fila es una secuencia. Las letras ambiguas, como N, se cuentan en su columna "
        "y no entran en el porcentaje. Por debajo del "
        f"{GC_LOW_PERCENT:.0f} % la marca es low-gc; por encima del "
        f"{GC_HIGH_PERCENT:.0f} %, high-gc."
    )
    if not stats:
        st.caption("Sin registros para analizar.")
        return
    defined = [item.gc_percent for item in stats if item.gc_percent is not None]
    flagged = sum(1 for item in stats if item.gc_flag)
    mean_gc = sum(defined) / len(defined) if defined else None
    with st.container(key="composition_metrics"):
        cards = st.columns(4)
    cards[0].metric("Secuencias", len(stats), border=True)
    cards[1].metric("Bases en total", f"{sum(item.length for item in stats):,}".replace(",", "."), border=True)
    cards[2].metric("GC medio", "—" if mean_gc is None else f"{mean_gc:.1f} %", border=True)
    cards[3].metric("Con aviso de GC", flagged, border=True)
    table = _record_rows(stats)
    results["Resumen GC"] = table
    registry, distribution = st.columns([3, 2], gap="medium")
    with registry, st.container(border=True):
        st.subheader("Registros")
        _records_table(table)
        with st.expander("Cálculo y umbrales de GC"):
            st.write(method)
        st.download_button(
            "Descargar CSV",
            data=stats_csv(stats),
            file_name="gc_resumen.csv",
            mime="text/csv",
        )
    with distribution, st.container(border=True):
        st.subheader("Distribución GC")
        chart_rows = [row for row in table if row["GC %"] is not None]
        if not chart_rows:
            st.info("Ninguna secuencia tiene un GC calculable, así que no hay barras.")
            return
        import plotly.express as px

        figure = px.bar(
            chart_rows,
            x="id",
            y="GC %",
            color="aviso",
            color_discrete_map={"": _BLUE, "low-gc": _AMBER, "high-gc": _INDIGO},
            labels={"id": "Secuencia", "aviso": "Aviso"},
        )
        for threshold in (GC_LOW_PERCENT, GC_HIGH_PERCENT):
            figure.add_hline(y=threshold, line_dash="dot", line_color="#758692",
                             annotation_text=f"Umbral {threshold:g} %", annotation_position="top left")
        figure.update_layout(yaxis_title="GC (%)", yaxis_range=[0, 100], xaxis_title="Secuencia",
                             showlegend=any(row["aviso"] for row in chart_rows),
                             bargap=.65, height=260, margin=dict(t=18, b=32))
        st.plotly_chart(figure, width="stretch")


def _tab_ventanas(records, results, params) -> None:
    st.header("GC por ventanas")
    st.write(
        "Una media de toda la secuencia esconde tramos ricos en GC o en AT. "
        "La ventana es el trozo que se mide. El paso es cuánto avanza la siguiente ventana."
    )
    if not records:
        st.info("Carga una secuencia en Entrada de datos.")
        return
    with st.expander("Secuencia, ventana y paso"):
        chosen, subset = _choose_record(records, "ventana")
        if chosen is None:
            return
        longest = len(chosen.seq)
        window, step = _window_controls(longest, "gc")
    params.update({"Registro · GC por ventanas": st.session_state["pick-ventana"], "Ventana · GC por ventanas (bases)": window, "Paso · GC por ventanas (bases)": step})
    rows = [row for row in gc_windows_for_records(subset, window, step) if row.id == chosen.id]
    results["GC por ventanas"] = [asdict(row) for row in rows]
    if not rows:
        st.warning("La secuencia es más corta que la ventana. Baja el tamaño de la ventana.")
        return
    if len(rows) > 3000:
        st.warning("Hay muchas ventanas. Sube el paso si el gráfico va lento.")
    import plotly.graph_objects as go

    gc_x = [row.start for row in rows if row.gc_percent is not None]
    gc_y = [row.gc_percent for row in rows if row.gc_percent is not None]
    cpg_x = [row.start for row in rows if row.cpg_oe is not None]
    cpg_y = [row.cpg_oe for row in rows if row.cpg_oe is not None]
    figure = go.Figure()
    figure.add_trace(go.Scatter(x=gc_x, y=gc_y, name="GC %", mode="lines+markers", line=dict(color=_BLUE)))
    figure.add_trace(
        go.Scatter(
            x=cpg_x,
            y=cpg_y,
            name="CpG O/E",
            mode="lines+markers",
            yaxis="y2",
            line=dict(color=_INDIGO),
        )
    )
    figure.update_layout(
        height=460,
        margin=dict(t=30),
        xaxis_title="Inicio de la ventana (base 1)",
        yaxis_title="GC %",
        yaxis2=dict(title="CpG observado/esperado", overlaying="y", side="right"),
        legend=dict(orientation="h"),
    )
    st.plotly_chart(figure, width="stretch")
    st.caption(
        "O/E = (número de dinucleótidos CG × longitud) / (número de C × número de G). "
        "Un cribado clásico de islas CpG mira tramos largos con GC > 50 % y O/E > 0,6. "
        "Esta gráfica no decide si una ventana es una isla."
    )


def _tab_skew(records, results, params) -> None:
    st.header("Sesgo GC")
    st.write(
        "El sesgo es (G − C) / (G + C). En muchos cromosomas bacterianos la suma acumulada "
        "tiene un mínimo cerca del origen de replicación y un máximo cerca del término. "
        "Son candidatos, no una demostración."
    )
    if not records:
        st.info("Carga una secuencia en Entrada de datos.")
        return
    with st.expander("Parámetros del sesgo GC"):
        chosen, subset = _choose_record(records, "sesgo")
        if chosen is None:
            return
        longest = len(chosen.seq)
        window, step = _window_controls(longest, "skew")
    params.update({"Registro · GC skew": st.session_state["pick-sesgo"], "Ventana · GC skew (bases)": window, "Paso · GC skew (bases)": step})
    rows = [row for row in skew_windows_for_records(subset, window, step) if row.id == chosen.id]
    results["GC skew"] = [asdict(row) for row in rows]
    if not rows:
        st.warning("La secuencia es más corta que la ventana. Baja el tamaño de la ventana.")
        return
    lowest, highest = skew_extrema(rows)
    if lowest is not None and highest is not None:
        left, right = st.columns(2)
        left.metric(
            "Mínimo acumulado (posible origen)",
            f"{lowest.cumulative_skew:.3f}",
            delta=f"base {lowest.start}",
            delta_color="off",
        )
        right.metric(
            "Máximo acumulado (posible término)",
            f"{highest.cumulative_skew:.3f}",
            delta=f"base {highest.start}",
            delta_color="off",
        )
    import plotly.graph_objects as go

    defined = [row for row in rows if row.gc_skew is not None]
    cumulative = [row for row in rows if row.cumulative_skew is not None]
    skew_figure = go.Figure()
    skew_figure.add_trace(
        go.Scatter(
            x=[row.start for row in defined],
            y=[row.gc_skew for row in defined],
            name="Sesgo de la ventana",
            mode="lines+markers",
            line=dict(color=_BLUE),
        )
    )
    skew_figure.update_layout(
        height=360,
        margin=dict(t=30),
        xaxis_title="Inicio de la ventana (base 1)",
        yaxis_title="(G − C) / (G + C)",
    )
    st.plotly_chart(skew_figure, width="stretch")

    cumulative_figure = go.Figure()
    cumulative_figure.add_trace(
        go.Scatter(
            x=[row.start for row in cumulative],
            y=[row.cumulative_skew for row in cumulative],
            name="Sesgo acumulado",
            mode="lines",
            line=dict(color=_BLUE),
        )
    )
    _mark_extremum(cumulative_figure, lowest, "Posible origen", _INDIGO)
    if highest is not lowest:
        _mark_extremum(cumulative_figure, highest, "Posible término", _AMBER)
    cumulative_figure.update_layout(
        height=360,
        margin=dict(t=30),
        xaxis_title="Inicio de la ventana (base 1)",
        yaxis_title="Suma acumulada",
        legend=dict(orientation="h"),
    )
    st.plotly_chart(cumulative_figure, width="stretch")


def _tab_tm(results, params) -> None:
    st.header("Temperatura de fusión de un oligo")
    st.write(
        "Esto es para un cebador corto, no para un ARNm. El cálculo es el de vecinos más próximos: "
        "suma la energía de cada dinucleótido, no solo el porcentaje de GC. "
        "Por eso dos oligos con el mismo GC pueden fundir a distinta temperatura."
    )
    oligo = st.text_input("Secuencia del oligo", value="CGTTCCAAAGATGTGGGCATGAGCTTAC")
    columns = st.columns(4)
    sodium = columns[0].number_input("Na+ (mM)", min_value=0.0, value=50.0, step=10.0)
    magnesium = columns[1].number_input("Mg2+ (mM)", min_value=0.0, value=0.0, step=0.5)
    dntps = columns[2].number_input("dNTPs (mM)", min_value=0.0, value=0.0, step=0.1)
    oligo_nm = columns[3].number_input("Oligo (nM, cada hebra)", min_value=0.1, value=25.0, step=5.0)
    params.update({"Oligo": oligo, "Na+ (mM)": sodium, "Mg2+ (mM)": magnesium,
                   "dNTPs (mM)": dntps, "Oligo por hebra (nM)": oligo_nm})
    cleaned = "".join(oligo.split())
    if len(cleaned) > OLIGO_GUIDANCE_LENGTH:
        st.warning(
            f"Este oligo tiene {len(cleaned)} bases. El modelo está pensado para unos "
            f"{OLIGO_GUIDANCE_LENGTH} nucleótidos o menos. El número se calcula igual, "
            "pero no lo uses como temperatura de PCR de un gen largo."
        )
    try:
        temperature = melting_temperature(
            oligo,
            MeltingConditions(
                na_mM=sodium,
                mg_mM=magnesium,
                dntps_mM=dntps,
                dnac1_nM=oligo_nm,
                dnac2_nM=oligo_nm,
            ),
        )
    except ValueError as exc:
        results["Tm de primers"] = [{"Error": explain_tm_error(exc)}]
        st.error(explain_tm_error(exc))
        return
    results["Tm de primers"] = [{"Tm (°C)": temperature}]
    st.metric("Tm", f"{temperature:.2f} °C")
    st.subheader("Por qué sale este número")
    st.markdown(
        """
La estabilidad del dúplex depende sobre todo del **apilamiento** entre pares vecinos.
Los pasos con G o C suelen apilarse con más fuerza, pero el orden concreto importa.
Los puentes de hidrógeno (tres en G-C, dos en A-T) aportan menos a la diferencia neta,
porque la base suelta también forma puentes con el agua.

Los fosfatos del esqueleto están cargados negativo, así que las hebras se repelen.
Na+ y Mg2+ apantallan esa repulsión y **suben** la Tm. Los dNTPs secuestran Mg2+ y
devuelven parte de esa repulsión, así que pueden **bajar** la Tm respecto al mismo
magnesio sin dNTPs. Más oligo en el tubo también favorece un poco el dúplex.

La U, si la escribes, se lee como T. Esto no pliega ARN ni puntúa pares G·U.
        """
    )


def _tab_aprende() -> None:
    st.header("Aprende")
    st.markdown(
        """
### Contenido en GC

Es la fracción de las bases contadas que son G o C. Resume la composición de una
secuencia. No mide, por sí solo, lo estable que es la doble hélice.

Las letras ambiguas, como N, no entran en el porcentaje: no sabemos si son G o C.
Sí cuentan en la longitud y en la columna de ambiguas. Una secuencia hecha solo de
N tiene el GC indefinido, no 0 %. Cero significaría "todo es A o T".

### Apilamiento y puentes de hidrógeno

En el dúplex, cada base está emparejada con la de enfrente y apilada sobre el par
de al lado. La estabilidad depende sobre todo de ese **apilamiento**. Los pasos que
contienen G o C suelen ser más estables, y por eso un ADN rico en GC tiende a
separar sus hebras a más temperatura. El orden importa: `GCGCGCGC` y `GGGGCCCC`
tienen el mismo GC y no la misma Tm.

Los puentes de hidrógeno existen (tres entre G y C, dos entre A y T). Aportan menos
a la diferencia neta de lo que sugiere ese recuento, porque al abrir la hélice la
base sigue formando puentes con el agua. Por eso la Tm se calcula sumando cada
dinucleótido, con un modelo de vecinos más próximos, y no con una fórmula que solo
mire el porcentaje de GC.

### Islas CpG

CpG es el dinucleótido C seguida de G. En muchos genomas de vertebrados ese par
es más raro de lo que predice el azar. El cociente observado/esperado compara los
CG que hay con los que tocarían si C y G se colocaran al azar. Un cribado clásico
busca tramos de más de 200 bp, con GC por encima del 50 % y observado/esperado por
encima de 0,6. La pestaña de ventanas muestra los números. No etiqueta islas sola.

### Sesgo GC y replicación bacteriana

El sesgo (G − C) / (G + C) dice si una ventana tiene más G que C. En muchas
bacterias la hebra líder de la replicación es más rica en G. La suma acumulada
baja a lo largo de un replichoro y sube a lo largo del otro. El mínimo suele caer
cerca del origen (oriC) y el máximo cerca del término (ter). Es una pista para
mirar el cromosoma, no una prueba, y no tiene sentido en un cebador de 20 bases.

### Sesgo técnico en PCR e Illumina

Un molde muy rico en GC cuesta más de abrir en una PCR y a menudo amplifica peor.
Uno muy rico en AT se abre con facilidad y puede amplificar con otra eficiencia.
En Illumina la cobertura suele bajar en los dos extremos. El aviso `low-gc` (por
debajo del 25 %) y `high-gc` (por encima del 65 %) es un corte práctico de este
programa, no un umbral del fabricante. Sirve para mirar otra vez esa secuencia
antes de diseñar cebadores o de fiarse de una cobertura.

### ARN

Si hay uracilo y no hay timina, la secuencia se marca como ARN. La U cuenta como
una base que no es GC, igual que la T. Los tallos del ARN también se estabilizan
por apilamiento, y además el ARN forma pares G·U, llamados wobble, que no son el
par Watson-Crick pero aparecen en muchos tallos. Esta aplicación no pliega el ARN.
        """
    )


def _choose_record(records, key: str):
    """Return the selected record and a one-item list, or ``(None, None)``."""
    usable = [record for record in records if len(record.seq) > 0]
    if not usable:
        st.warning("Las secuencias están vacías.")
        return None, None
    labels = record_labels(usable)
    choice = st.selectbox("Secuencia", labels, key=f"pick-{key}")
    index = labels.index(choice)
    return usable[index], [usable[index]]


def _window_controls(longest: int, key: str) -> tuple[int, int]:
    if longest == 1:
        st.caption("Ventana y paso: 1 base.")
        return 1, 1
    upper = max(longest, 1)
    window_key = f"window-{key}"
    step_key = f"step-{key}"
    _keep_in_range(window_key, min(100, upper), upper)
    _keep_in_range(step_key, min(50, upper), upper)
    columns = st.columns(2)
    window = columns[0].slider("Tamaño de la ventana (bases)", min_value=1, max_value=upper, key=window_key)
    step = columns[1].slider("Paso (bases)", min_value=1, max_value=upper, key=step_key)
    return int(window), int(step)


def _keep_in_range(key: str, default: int, upper: int) -> None:
    """Stop a slider from keeping a value that no longer fits the sequence."""
    if key not in st.session_state:
        st.session_state[key] = default
        return
    st.session_state[key] = max(1, min(int(st.session_state[key]), upper))


def _mark_extremum(figure, row, name: str, color: str) -> None:
    if row is None or row.cumulative_skew is None:
        return
    figure.add_trace(
        {
            "type": "scatter",
            "x": [row.start],
            "y": [row.cumulative_skew],
            "mode": "markers+text",
            "name": name,
            "text": [name],
            "textposition": "top center",
            "marker": {"size": 12, "color": color},
        }
    )


if __name__ == "__main__":
    main()

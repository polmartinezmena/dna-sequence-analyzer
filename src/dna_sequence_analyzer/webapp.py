"""Streamlit interface for the DNA sequence analyzer.

Run it with ``dna-analyzer web`` or ``streamlit run app/streamlit_app.py``.
Every number comes from the same functions as the command line.
"""

from __future__ import annotations

import streamlit as st
from html import escape
from dataclasses import asdict
from hashlib import sha256
from os import getenv
from time import monotonic
from urllib.parse import quote
from dna_sequence_analyzer.blast import BlastOptions, search, session_wait_seconds
from dna_sequence_analyzer.lab import quality, run_metadata, report, result_csv

from dna_sequence_analyzer.plots import lab_template
from dna_sequence_analyzer.ui_theme import ACCENT, SERIES_2, format_es, stylesheet
from dna_sequence_analyzer.viewer import render_sequence_html
from dna_sequence_analyzer.translate import TABLES, transcribe, translate_frames
from dna_sequence_analyzer.gc import GC_HIGH_PERCENT, GC_LOW_PERCENT
from dna_sequence_analyzer.orfs import export_orfs, find_orfs
from dna_sequence_analyzer.tm import OLIGO_GUIDANCE_LENGTH, MeltingConditions, melting_temperature
from dna_sequence_analyzer.web_input import (
    bundled_example_fasta,
    describe_records,
    explain_tm_error,
    parse_sequence_text,
    looks_like_fasta,
    read_fasta_text,
    record_labels,
    text_from_upload,
)
from dna_sequence_analyzer.windows import (
    gc_windows_for_records,
    skew_extrema,
    skew_windows_for_records,
)

_MODULES = ("Resumen", "Secuencia y traducción", "ORFs", "BLAST", "Exportar", "Aprende")
_TYPE = {"dna": "ADN", "rna": "ARN", "mixed": "Mixto", "unspecified": "Indeterminado"}


def _header(record, item, checks):
    title = escape(record.id) if record else "Analizador de secuencias"
    description = escape(record.description.removeprefix(record.id).strip()) if record else ""
    details = "" if item is None else " · ".join((
        format_es(item.length, unit="bases"), _TYPE.get(item.nucleic_acid, item.nucleic_acid),
        f"GC {format_es(item.gc_percent, 2, '%')}",
        f"Ambiguas {format_es(item.ambiguous_bases)}",
    ))
    heading, qc = st.columns([5, 1])
    with heading:
        description_html = f"<p>{description}</p>" if description else ""
        st.html(f'<header class="lab-record"><h1>{title}</h1>{description_html}<p>{details}</p></header>')
    with qc:
        severity = next((level for level in ("Error", "Aviso") if any(
            state == level for state, _ in checks)), "OK")
        icon = {"OK": ":material/check_circle:", "Aviso": ":material/warning:", "Error": ":material/error:"}
        with st.popover(severity, icon=icon[severity], help="Control de calidad de la muestra"):
            for level, message in checks:
                st.write(f"**{level}** · {message}")


def main() -> None:
    st.set_page_config(page_title="Analizador de secuencias", layout="wide")
    mode = st.context.theme.type or "light"
    st.html(stylesheet(mode))
    template = lab_template(mode)
    records, error = _sidebar_records()
    records = records or []
    stats = describe_records(records)
    checks = quality(records, stats, error, comments=st.session_state.fasta_comments) if records or error else [("Aviso", "Sin muestra cargada.")]
    record = item = None
    if records:
        labels = record_labels(records)
        if st.session_state.get("registro_activo") not in labels:
            st.session_state.registro_activo = labels[0]
        descriptions = {
            label: f"{label} · {format_es(entry.length, unit='bases')} · "
                   f"{format_es(entry.gc_percent, 2, '%')}"
            for label, entry in zip(labels, stats)
        }
        choice = st.sidebar.selectbox("Registro activo", labels, key="registro_activo",
                                     format_func=descriptions.__getitem__)
        index = labels.index(choice)
        record, item = records[index], stats[index]
    else:
        st.session_state.pop("registro_activo", None)
    _header(record, item, checks)
    if error:
        st.error(error)
    elif not records:
        st.info("Carga un FASTA, pega una secuencia o selecciona Ejemplo en Muestra.")
    params = {"Origen": st.session_state.source}
    if record:
        params["Registro activo"] = st.session_state.registro_activo
    results = {}
    tabs = dict(zip(_MODULES, st.tabs(_MODULES, key="modulo_activo", on_change="rerun")))
    active = st.session_state.get("modulo_activo", _MODULES[0])
    sidebar_parameters = st.sidebar.container()
    sidebar_parameters.header("Parámetros")
    if active in ("Exportar", "Aprende") or (not record and active != "Oligo Tm"):
        sidebar_parameters.caption("Sin parámetros para esta vista.")
    with tabs["Resumen"]:
        _tab_resumen(stats, results)
        if record and len(record.seq):
            settings = sidebar_parameters if active == "Resumen" else st.container()
            with settings:
                if active == "Resumen":
                    window, step = _window_controls(len(record.seq), "composition")
                else:
                    upper = len(record.seq)
                    window = min(st.session_state.get("saved_window", 100), upper)
                    step = min(st.session_state.get("saved_step", 50), upper)
            _tab_ventanas(record, results, params, window, step, template)
            _tab_skew(record, results, params, window, step, template)
        _tab_tm(results, params)
    with tabs["Secuencia y traducción"]:
        _tab_traduccion(record, results, params, sidebar_parameters, active)
    with tabs["ORFs"]:
        _tab_orfs(record, results, params, sidebar_parameters, active, template)
    with tabs["BLAST"]:
        _tab_blast(record, results, params, sidebar_parameters, active)
    with tabs["Aprende"]:
        _tab_aprende()
    metadata = run_metadata(records, params, st.session_state.get("metadata"))
    st.session_state.metadata = metadata
    with tabs["Exportar"]:
        st.caption("Descarga el informe y las tablas desde el panel lateral.")
        with st.popover("Metadatos"):
            st.json(metadata, expanded=False)
    _exports(metadata, checks, results)


def _exports(metadata, checks, results):
    st.sidebar.header("Exportar")
    document = report(metadata, checks, results)
    document = document.replace("<title>", stylesheet() + "<title>", 1)
    document = document.replace("<h1>", '<main class="lab-report"><h1>', 1) + "</main>"
    run = metadata["ID"]
    st.sidebar.download_button("Informe HTML", document, f"analisis-{run}.html",
                               "text/html", type="primary")
    available = [title for title, rows in results.items() if rows]
    if available:
        with st.sidebar.popover("CSV de resultados"):
            selected = st.selectbox("Tabla", available)
            st.download_button("Descargar CSV", result_csv(results[selected]),
                               f"{selected}-{run}.csv", "text/csv")
    if results.get("ORFs"):
        with st.sidebar.popover("Archivos de ORFs"):
            for kind, extension in (("fasta", "fasta"), ("gff", "gff3")):
                hits = st.session_state.orf_result[2]
                st.download_button(f"Descargar {extension.upper()}", export_orfs(hits, kind),
                                   f"orfs-{run}.{extension}", "text/plain")
    st.sidebar.html(f'<div class="lab-run">ID de ejecución {escape(run[:12])}<br>'
                    f'SHA-256 {escape(metadata["SHA-256 de la muestra"][:12])}</div>')


def _tab_traduccion(record, results, params, settings, active):
    if record is None or not len(record.seq):
        st.info("Selecciona un registro con bases para ver su traducción.")
        return
    table = st.session_state.get("saved_genetic_table", 1)
    if active == "Secuencia y traducción":
        table = settings.selectbox("Tabla del código genético (NCBI)", TABLES,
                                   index=TABLES.index(table))
        st.session_state.saved_genetic_table = table
    params["Código genético"] = table
    sequence = str(record.seq)
    try:
        frames = translate_frames(sequence, table)
    except ValueError as exc:
        results["Traducción"] = [{"Error": str(exc)}]
        st.error(str(exc))
        st.html(render_sequence_html(sequence, 1, min(len(sequence), 60), []))
        return
    results["Traducción"] = [asdict(frame) for frame in frames]
    results["Transcripción"] = [{"ARN": transcribe(sequence)}]
    st.caption("Cada aminoácido aparece centrado sobre su codón. * indica una parada.")
    first, second, third = st.columns([2, 1, 1])
    start = int(first.number_input("Ir a posición", 1, len(sequence), 1))
    complement = second.checkbox("Hebra complementaria", value=True)
    rna = third.checkbox("Mostrar ARN")
    end = min(start + 119, len(sequence))
    st.html(render_sequence_html(sequence, start, end, frames, table, complement, rna))
    frame_name = st.selectbox("Marco de proteína", [frame.frame for frame in frames])
    protein = next(frame.protein for frame in frames if frame.frame == frame_name)
    if protein:
        st.subheader(f"Proteína del marco {frame_name}")
        lines = []
        for offset in range(0, len(protein), 60):
            fragment = protein[offset:offset + 60]
            groups = " ".join(fragment[i:i + 10] for i in range(0, len(fragment), 10))
            lines.append(f"{format_es(offset + 1):>8}  {groups}")
        st.html(f'<pre class="lab-protein">{escape(chr(10).join(lines))}</pre>')
        identifier = record.id.replace(chr(10), "").replace(chr(13), "")
        fasta = f">{identifier} marco={frame_name} tabla={table}\n" + "\n".join(
            protein[i:i + 60] for i in range(0, len(protein), 60)
        ) + "\n"
        st.download_button("Descargar proteína FASTA", fasta, "proteina.fasta", "text/plain")
    else:
        st.caption("El marco no tiene codones completos.")


def _orf_key(record):
    return record.id, sha256(str(record.seq).encode()).hexdigest()


def _tab_orfs(record, results, params, settings, active, template):
    if record is None or not len(record.seq):
        st.info("Selecciona un registro con bases para buscar ORFs.")
        return
    st.caption("Localiza tramos candidatos a codificar proteínas en los seis marcos.")
    key = _orf_key(record)
    saved = st.session_state.get("orf_result")
    if active == "ORFs":
        with settings:
            min_aa = st.number_input("Longitud mínima (aa)", min_value=1, value=30)
            table = st.selectbox("Código genético de ORFs (NCBI)", TABLES)
            nested = st.checkbox("Incluir inicios anidados")
            alternative = st.checkbox("Inicios alternativos de la tabla NCBI")
            terminal = st.checkbox("Incluir ORFs sin parada")
        if st.button("Buscar ORFs", type="primary"):
            try:
                hits = find_orfs(str(record.seq), table=table, min_aa=min_aa,
                                 include_nested=nested, allow_no_stop=terminal,
                                 alternative_starts=alternative)
            except ValueError as error:
                st.error(str(error))
                return
            parameters = {"Tabla de ORFs": table, "Longitud mínima de ORF (aa)": min_aa,
                          "Inicios alternativos": alternative, "Inicios anidados": nested, "ORFs sin parada": terminal}
            saved = (key, parameters, [(record.id, hit) for hit in hits])
            st.session_state.orf_result = saved
    if saved is None or saved[0] != key:
        st.info("Ajusta los parámetros y pulsa Buscar ORFs.")
        return
    params.update(saved[1])
    hits = [hit for _, hit in saved[2]]
    if not hits:
        st.info("No se encontraron ORFs con estos parámetros.")
        return
    results["ORFs"] = [
        {"id": record.id, **asdict(hit), "longitud": len(hit.protein),
         "Hebra": hit.strand.replace("-", "−"), "Posición del codón de inicio": hit.start_codon_position}
        for hit in hits
    ]
    st.dataframe(results["ORFs"], hide_index=True, width="stretch",
                 column_config={"id": "Registro", "frame": "Marco", "start": "Inicio",
                                "end": "Fin", "protein": "Proteína", "longitud": "Longitud (aa)",
                                "has_stop": "Parada", "start_codon": "Codón de inicio"})
    import plotly.graph_objects as go

    figure = go.Figure()
    for strand, color in (("+", ACCENT), ("-", SERIES_2)):
        coordinates, frames = [], []
        for hit in hits:
            if hit.frame.startswith(strand):
                coordinates.extend((hit.start, hit.end, None))
                frames.extend((hit.frame, hit.frame, None))
        figure.add_scatter(x=coordinates, y=frames, mode="lines+markers",
                           name=f"Hebra {strand}", line={"color": color, "width": 5})
    figure.update_layout(template=template, xaxis_title="Posición (base 1)",
                         yaxis={"categoryorder": "array",
                                "categoryarray": ["-3", "-2", "-1", "+3", "+2", "+1"]},
                         hovermode="closest")
    st.plotly_chart(figure, theme=None, width="stretch")


def _tab_blast(record, results, params, settings, active):
    if record is None or not len(record.seq):
        st.info("Selecciona un registro público con bases para consultar NCBI.")
        return
    key = _orf_key(record)
    saved = st.session_state.get("blast_result")
    if active == "BLAST":
        st.caption("Compara el registro activo o una proteína de sus ORFs con secuencias de NCBI.")
        with settings:
            program = st.segmented_control("Programa", ["blastn", "blastp"], default="blastn") or "blastn"
            databases = ["nt", "refseq_rna"] if program == "blastn" else ["nr"]
            database = st.selectbox("Base de datos", databases)
            evalue = st.number_input("Valor E máximo", min_value=0.000000001, value=0.01, format="%.9g")
            max_hits = st.number_input("Máximo de coincidencias", min_value=1, value=20)
            megablast = program == "blastn" and st.checkbox("MegaBLAST", value=True)
            timeout = st.number_input("Tiempo máximo (s)", min_value=1, value=300)
        query = str(record.seq)
        label = record.id
        if program == "blastp":
            orfs = st.session_state.get("orf_result")
            if not orfs or orfs[0] != key or not orfs[2]:
                st.info("Busca primero los ORFs del registro activo.")
                return
            proteins = {
                f"{hit.frame} · {format_es(hit.start)}–{format_es(hit.end)} · {format_es(len(hit.protein))} aa":
                hit.protein for _, hit in orfs[2]
            }
            label = st.selectbox("Proteína de un ORF", list(proteins))
            query = proteins[label]
        email = getenv("NCBI_EMAIL")
        if not email:
            st.info("Configura NCBI_EMAIL antes de iniciar el servidor para activar BLAST remoto.")
        confirmed = st.checkbox("La consulta es pública y autorizo enviarla a NCBI.")
        maximum, units = (10_000, "nt") if program == "blastn" else (2_000, "aa")
        oversized = len("".join(query.split())) > maximum
        if oversized:
            st.error(f"La consulta supera el máximo de {format_es(maximum, unit=units)} para {program}.")
        wait = session_wait_seconds(st.session_state.get("last_blast_search"), monotonic())
        if wait:
            st.info(f"Faltan {format_es(wait, unit='s')} para la próxima búsqueda BLAST.")
            st.button("Actualizar espera")
        if st.button("Lanzar BLAST", type="primary", disabled=not confirmed or not email or oversized or bool(wait)):
            options = BlastOptions(program, database, evalue, max_hits, megablast)
            st.session_state.last_blast_search = monotonic()
            with st.status("Preparando consulta…", expanded=True) as status:
                try:
                    hits = search(query, options, email=email, timeout=timeout, progress=status.write)
                except (ValueError, OSError) as error:
                    status.update(label="La consulta BLAST falló", state="error")
                    st.error(str(error))
                    return
                status.update(label="Consulta completada", state="complete", expanded=False)
            saved = (key, {"Consulta BLAST": label, "Programa BLAST": program,
                           "Base de datos BLAST": database, "Valor E máximo": evalue,
                           "Máximo de coincidencias": max_hits, "MegaBLAST": megablast if program == "blastn" else None}, [asdict(hit) for hit in hits])
            st.session_state.blast_result = saved
    if saved is None or saved[0] != key:
        return
    _, parameters, rows = saved
    params.update(parameters)
    results["BLAST"] = rows
    if not rows:
        st.info("NCBI no devolvió coincidencias.")
        return
    section = "protein" if parameters["Programa BLAST"] == "blastp" else "nuccore"
    view = [
        {**row, "acceso_ncbi": f"https://www.ncbi.nlm.nih.gov/{section}/{quote(row['accession'], safe='')}"}
        for row in rows
    ]
    st.caption(f"Última consulta · {parameters['Programa BLAST']} · {parameters['Consulta BLAST']}")
    st.dataframe(view, hide_index=True, width="stretch",
                 column_order=("acceso_ncbi", "description", "identity", "coverage", "evalue", "bitscore"),
                 column_config={
                     "acceso_ncbi": st.column_config.LinkColumn("Accesión", display_text=r"/([^/]+)$"),
                     "description": "Descripción",
                     "identity": st.column_config.NumberColumn("Identidad (%)", format="localized"),
                     "coverage": st.column_config.NumberColumn("Cobertura (%)", format="localized"),
                     "evalue": st.column_config.NumberColumn("Valor E", format="scientific"),
                     "bitscore": st.column_config.NumberColumn("Puntuación de bits", format="localized"),
                 })


def _sidebar_records():
    """Read the sidebar and return ``(records, error_message)``."""
    st.session_state.fasta_comments = 0
    st.sidebar.header("Muestra")
    source = st.sidebar.segmented_control(
        "Origen", ["Ejemplo", "Pegar", "Subir FASTA"], default="Subir FASTA",
        key="source", label_visibility="collapsed",
    ) or "Subir FASTA"
    try:
        if source == "Ejemplo":
            path = bundled_example_fasta()
            st.sidebar.caption(f"Archivo: `{path.name}`")
            return _parse_sample(path.read_text(encoding="utf-8")), None
        if source == "Pegar":
            pasted = st.sidebar.text_area(
                "FASTA o secuencia",
                height=220,
                placeholder=">mi_gen\nATGCATGC\n\no solo ATGCATGC",
            )
            if not pasted.strip():
                return None, "Pega un FASTA o una secuencia de bases en el panel de la izquierda."
            return _parse_sample(pasted), None
        uploaded = st.sidebar.file_uploader("Archivo FASTA", type=["fasta", "fa", "fna", "txt"])
        if uploaded is None:
            return None, None
        return _parse_sample(text_from_upload(uploaded.getvalue())), None
    except (ValueError, FileNotFoundError, OSError) as exc:
        return None, str(exc)


def _parse_sample(text):
    if looks_like_fasta(text):
        records, st.session_state.fasta_comments = read_fasta_text(text)
        return records
    return parse_sequence_text(text)


def _tab_resumen(stats, results) -> None:
    st.header("Resumen de GC")
    st.caption("Compara los registros por longitud, contenido en GC y bases ambiguas.")
    if not stats:
        st.info("Cuando haya una secuencia válida, aquí verás la tabla de registros.")
        return
    defined = [item.gc_percent for item in stats if item.gc_percent is not None]
    flagged = sum(1 for item in stats if item.gc_flag)
    mean_gc = sum(defined) / len(defined) if defined else None
    cards = st.columns(4)
    cards[0].metric("Secuencias", len(stats))
    cards[1].metric("Bases en total", format_es(sum(item.length for item in stats)))
    cards[2].metric("GC medio", format_es(mean_gc, 1, "%"))
    cards[3].metric("Con aviso de GC", flagged)
    table = [
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
    results["Resumen GC"] = table
    st.dataframe(table, row_height=32, width="stretch", hide_index=True, column_config={
        "id": st.column_config.TextColumn("Identificador", help="Identificador FASTA"),
        "longitud": st.column_config.NumberColumn("Longitud (bases)", format="localized"),
        "GC %": st.column_config.ProgressColumn("GC (%)", min_value=0, max_value=100, format="localized"),
        "ambiguas": st.column_config.NumberColumn("Bases ambiguas", format="localized"),
        "tipo": "Tipo", "aviso": "Aviso",
    })


def _tab_ventanas(record, results, params, window, step, template) -> None:
    st.subheader("GC por ventanas y CpG")
    st.caption("El perfil compara el contenido de GC y el cociente CpG entre ventanas.")
    params.update({"Ventana (bases)": window, "Paso (bases)": step})
    rows = gc_windows_for_records([record], window, step)
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
    figure.add_trace(go.Scatter(x=gc_x, y=gc_y, name="GC %", mode="lines+markers", line=dict(color=ACCENT), fill="tozeroy", fillcolor="rgba(76,91,212,.14)"))
    figure.add_trace(
        go.Scatter(
            x=cpg_x,
            y=cpg_y,
            name="CpG O/E",
            mode="lines+markers",
            yaxis="y2",
            line=dict(color=SERIES_2, dash="dash"),
        )
    )
    figure.update_layout(
        template=template,
        margin=dict(t=30),
        xaxis_title="Inicio de la ventana (base 1)",
        yaxis_title="GC %",
        yaxis2=dict(title="CpG observado/esperado", overlaying="y", side="right"),
        legend=dict(orientation="h"),
    )
    for threshold in (GC_LOW_PERCENT, GC_HIGH_PERCENT):
        figure.add_hline(y=threshold, line_dash="dash", line_color="#8A94A6",
                         annotation_text=format_es(threshold, unit="%"))
    st.plotly_chart(figure, theme=None, width="stretch", config={"locale": "es"})


def _tab_skew(record, results, params, window, step, template) -> None:
    st.subheader("Sesgo GC")
    st.caption("Muestra el equilibrio entre G y C por ventana y su suma acumulada.")
    rows = skew_windows_for_records([record], window, step)
    results["GC skew"] = [asdict(row) for row in rows]
    if not rows:
        st.warning("La secuencia es más corta que la ventana. Baja el tamaño de la ventana.")
        return
    lowest, highest = skew_extrema(rows)
    if lowest is not None and highest is not None:
        left, right = st.columns(2)
        left.metric(
            "Mínimo acumulado (posible origen)",
            format_es(lowest.cumulative_skew, 3),
            delta=f"base {format_es(lowest.start)}",
            delta_color="off",
        )
        right.metric(
            "Máximo acumulado (posible término)",
            format_es(highest.cumulative_skew, 3),
            delta=f"base {format_es(highest.start)}",
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
            line=dict(color=ACCENT),
        )
    )
    skew_figure.update_layout(
        template=template,
        margin=dict(t=30),
        xaxis_title="Inicio de la ventana (base 1)",
        yaxis_title="(G − C) / (G + C)",
    )
    st.plotly_chart(skew_figure, theme=None, width="stretch", config={"locale": "es"})

    cumulative_figure = go.Figure()
    cumulative_figure.add_trace(
        go.Scatter(
            x=[row.start for row in cumulative],
            y=[row.cumulative_skew for row in cumulative],
            name="Sesgo acumulado",
            mode="lines",
            line=dict(color=ACCENT),
        )
    )
    _mark_extremum(cumulative_figure, lowest, "Posible origen", SERIES_2)
    if highest is not lowest:
        _mark_extremum(cumulative_figure, highest, "Posible término", "#B26B00")
    cumulative_figure.update_layout(
        template=template,
        margin=dict(t=30),
        xaxis_title="Inicio de la ventana (base 1)",
        yaxis_title="Suma acumulada",
        legend=dict(orientation="h"),
    )
    st.plotly_chart(cumulative_figure, theme=None, width="stretch", config={"locale": "es"})


def _tab_tm(results, params) -> None:
    st.header("Temperatura de fusión de un oligo")
    st.caption("Estima la temperatura de fusión de un oligo corto según su secuencia y concentraciones.")
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
    st.metric("Tm", format_es(temperature, 2, "°C"))


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


def _window_controls(longest: int, key: str) -> tuple[int, int]:
    if longest == 1:
        st.caption("Ventana y paso: 1 base.")
        return 1, 1
    upper = max(longest, 1)
    window_key = f"window-{key}"
    step_key = f"step-{key}"
    _keep_in_range(window_key, min(st.session_state.get("saved_window", 100), upper), upper)
    _keep_in_range(step_key, min(st.session_state.get("saved_step", 50), upper), upper)
    columns = st.columns(2)
    window = columns[0].slider("Tamaño de la ventana (bases)", min_value=1, max_value=upper, key=window_key)
    step = columns[1].slider("Paso (bases)", min_value=1, max_value=upper, key=step_key)
    st.session_state.saved_window = int(window)
    st.session_state.saved_step = int(step)
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

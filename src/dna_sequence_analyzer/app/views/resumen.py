"""Vista resumen."""

import streamlit as st
from dataclasses import asdict

from dna_sequence_analyzer.gc import GC_HIGH_PERCENT, GC_LOW_PERCENT
from dna_sequence_analyzer.plots import lab_template
from dna_sequence_analyzer.ui_theme import ACCENT, SERIES_2, format_es, spanish_table
from dna_sequence_analyzer.windows import gc_windows_for_records, skew_extrema, skew_windows_for_records


_SEQUENCE_TYPES = {"dna": "ADN", "rna": "ARN", "mixed": "Mixto", "unspecified": "Indeterminado"}
_GC_FLAGS = {"low-gc": "GC bajo", "high-gc": "GC alto"}


def render(record, item, stats, checks, results, params, settings, exports):
    with st.container(border=True, key="tarjeta-resumen"):
        _tab_resumen(stats, results)
    if not len(record.seq):
        return
    with settings:
        window, step = _window_controls(len(record.seq), "composition")
    template = lab_template("light")
    with st.container(border=True, key="tarjeta-ventanas"):
        _tab_ventanas(record, results, params, window, step, template)
    with st.container(border=True, key="tarjeta-sesgo"):
        _tab_skew(record, results, params, window, step, template)


def _tab_resumen(stats, results) -> None:
    st.header("Resumen de GC")
    st.caption("Compara los registros por longitud, contenido en GC y bases ambiguas.")
    if not stats:
        st.info("Cuando haya una secuencia válida, aquí verás la tabla de registros.")
        return
    table = [
        {
            "id": item.id,
            "longitud": item.length,
            "GC %": None if item.gc_percent is None else round(item.gc_percent, 2),
            "perfil": item.gc_percent,
            "ambiguas": item.ambiguous_bases,
            "tipo": _SEQUENCE_TYPES[item.nucleic_acid],
            "aviso": _GC_FLAGS.get(item.gc_flag, ""),
        }
        for item in stats
    ]
    results["Resumen GC"] = [{k: v for k, v in row.items() if k != "perfil"} for row in table]
    st.dataframe(spanish_table(table, {"GC %": 2}), row_height=32, width="stretch", hide_index=True, column_config={
        "id": st.column_config.TextColumn("Identificador", help="Identificador FASTA"),
        "longitud": "Longitud (bases)",
        "GC %": "GC (%)",
        "perfil": st.column_config.ProgressColumn("Perfil GC", min_value=0, max_value=100),
        "ambiguas": "Bases ambiguas",
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
    figure.add_trace(go.Scatter(x=gc_x, y=gc_y, name="GC %", mode="lines+markers", line=dict(color=ACCENT), fill="tozeroy", fillcolor="rgba(224,138,46,.14)"))
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
    st.plotly_chart(figure, theme=None, width="stretch")


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
    st.plotly_chart(skew_figure, theme=None, width="stretch")

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
    st.plotly_chart(cumulative_figure, theme=None, width="stretch")


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

"""Streamlit interface for the DNA sequence analyzer.

Run it with ``dna-analyzer web`` or ``streamlit run app/streamlit_app.py``.
Every number comes from the same functions as the command line.
"""

from __future__ import annotations

import streamlit as st

from dna_sequence_analyzer import __version__
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

_TEAL = "#0B6E6A"
_SAND = "#C4A574"
_CORAL = "#C65B4A"


def main() -> None:
    """Draw the whole page. Streamlit runs this again on every interaction."""
    st.set_page_config(
        page_title="Analizador de ADN",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded",
    )
    st.title("🧬 Analizador de secuencias de ADN")
    st.caption(
        f"Versión {__version__}. Contenido en GC, ventanas, sesgo y temperatura de fusión, "
        "sin usar la terminal."
    )
    records, error = _sidebar_records()
    if error:
        st.error(error)
    tabs = st.tabs(
        ["Resumen GC", "GC por ventanas", "GC skew", "Tm de primers", "Aprende"]
    )
    with tabs[0]:
        _tab_resumen(records)
    with tabs[1]:
        _tab_ventanas(records)
    with tabs[2]:
        _tab_skew(records)
    with tabs[3]:
        _tab_tm()
    with tabs[4]:
        _tab_aprende()


def _sidebar_records():
    """Read the sidebar and return ``(records, error_message)``."""
    st.sidebar.header("Secuencias")
    st.sidebar.caption(
        "Las tres primeras pestañas usan esta entrada. La pestaña de Tm tiene su propio oligo."
    )
    source = st.sidebar.radio(
        "Origen",
        ["Ejemplo del repositorio", "Pegar texto", "Subir un FASTA"],
        index=0,
    )
    try:
        if source == "Ejemplo del repositorio":
            path = bundled_example_fasta()
            st.sidebar.caption(f"Archivo: `{path.name}`")
            return parse_sequence_text(path.read_text(encoding="utf-8")), None
        if source == "Pegar texto":
            pasted = st.sidebar.text_area(
                "FASTA o secuencia",
                height=220,
                placeholder=">mi_gen\nATGCATGC\n\no solo ATGCATGC",
            )
            if not pasted.strip():
                return None, "Pega un FASTA o una secuencia de bases en el panel de la izquierda."
            return parse_sequence_text(pasted), None
        uploaded = st.sidebar.file_uploader("Archivo FASTA", type=["fasta", "fa", "fna", "txt"])
        if uploaded is None:
            return None, "Sube un archivo FASTA para empezar."
        return parse_sequence_text(text_from_upload(uploaded.getvalue())), None
    except (ValueError, FileNotFoundError, OSError) as exc:
        return None, str(exc)


def _tab_resumen(records) -> None:
    st.header("Resumen de GC")
    st.write(
        "Cada fila es una secuencia. Las letras ambiguas, como N, se cuentan en su columna "
        "y no entran en el porcentaje. Por debajo del "
        f"{GC_LOW_PERCENT:.0f} % la marca es low-gc; por encima del "
        f"{GC_HIGH_PERCENT:.0f} %, high-gc."
    )
    if not records:
        st.info("Cuando haya una secuencia válida, aquí verás la tabla y el gráfico.")
        return
    stats = describe_records(records)
    defined = [item.gc_percent for item in stats if item.gc_percent is not None]
    flagged = sum(1 for item in stats if item.gc_flag)
    mean_gc = sum(defined) / len(defined) if defined else None
    cards = st.columns(4)
    cards[0].metric("Secuencias", len(stats))
    cards[1].metric("Bases en total", f"{sum(item.length for item in stats):,}".replace(",", "."))
    cards[2].metric("GC medio", "—" if mean_gc is None else f"{mean_gc:.1f} %")
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
    st.dataframe(table, width="stretch", hide_index=True)
    st.download_button(
        "Descargar CSV",
        data=stats_csv(stats),
        file_name="gc_resumen.csv",
        mime="text/csv",
    )
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
        color_discrete_map={"": _TEAL, "low-gc": _SAND, "high-gc": _CORAL},
        labels={"id": "Secuencia", "aviso": "Aviso"},
    )
    figure.update_layout(yaxis_title="GC %", xaxis_title="", height=420, margin=dict(t=30))
    st.plotly_chart(figure, width="stretch")


def _tab_ventanas(records) -> None:
    st.header("GC por ventanas")
    st.write(
        "Una media de toda la secuencia esconde tramos ricos en GC o en AT. "
        "La ventana es el trozo que se mide. El paso es cuánto avanza la siguiente ventana."
    )
    if not records:
        st.info("Carga una secuencia en el panel de la izquierda.")
        return
    chosen, subset = _choose_record(records, "ventana")
    if chosen is None:
        return
    longest = len(chosen.seq)
    window, step = _window_controls(longest, "gc")
    rows = [row for row in gc_windows_for_records(subset, window, step) if row.id == chosen.id]
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
    figure.add_trace(go.Scatter(x=gc_x, y=gc_y, name="GC %", mode="lines+markers", line=dict(color=_TEAL)))
    figure.add_trace(
        go.Scatter(
            x=cpg_x,
            y=cpg_y,
            name="CpG O/E",
            mode="lines+markers",
            yaxis="y2",
            line=dict(color=_CORAL),
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


def _tab_skew(records) -> None:
    st.header("Sesgo GC")
    st.write(
        "El sesgo es (G − C) / (G + C). En muchos cromosomas bacterianos la suma acumulada "
        "tiene un mínimo cerca del origen de replicación y un máximo cerca del término. "
        "Son candidatos, no una demostración."
    )
    if not records:
        st.info("Carga una secuencia en el panel de la izquierda.")
        return
    chosen, subset = _choose_record(records, "sesgo")
    if chosen is None:
        return
    longest = len(chosen.seq)
    window, step = _window_controls(longest, "skew")
    rows = [row for row in skew_windows_for_records(subset, window, step) if row.id == chosen.id]
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
            line=dict(color=_TEAL),
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
            line=dict(color=_TEAL),
        )
    )
    _mark_extremum(cumulative_figure, lowest, "Posible origen", _CORAL)
    if highest is not lowest:
        _mark_extremum(cumulative_figure, highest, "Posible término", _SAND)
    cumulative_figure.update_layout(
        height=360,
        margin=dict(t=30),
        xaxis_title="Inicio de la ventana (base 1)",
        yaxis_title="Suma acumulada",
        legend=dict(orientation="h"),
    )
    st.plotly_chart(cumulative_figure, width="stretch")


def _tab_tm() -> None:
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
        st.error(explain_tm_error(exc))
        return
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

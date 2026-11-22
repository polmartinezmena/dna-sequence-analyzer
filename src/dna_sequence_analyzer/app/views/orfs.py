"""Vista orfs."""

import streamlit as st
from dataclasses import asdict
from hashlib import sha256

from dna_sequence_analyzer.orfs import find_orfs
from dna_sequence_analyzer.plots import lab_template
from dna_sequence_analyzer.translate import TABLES
from dna_sequence_analyzer.ui_theme import ACCENT, SERIES_2, spanish_table


def record_key(record):
    return record.id, sha256(str(record.seq).encode()).hexdigest()


def render(record, item, stats, checks, results, params, settings, exports):
    if record is None or not len(record.seq):
        st.info("Selecciona un registro con bases para buscar ORFs.")
        return
    st.caption("Localiza tramos candidatos a codificar proteínas en los seis marcos.")
    key = record_key(record)
    saved = st.session_state.get("orf_result")
    with settings:
        min_aa = st.number_input("Longitud mínima (aa)", min_value=1, value=30)
        table = st.selectbox("Código genético de ORFs (NCBI)", TABLES)
        nested = st.checkbox("Incluir inicios anidados")
        alternative = st.checkbox("Inicios alternativos de la tabla NCBI")
        terminal = st.checkbox("Incluir ORFs sin parada")
        search_clicked = st.button("Buscar ORFs", type="primary", key="buscar-orfs")
    if search_clicked:
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
    st.dataframe(spanish_table(results["ORFs"]), hide_index=True, width="stretch",
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
    figure.update_layout(template=lab_template("light"), xaxis_title="Posición (base 1)",
                         yaxis={"categoryorder": "array",
                                "categoryarray": ["-3", "-2", "-1", "+3", "+2", "+1"]},
                         hovermode="closest")
    st.plotly_chart(figure, theme=None, width="stretch")

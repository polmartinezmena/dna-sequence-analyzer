"""Vista traduccion."""

import streamlit as st
from dataclasses import asdict
from html import escape

from dna_sequence_analyzer.translate import TABLES, transcribe, translate_frames
from dna_sequence_analyzer.ui_theme import format_es
from dna_sequence_analyzer.viewer import render_sequence_html


def render(record, item, stats, checks, results, params, settings, exports):
    if record is None or not len(record.seq):
        st.info("Selecciona un registro con bases para ver su traducción.")
        return
    with settings:
        table = st.selectbox("Tabla del código genético (NCBI)", TABLES)
        complement = st.checkbox("Hebra complementaria", value=True)
        rna = st.checkbox("Mostrar ARN")
    params.update({"Código genético": table, "Hebra complementaria": complement, "Mostrar ARN": rna})
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
    start = st.number_input("Ir a posición", 1, len(sequence), 1)
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
        exports["Proteína FASTA"] = (fasta, "proteina.fasta", "text/plain")
    else:
        st.caption("El marco no tiene codones completos.")

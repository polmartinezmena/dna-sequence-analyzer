"""Vista visor."""

from html import escape

import streamlit as st
from dna_sequence_analyzer.viewer import render_sequence_html


def render(record, item, stats, checks, results, params, settings, exports):
    settings.html('<h2 class="section-label">Control de calidad</h2>')
    for level, message in checks:
        settings.html(f"<p>{escape(level)} · {escape(message)}</p>")
    if record.seq:
        st.html(render_sequence_html(str(record.seq), 1, min(100, item.length)))
    else:
        st.info("El registro está vacío.")

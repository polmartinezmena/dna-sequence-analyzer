"""Vista blast."""

import streamlit as st
from dataclasses import asdict
from os import getenv
from time import monotonic
from urllib.parse import quote

from dna_sequence_analyzer.blast import BlastOptions, search, session_wait_seconds
from dna_sequence_analyzer.ui_theme import format_es, spanish_table

from .orfs import record_key


def render(record, item, stats, checks, results, params, settings, exports):
    if record is None or not len(record.seq):
        st.info("Selecciona un registro público con bases para consultar NCBI.")
        return
    key = record_key(record)
    saved = st.session_state.get("blast_result")
    st.caption("Compara el registro activo o una proteína de sus ORFs con secuencias de NCBI.")
    with settings:
        program = st.segmented_control("Programa", ["blastn", "blastp"], default="blastn") or "blastn"
        databases = ["nt", "refseq_rna"] if program == "blastn" else ["nr"]
        database = st.selectbox("Base de datos", databases)
        evalue = st.number_input("Valor E máximo", min_value=0.000000001, value=0.01, format="%.9g")
        max_hits = st.number_input("Máximo de coincidencias", min_value=1, value=20)
        megablast = program == "blastn" and st.checkbox("MegaBLAST", value=True)
        timeout = st.number_input("Tiempo máximo (s)", min_value=1, max_value=300, value=300)
    query, label = str(record.seq), record.id
    if program == "blastp":
        orfs = st.session_state.get("orf_result")
        if not orfs or orfs[0] != key or not orfs[2]:
            settings.info("Busca primero los ORFs del registro activo.")
            query = None
        else:
            proteins = {
                f"{hit.frame} · {format_es(hit.start)}–{format_es(hit.end)} · {format_es(len(hit.protein))} aa":
                hit.protein for _, hit in orfs[2]
            }
            label = settings.selectbox("Proteína de un ORF", list(proteins))
            query = proteins[label]
    email = getenv("NCBI_EMAIL")
    if not email:
        settings.info("Configura NCBI_EMAIL antes de iniciar el servidor para activar BLAST remoto.")
    confirmed = settings.checkbox("La consulta es pública y autorizo enviarla a NCBI.")
    maximum, units = (10_000, "nt") if program == "blastn" else (2_000, "aa")
    oversized = query is not None and len("".join(query.split())) > maximum
    if oversized:
        settings.error(f"La consulta supera el máximo de {format_es(maximum, unit=units)} para {program}.")
    wait = session_wait_seconds(st.session_state.get("last_blast_search"), monotonic())
    if wait:
        settings.info(f"Faltan {format_es(wait, unit='s')} para la próxima búsqueda BLAST.")
        settings.button("Actualizar espera")
    can_launch = bool(query and confirmed and email and not oversized and not wait)
    if settings.button("Lanzar BLAST", key="lanzar-blast", type="primary", disabled=not can_launch) and can_launch:
        options = BlastOptions(program, database, evalue, max_hits, megablast)
        st.session_state.last_blast_search = monotonic()
        st.session_state.pop("blast_error", None)
        with st.status("Preparando consulta…", expanded=True) as status:
            try:
                hits = search(query, options, email=email, timeout=timeout, progress=status.write)
            except (ValueError, OSError) as error:
                status.update(label="La consulta BLAST falló", state="error")
                st.session_state.blast_error = key, str(error)
            else:
                status.update(label="Consulta completada", state="complete", expanded=False)
                parameters = {
                    "Consulta BLAST": label, "Programa BLAST": program, "Base de datos BLAST": database,
                    "Valor E máximo": evalue, "Máximo de coincidencias": max_hits, "Tiempo máximo (s)": timeout,
                    "MegaBLAST": megablast if program == "blastn" else None,
                }
                st.session_state.blast_result = key, parameters, [asdict(hit) for hit in hits]
        st.rerun()
    if error := st.session_state.get("blast_error"):
        if error[0] == key:
            st.error(error[1])
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
    st.dataframe(spanish_table(view, {"identity": 2, "coverage": 2, "bitscore": 1}), hide_index=True, width="stretch",
                 column_order=("acceso_ncbi", "description", "identity", "coverage", "evalue", "bitscore"),
                 column_config={
                     "acceso_ncbi": st.column_config.LinkColumn("Accesión", display_text=r"/([^/]+)$"),
                     "description": "Descripción",
                     "identity": "Identidad (%)",
                     "coverage": "Cobertura (%)",
                     "evalue": st.column_config.NumberColumn("Valor E", format="scientific"),
                     "bitscore": "Puntuación de bits",
                 })

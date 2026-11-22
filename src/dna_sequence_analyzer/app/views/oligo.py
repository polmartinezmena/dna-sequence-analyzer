"""Vista oligo."""

import streamlit as st
from dna_sequence_analyzer.tm import OLIGO_GUIDANCE_LENGTH, MeltingConditions, melting_temperature
from dna_sequence_analyzer.ui_theme import format_es
from dna_sequence_analyzer.web_input import explain_tm_error


def render(record, item, stats, checks, results, params, settings, exports):
    st.header("Temperatura de fusión de un oligo")
    st.caption("Estima la temperatura de fusión de un oligo corto según su secuencia y concentraciones.")
    with settings, st.form("form-oligo", border=False):
        oligo = st.text_input("Secuencia del oligo", value="CGTTCCAAAGATGTGGGCATGAGCTTAC")
        sodium = st.number_input("Na+ (mM)", min_value=0.0, value=50.0, step=10.0)
        magnesium = st.number_input("Mg²+ (mM)", min_value=0.0, value=0.0, step=0.5)
        dntps = st.number_input("dNTPs (mM)", min_value=0.0, value=0.0, step=0.1)
        oligo_nm = st.number_input("Oligo (nM, cada hebra)", min_value=0.1, value=25.0, step=5.0)
        st.form_submit_button("Calcular Tm", type="primary")
    params.update({"Oligo": oligo, "Na+ (mM)": sodium, "Mg2+ (mM)": magnesium,
                   "dNTPs (mM)": dntps, "Oligo por hebra (nM)": oligo_nm})
    cleaned = "".join(oligo.split())
    if len(cleaned) > OLIGO_GUIDANCE_LENGTH:
        st.warning(
            f"Este oligo tiene {format_es(len(cleaned))} bases. El modelo está pensado para unos "
            f"{format_es(OLIGO_GUIDANCE_LENGTH)} nucleótidos o menos. El número se calcula igual, "
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

import sys
import shutil
import subprocess
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

st.set_page_config(layout="wide")
st.title("INDRA Projection")

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
OUTPUT_DIR = BASE_DIR / "outputs"

DATA_DIR.mkdir(exist_ok=True)
OUTPUT_DIR.mkdir(exist_ok=True)

uploaded_file = st.file_uploader(
    "Compendium Excel-Datei hochladen",
    type=["xlsx"]
)

def run_script(script_name):
    return subprocess.run(
        [sys.executable, script_name],
        cwd=BASE_DIR,
        capture_output=True,
        text=True
    )

if uploaded_file is not None:
    input_path = DATA_DIR / "compendium.xlsx"

    with open(input_path, "wb") as f:
        f.write(uploaded_file.getbuffer())

    st.success("Datei hochgeladen.")

    if st.button("Diagramm erzeugen"):
        with st.spinner("A1 läuft: Daten berechnen und Cartesian Product erzeugen..."):
            result_a1 = run_script("A1.py")

        if result_a1.returncode != 0:
            st.error("Fehler in A1.py")
            st.code(result_a1.stderr)
            st.stop()

        st.success("A1 fertig.")

        with st.expander("A1 Ausgabe anzeigen"):
            st.code(result_a1.stdout)

        cartesian_output = OUTPUT_DIR / "CartesianProduct_constraints.xlsx"
        processed_output = OUTPUT_DIR / "compendium_processed.xlsx"

        if not cartesian_output.exists():
            st.error(f"Nicht gefunden: {cartesian_output}")
            st.stop()

        if not processed_output.exists():
            st.error(f"Nicht gefunden: {processed_output}")
            st.stop()

        shutil.copy(cartesian_output, DATA_DIR / "CartesianProduct_constraints.xlsx")
        shutil.copy(processed_output, DATA_DIR / "compendium_processed.xlsx")

        with st.spinner("A2 läuft: Diagramm erzeugen..."):
            result_a2 = run_script("A2.py")

        if result_a2.returncode != 0:
            st.error("Fehler in A2.py")
            st.code(result_a2.stderr)
            st.stop()

        st.success("Diagramm erzeugt.")

        with st.expander("A2 Ausgabe anzeigen"):
            st.code(result_a2.stdout)

        html_path = OUTPUT_DIR / "Metanumber_Plot_Ca_HCO3_Bands.html"

        if html_path.exists():
            html = html_path.read_text(encoding="utf-8")
            components.html(html, height=950, scrolling=True)

            with open(html_path, "rb") as f:
                st.download_button(
                    "HTML-Plot herunterladen",
                    data=f,
                    file_name="INDRA_projection.html",
                    mime="text/html"
                )
        else:
            st.error("HTML-Plot wurde nicht gefunden.")

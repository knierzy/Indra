import streamlit as st
from pathlib import Path
import subprocess
import shutil
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

if uploaded_file is not None:
    input_path = DATA_DIR / "compendium.xlsx"

    with open(input_path, "wb") as f:
        f.write(uploaded_file.getbuffer())

    st.success("Datei hochgeladen.")

    if st.button("Diagramm erzeugen"):
        with st.spinner("A1 läuft: Daten berechnen und Cartesian Product erzeugen..."):
            result_a1 = subprocess.run(
                ["python", "A1.py"],
                cwd=BASE_DIR,
                capture_output=True,
                text=True
            )

        if result_a1.returncode != 0:
            st.error("Fehler in A1.py")
            st.code(result_a1.stderr)
            st.stop()

        st.success("A1 fertig.")

        # A1 schreibt nach outputs, A2 erwartet aktuell data.
        # Deshalb kopieren wir die erzeugten Dateien für A2.
        shutil.copy(
            OUTPUT_DIR / "CartesianProduct_constraints.xlsx",
            DATA_DIR / "CartesianProduct_constraints.xlsx"
        )

        shutil.copy(
            OUTPUT_DIR / "compendium_processed.xlsx",
            DATA_DIR / "compendium_processed.xlsx"
        )

        with st.spinner("A2 läuft: Diagramm erzeugen..."):
            result_a2 = subprocess.run(
                ["python", "A2.py"],
                cwd=BASE_DIR,
                capture_output=True,
                text=True
            )

        if result_a2.returncode != 0:
            st.error("Fehler in A2.py")
            st.code(result_a2.stderr)
            st.stop()

        st.success("Diagramm erzeugt.")

        html_path = OUTPUT_DIR / "Metanumber_Plot_Ca_HCO3_Bands.html"

        if html_path.exists():
            html = html_path.read_text(encoding="utf-8")
            components.html(html, height=950, scrolling=True)
        else:
            st.error("HTML-Plot wurde nicht gefunden.")

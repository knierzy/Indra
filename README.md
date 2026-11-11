
# INDRA Projection

INDRA Projection is a deterministic two-dimensional visualization framework for multivariate compositional data, developed to improve subgroup discrimination and demonstrated in its present form using hydrochemical datasets.

The method uses integer-encoded compositions and adjustable irrational numerical bases (e.g., e + 14) to generate plotting coordinates.

## Installation

Requires Python 3.12.

```bash
git clone https://github.com/knierzy/Indra.git
cd Indra
pip install -r requirements.txt
python -m playwright install chromium
```

## Files

- **A1.py** – Preprocessing of hydrochemical input data, including ion normalization and metanumber generation.
- **A2.py** – Generation and visualization of the INDRA Projection.
- **indra_pro.py** – Interactive Streamlit implementation.
- **B1.xlsx** – Hydrochemical input dataset used for the INDRA analyses.
- **B2.xlsx** – Log-Euclidean distance matrix.
- **B3.xlsx** – Aitchison distance matrix.
- **B4.xlsx** – Mahalanobis distance matrix.

## Usage

Run the preprocessing script followed by the visualization script:

```bash
python A1.py
python A2.py
```

A1.py reads `data/B1.xlsx` and generates `compendium_processed.xlsx` and `CartesianProduct_constraints.xlsx` in the `outputs` directory.

A2.py reads these generated files and creates the INDRA Projection as interactive HTML and PNG files in the same directory.

## Web Application

An interactive version is available at:

https://indrapro.streamlit.app

The application supports Excel uploads, automatic data processing, adjustable transformation bases, and interactive visualization.
The included B1.xlsx dataset can be uploaded to the web application to test the INDRA Projection.

Note: The application may take a moment to start after a period of inactivity.

## Author

Wolfgang Knierzinger  
Federal Agency for Water Management, Austria

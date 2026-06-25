"""
This script processes hydrochemical Excel data and generates:

- cleaned ion datasets
- HCO3 calculations from ANC
- meq/L and percentage compositions
- 5–95 % filtered typical hydrochemical ranges
  (except for DA_* reference groups)
- constrained Cartesian ion combinations
- hydrochemical meta numbers

The generated meta-number datasets are later used
for INDRA projection and hydrochemical pattern analysis.
"""

import re
import numpy as np
import pandas as pd
import itertools

import streamlit as st
from pathlib import Path
import tempfile

st.set_page_config(layout="wide")
st.markdown(
    """
    <style>
    .block-container {
        max-width: none;
        padding-left: 1rem;
        padding-right: 1rem;
    }
    </style>
    """,
    unsafe_allow_html=True
)
st.title("INDRA Projection")

BASE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = BASE_DIR / "outputs"
OUTPUT_DIR.mkdir(exist_ok=True)

uploaded_file = st.file_uploader(
    "Compendium Excel-Datei hochladen",
    type=["xlsx"]
)

if uploaded_file is None:
    st.info("Bitte eine Excel-Datei hochladen.")
    st.stop()

with tempfile.NamedTemporaryFile(delete=False, suffix=".xlsx") as tmp:
    tmp.write(uploaded_file.getbuffer())
    input_file = Path(tmp.name)

preferred_sheet = "Sheet1"

output_file = OUTPUT_DIR / "compendium_processed.xlsx"
output_file_cartesian = OUTPUT_DIR / "CartesianProduct_constraints.xlsx"

if not st.button("Diagramm erzeugen"):
    st.stop()
# Helper Functions

def find_header_row(xls_path, sheet, probes=5):
    """
    Detects the most likely header row in the Excel sheet
    by searching the first rows for hydrochemical keywords.
    """

    probe = pd.read_excel(
        xls_path,
        sheet_name=sheet,
        header=None,
        nrows=probes
    )

    patterns = [
        r'calcium',
        r'magnesium',
        r'sodium',
        r'potassium',
        r'nitrate',
        r'chloride',
        r'sulfate',
        r'bicarbonate|hydrogencarbonate|hco3',
        r'acid\s*neutralizing|anc',
        r'\bph\b'
    ]

    for i in range(len(probe)):
        row_str = " | ".join(
            probe.iloc[i].fillna("").map(str).tolist()
        ).lower()

        if any(re.search(p, row_str) for p in patterns):
            return i

    return 0

# Find matching column names

def pick(colnames, *regexes):
    """
    Searches column names using multiple regex patterns
    and returns the first matching column.
    """

    for r in regexes:

        rx = re.compile(r, flags=re.IGNORECASE)

        for c in colnames:

            if isinstance(c, str) and rx.search(c):
                return c

    return None


# Convert laboratory text values to numeric values

def to_num(s):
    """
    Converts textual laboratory values into numeric values.

    Handles:
    - detection limit values like '<1'
    - values above limits like '>5'
    - bracketed values '[0.3]'
    - ranges like '0.2-0.4'

    Rules:
    - values below detection limit (<x) → x / 2
    - values above threshold (>x) → x
    - ranges → arithmetic mean
    """

    if pd.isna(s):
        return np.nan

    if isinstance(s, str):

        # Standardize decimal separator
        s = s.strip().replace(",", ".")

        # Below detection limit: <x → x/2
        if s.startswith("<"):

            try:
                return float(s[1:]) / 2

            except:
                return np.nan

        # Above threshold: >x → x
        if s.startswith(">"):

            try:
                return float(s[1:])

            except:
                return np.nan

        # Remove brackets
        s = re.sub(r'[\[\]\(\)]', '', s)

        # Range value: 0.1-0.3 → mean value
        if "-" in s:

            parts = s.split("-")

            try:
                vals = [
                    float(p)
                    for p in parts
                    if p.strip() != ""
                ]

                if len(vals) == 2:
                    return sum(vals) / 2

                elif len(vals) == 1:
                    return vals[0]

            except:
                return np.nan

    # Standard numeric conversion
    return pd.to_numeric(s, errors="coerce")


def hco3_quick_from_sbv(sbv_mmolL):
    return sbv_mmolL * 61.016

def hco3_precise_from_sbv_ph(sbv_mmolL, pH, pK1=6.35, pK2=10.33):
    H  = np.power(10.0, -pH)
    OH = np.power(10.0,  pH - 14.0)
    a1 = np.power(10.0, (pH - pK1))
    a2 = np.power(10.0, (2.0*pH - pK1 - pK2))
    den = 1.0 + a1 + a2
    alpha1 = a1 / den
    alpha2 = a2 / den
    CT_mmolL = (sbv_mmolL - OH*1000.0 + H*1000.0) / (alpha1 + 2.0*alpha2)
    hco3_mmolL = alpha1 * CT_mmolL
    co3_mmolL  = alpha2 * CT_mmolL
    hco3_mgL = hco3_mmolL * 61.016
    co3_mgL  = co3_mmolL  * 60.008
    return hco3_mgL, co3_mgL

def mgL_to_meqL(series, mol_mass, charge):
    return (series / mol_mass) * charge

def norm_col(c):
    if not isinstance(c, str):
        return c
    return c.replace('\xa0', ' ').strip()


def normalize_row(row):
    vals = [row.get(f'Anteil_%_{ion}', np.nan) for ion in ionen]
    if any(pd.isna(vals)) or row.get('Summe_Gesamt_meq_L', np.nan) <= 0:
        return pd.Series([np.nan]*len(ionen), index=[f'Anteil_int_%_{ion}' for ion in ionen])
    arr = np.array(vals, dtype=float)
    flo = np.floor(arr).astype(int)
    rest = int(100 - flo.sum())
    if rest != 0:
        diffs = arr - flo
        order = np.argsort(-diffs)  # größte Nachkommastellen zuerst
        for i in range(min(abs(rest), len(order))):
            flo[order[i]] += 1 if rest > 0 else -1
    return pd.Series(flo, index=[f'Anteil_int_%_{ion}' for ion in ionen])



# Read Excel file

xls = pd.ExcelFile(input_file)
sheet = preferred_sheet if preferred_sheet in xls.sheet_names else xls.sheet_names[0]
hrow = find_header_row(input_file, sheet)

df_all = pd.read_excel(input_file, sheet_name=sheet, header=None)
header = df_all.iloc[hrow].astype(str).map(norm_col)

df = df_all.iloc[hrow + 1:].copy()
df.columns = header
df.columns = [norm_col(c) for c in df.columns]

df = df.dropna(how="all").copy()
df = df.reset_index(drop=True)


# Create ID column from the first column

df['ID'] = df.iloc[:, 0].astype(str).str.strip()
print("✅ ID extracted from the first column.")


# Create measurement station column from the second column

station_col_raw = df.columns[1]
df['Messstation'] = df.iloc[:, 1].astype(str).str.strip()

print(f"✅ Measurement station extracted from column '{station_col_raw}'.")


# Create municipality name column from the fourth column

gemeinde_col_raw = df.columns[3]
df['Gemeindename'] = df.iloc[:, 3].astype(str).str.strip()

print(f"✅ Municipality name extracted from column '{gemeinde_col_raw}'.")


print(f"📑 Sheets: {xls.sheet_names}")
print(f"➡️ Used sheet: {sheet}")
print(f"🧭 Used header row: {hrow}")


# Find relevant columns

cols = list(df.columns)

mapping = {
    'SAMPLING_DATE': pick(cols, r'G102|sampling|entnahme|date'),
    'CALCIUM_mg_L': pick(cols, r'G134|calcium'),
    'MAGNESIUM_mg_L': pick(cols, r'G135|magnesium'),
    'SODIUM_mg_L': pick(cols, r'G136|sodium|natrium'),
    'POTASSIUM_mg_L': pick(cols, r'G137|potassium|kalium'),
    'NITRATE_mg_L': pick(cols, r'G154|nitrate|nitrat'),
    'CHLORIDE_mg_L': pick(cols, r'G155|chloride|chlorid'),
    'SULFATE_mg_L': pick(cols, r'G156|sulfate|sulfat'),
    'BICARBONATE_mg_L': pick(cols, r'G157|bicarbonate|hydrogencarbonate|hydrogenkarbonat|hco3'),
    'ACID_NEUTRALIZING_CAPACITY': pick(cols, r'G158|acid\s*neutralizing|anc|sbv'),
    'pH': pick(cols, r'\bph\b'),
}

print("\n🔎 Column mapping:")
for k, v in mapping.items():
    print(f"{k:30s} -> {v}")

required_mapping = [
    "CALCIUM_mg_L",
    "MAGNESIUM_mg_L",
    "SODIUM_mg_L",
    "POTASSIUM_mg_L",
    "NITRATE_mg_L",
    "CHLORIDE_mg_L",
    "SULFATE_mg_L",
    "BICARBONATE_mg_L",
]

missing_mapping = [k for k in required_mapping if mapping[k] is None]
if missing_mapping:
    raise ValueError(f"❌ Missing column mapping: {missing_mapping}")
# Create Art column from column B

df['Art'] = df.iloc[:, 1].astype(str).str.strip()
print("✅ Group column created.")


print("\n🔎 Group size by Art:")
print(df.groupby('Art').size().describe())



print("\n🔎 Group sizes for percentile filtering:")


print("\n🔎 Column mapping:")
for k, v in mapping.items():
    print(f"  {k:15s} → {v if (isinstance(v, str) and v in df.columns) else str(v)}")



# Calculate bicarbonate from ANC or use existing bicarbonate column

anc_col = mapping['ACID_NEUTRALIZING_CAPACITY']
hco3_col = mapping['BICARBONATE_mg_L']

if anc_col is not None:
    df['ANC_mmol_L'] = df[anc_col].apply(to_num)
else:
    df['ANC_mmol_L'] = np.nan


# Use existing bicarbonate values if available

if hco3_col:
    df['HCO3_mg_L_original'] = df[hco3_col].apply(to_num)
else:
    df['HCO3_mg_L_original'] = np.nan


# Calculate HCO3 only if ANC is available

if anc_col is not None:

    df['HCO3_mg_L_quick'] = hco3_quick_from_sbv(
        df['ANC_mmol_L']
    )

    if mapping['pH']:

        df['pH_num'] = df[mapping['pH']].apply(to_num)

        hco3_precise = []
        co3_precise = []

        for anc, pH in zip(
            df['ANC_mmol_L'],
            df['pH_num']
        ):

            if pd.notna(anc) and pd.notna(pH):

                hco3_val, co3_val = (
                    hco3_precise_from_sbv_ph(
                        anc,
                        pH
                    )
                )

                hco3_precise.append(hco3_val)
                co3_precise.append(co3_val)

            else:

                hco3_precise.append(np.nan)
                co3_precise.append(np.nan)

        df['HCO3_mg_L_precise'] = hco3_precise
        df['CO3_mg_L_precise'] = co3_precise

    else:

        df['HCO3_mg_L_precise'] = np.nan

else:

    df['HCO3_mg_L_quick'] = np.nan
    df['HCO3_mg_L_precise'] = np.nan


# Select final bicarbonate value

df['HCO3_mg_L_final'] = df['HCO3_mg_L_original']

df['HCO3_mg_L_final'] = (
    df['HCO3_mg_L_final']
    .fillna(df['HCO3_mg_L_precise'])
)

df['HCO3_mg_L_final'] = (
    df['HCO3_mg_L_final']
    .fillna(df['HCO3_mg_L_quick'])
)


# Extract ion concentrations in mg/L

df['Ca_mg_L'] = df[mapping['CALCIUM_mg_L']].apply(to_num) if mapping['CALCIUM_mg_L'] else np.nan
df['Mg_mg_L'] = df[mapping['MAGNESIUM_mg_L']].apply(to_num) if mapping['MAGNESIUM_mg_L'] else np.nan
df['Na_mg_L'] = df[mapping['SODIUM_mg_L']].apply(to_num) if mapping['SODIUM_mg_L'] else np.nan
df['K_mg_L'] = df[mapping['POTASSIUM_mg_L']].apply(to_num) if mapping['POTASSIUM_mg_L'] else np.nan
df['Cl_mg_L'] = df[mapping['CHLORIDE_mg_L']].apply(to_num) if mapping['CHLORIDE_mg_L'] else np.nan
df['SO4_mg_L'] = df[mapping['SULFATE_mg_L']].apply(to_num) if mapping['SULFATE_mg_L'] else np.nan

no3_raw = df[mapping['NITRATE_mg_L']].apply(to_num) if mapping['NITRATE_mg_L'] else np.nan


# Convert nitrate-nitrogen to nitrate only for lake samples

df['NO3_mg_L'] = np.where(
    df['ID'].str.startswith("Lake", na=False),
    no3_raw * (62.0 / 14.0),
    no3_raw
)


# Use final HCO3 value as standardized HCO3 concentration

df['HCO3_mg_L'] = df['HCO3_mg_L_final']


# Keep only complete ion analyses

required_ions = [
    'Ca_mg_L',
    'Mg_mg_L',
    'Na_mg_L',
    'K_mg_L',
    'NO3_mg_L',
    'Cl_mg_L',
    'SO4_mg_L',
    'HCO3_mg_L'
]

print("\n🔎 Availability before complete-ion filter:")
for c in required_ions:
    print(c, df[c].notna().sum())

print("\n🔎 Bicarbonate sources:")
print("HCO3 original:", df['HCO3_mg_L_original'].notna().sum())
print("ANC:", df['ANC_mmol_L'].notna().sum())
print("HCO3 quick:", df['HCO3_mg_L_quick'].notna().sum())
print("HCO3 final:", df['HCO3_mg_L_final'].notna().sum())





n_before = len(df)

df = df.dropna(subset=required_ions).copy()
if df.empty:
    raise ValueError(
        "No complete ion analyses found. Check the column mapping above: "
        "at least one required ion column was not detected."
    )


n_after = len(df)

print(
    f"🧪 Complete ion analyses: "
    f"{n_after} of {n_before} rows retained "
    f"({n_after / n_before:.1%})"
)


# Convert ion concentrations from mg/L to meq/L

df['meq_L_Ca2+'] = mgL_to_meqL(df['Ca_mg_L'], 40.078, 2)
df['meq_L_Mg2+'] = mgL_to_meqL(df['Mg_mg_L'], 24.305, 2)
df['meq_L_Na+'] = mgL_to_meqL(df['Na_mg_L'], 22.990, 1)
df['meq_L_K+'] = mgL_to_meqL(df['K_mg_L'], 39.098, 1)

df['meq_L_Cl-'] = mgL_to_meqL(df['Cl_mg_L'], 35.453, 1)
df['meq_L_SO4_2-'] = mgL_to_meqL(df['SO4_mg_L'], 96.06, 2)
df['meq_L_NO3-'] = mgL_to_meqL(df['NO3_mg_L'], 62.0049, 1)
df['meq_L_HCO3-'] = mgL_to_meqL(df['HCO3_mg_L'], 61.016, 1)


# Calculate cation and anion sums

df['Sum_Kationen_meq_L'] = df[
    ['meq_L_Ca2+', 'meq_L_Mg2+', 'meq_L_Na+', 'meq_L_K+']
].sum(axis=1)

df['Sum_Anionen_meq_L'] = df[
    ['meq_L_Cl-', 'meq_L_SO4_2-', 'meq_L_NO3-', 'meq_L_HCO3-']
].sum(axis=1)


# Calculate charge balance error

den = df['Sum_Kationen_meq_L'] + df['Sum_Anionen_meq_L']

df['Bilanzfehler_%'] = np.where(
    den > 0,
    np.abs(df['Sum_Kationen_meq_L'] - df['Sum_Anionen_meq_L']) / den * 100,
    np.nan
)


# Create sheet data with all measurements in meq/L and percent

df_meq_pct = df.copy()

df_meq_pct['Summe_Gesamt_meq_L'] = (
    df_meq_pct['meq_L_Ca2+'] +
    df_meq_pct['meq_L_Mg2+'] +
    df_meq_pct['meq_L_Na+'] +
    df_meq_pct['meq_L_K+'] +
    df_meq_pct['meq_L_Cl-'] +
    df_meq_pct['meq_L_SO4_2-'] +
    df_meq_pct['meq_L_NO3-'] +
    df_meq_pct['meq_L_HCO3-']
)

ionen = [
    'Ca2+',
    'Mg2+',
    'Na+',
    'K+',
    'Cl-',
    'SO4_2-',
    'NO3-',
    'HCO3-'
]


# Calculate percentage share of each ion

for ion in ionen:

    num = df_meq_pct[f'meq_L_{ion}']
    den = df_meq_pct['Summe_Gesamt_meq_L']

    df_meq_pct[f'Anteil_%_{ion}'] = np.where(
        den > 0,
        (num / den) * 100,
        np.nan
    )


# Calculate integer percentage shares using the largest remainder method

df_meq_pct[[f'Anteil_int_%_{ion}' for ion in ionen]] = (
    df_meq_pct.apply(normalize_row, axis=1)
)

# Apply 5–95 percent filter per GZÜV ID

ions_mgL_for_filter = [
    'Ca_mg_L',
    'Mg_mg_L',
    'Na_mg_L',
    'K_mg_L',
    'Cl_mg_L',
    'SO4_mg_L',
    'NO3_mg_L',
    'HCO3_mg_L'
]

df['__keep__'] = True

for gid, g in df.groupby('Art', dropna=False):

    # Case 1: DA_* groups are not filtered

    if isinstance(gid, str) and gid.startswith("DA_"):

        df.loc[g.index, '__keep__'] = True
        continue

    # Case 2: standard groups are filtered by the 5–95 percent range

    mask_g = pd.Series(True, index=g.index)

    for col in ions_mgL_for_filter:

        if col in g.columns:

            n_nonnull = g[col].notna().sum()

            if n_nonnull >= 3:

                lo = g[col].quantile(0.05)
                hi = g[col].quantile(0.95)

                mask_g &= (g[col].between(lo, hi)) | (g[col].isna())

            else:

                mask_g &= True

    df.loc[g.index, '__keep__'] = mask_g


df_typisch = df[df['__keep__']].drop(columns='__keep__').copy()

print(
    f"🧪 Typical data per Art, 5–95 percent range: "
    f"{len(df_typisch)} of {len(df)} rows retained"
)


# Add percentage shares based on meq/L

df_typisch['Summe_Gesamt_meq_L'] = (
    df_typisch['meq_L_Ca2+'] +
    df_typisch['meq_L_Mg2+'] +
    df_typisch['meq_L_Na+'] +
    df_typisch['meq_L_K+'] +
    df_typisch['meq_L_Cl-'] +
    df_typisch['meq_L_SO4_2-'] +
    df_typisch['meq_L_NO3-'] +
    df_typisch['meq_L_HCO3-']
)

ionen = [
    'Ca2+',
    'Mg2+',
    'Na+',
    'K+',
    'Cl-',
    'SO4_2-',
    'NO3-',
    'HCO3-'
]

for ion in ionen:

    num = df_typisch[f'meq_L_{ion}']
    den_all = df_typisch['Summe_Gesamt_meq_L']

    df_typisch[f'Anteil_%_{ion}'] = np.where(
        den_all > 0,
        (num / den_all) * 100,
        np.nan
    )


# Calculate integer percentage shares using the largest remainder method

df_typisch[[f'Anteil_int_%_{ion}' for ion in ionen]] = (
    df_typisch.apply(normalize_row, axis=1)
)


# Calculate correlation matrix values per Art

ions = [
    "Ca2+",
    "Mg2+",
    "Na+",
    "K+",
    "Cl-",
    "SO4_2-",
    "NO3-",
    "HCO3-"
]

corr_results = []

for gid, g in df_typisch.groupby("Art"):

    cols = [f"Anteil_%_{ion}" for ion in ions]
    df_sub = g[cols].dropna()

    if len(df_sub) < 3:
        continue

    corr = df_sub.corr()

    for i in range(len(cols)):

        for j in range(i + 1, len(cols)):

            ionA = cols[i].replace("Anteil_%_", "")
            ionB = cols[j].replace("Anteil_%_", "")

            corr_val = corr.iloc[i, j]

            corr_results.append({
                "Art": gid,
                "Ion_A": ionA,
                "Ion_B": ionB,
                "Correlation": round(corr_val, 3)
            })


df_corr = pd.DataFrame(corr_results)

# Create Sheet 3 with min/max values of integer ion percentages per locality

df_typisch = df_typisch.copy()



# No depth filter is applied, all typical data are used

df_typ_bis10 = df_typisch.copy()


# Check whether the ID column exists

if 'ID' not in df_typ_bis10.columns:
    raise ValueError("❌ Column 'ID' not found. Please check whether it exists in the input file.")


ionen = [
    'Ca2+',
    'Mg2+',
    'Na+',
    'K+',
    'Cl-',
    'SO4_2-',
    'NO3-',
    'HCO3-'
]

int_cols = [f'Anteil_int_%_{ion}' for ion in ionen]

for c in int_cols:

    if c not in df_typ_bis10.columns:
        df_typ_bis10[c] = np.nan

    df_typ_bis10[c] = pd.to_numeric(df_typ_bis10[c], errors='coerce')


# Calculate min/max values of integer percentage shares per measurement station

agg_dict_int = {c: ['min', 'max'] for c in int_cols}

minmax_typisch = (
    df_typ_bis10
    .groupby('Art', dropna=False)
    .agg(agg_dict_int)
)


# Flatten multi-level column names

minmax_typisch.columns = [
    f"{c}_{stat}" for c, stat in minmax_typisch.columns
]

minmax_typisch = minmax_typisch.reset_index()


# Add municipality name back to each measurement station

gemeinde_map = (
    df_typ_bis10[['Art', 'Gemeindename']]
    .dropna()
    .drop_duplicates(subset='Art')
)

minmax_typisch = minmax_typisch.merge(
    gemeinde_map,
    on='Art',
    how='left'
)


# Move municipality name to the first column

cols = ['Gemeindename'] + [
    c for c in minmax_typisch.columns
    if c != 'Gemeindename'
]

minmax_typisch = minmax_typisch[cols]


# Count number of typical records per Art

counts = (
    df_typ_bis10
    .groupby('Art', dropna=False)
    .size()
    .rename('n_typisch')
)


# Calculate mean charge balance error per Art

bilanz_mean = (
    df_typ_bis10
    .groupby('Art', dropna=False)['Bilanzfehler_%']
    .mean()
    .rename('Bilanzfehler_mean_%')
)


# Add record count and mean charge balance error to min/max table

minmax_typisch = (
    minmax_typisch
    .merge(counts, on='Art', how='left')
    .merge(bilanz_mean, on='Art', how='left')
)


# Export processed data to Excel

with pd.ExcelWriter(output_file, engine="openpyxl") as writer:

    df.to_excel(
        writer,
        sheet_name="Calculated_Data",
        index=False
    )

    df_typisch.to_excel(
        writer,
        sheet_name="Typical_Data_5_95",
        index=False
    )

    minmax_typisch.to_excel(
        writer,
        sheet_name="MinMax_Typical",
        index=False
    )

    df_meq_pct.to_excel(
        writer,
        sheet_name="All_meq_Percent",
        index=False
    )

    df_corr.to_excel(
        writer,
        sheet_name="Correlation_Pairs",
        index=False
    )


print(
    f"\n✅ File saved with calculated data, grouped data, "
    f"min/max values, meq/percentage data, and correlation pairs: {output_file}"
)


print("\n📦 Creating constrained Cartesian product and meta numbers ...")


# Define ion sets

ion_pairs_kat = ["Ca2+", "Mg2+", "Na+", "K+"]
ion_pairs_ani = ["HCO3-", "SO4_2-", "Cl-", "NO3-"]

all_ions = ion_pairs_kat + ion_pairs_ani


# Learn constraint bands between ion pairs

def lern_constraint_bands(df_typ):
    """
    Learns constraints between ion pairs.

    The constraints are based on:
    - floating-point percentage values
    - robust 95% quantile limits instead of absolute maxima
    - positive correlations between ion pairs
    """

    ions = [
        "Ca2+",
        "Mg2+",
        "Na+",
        "K+",
        "HCO3-",
        "SO4_2-",
        "Cl-",
        "NO3-"
    ]

    constraints = {}

    for gid, g in df_typ.groupby("Art"):

        c = {"Art": gid}

        for i in range(len(ions)):

            for j in range(i + 1, len(ions)):

                ionA = ions[i]
                ionB = ions[j]

                # Use floating-point percentages instead of integer percentages

                colA = f"Anteil_%_{ionA}"
                colB = f"Anteil_%_{ionB}"

                key = f"{ionA}__{ionB}"

                if colA not in g or colB not in g:

                    c[f"{key}_lo"] = None
                    c[f"{key}_hi"] = None
                    continue

                df_pair = g[[colA, colB]].dropna()

                if df_pair.shape[0] < 5:

                    c[f"{key}_lo"] = None
                    c[f"{key}_hi"] = None
                    continue

                # Calculate absolute difference between paired ion percentages

                diff = (df_pair[colA] - df_pair[colB]).abs()

                # Calculate correlation between both ions

                corr = df_pair.corr().iloc[0, 1]

                # Report strong negative correlations, but do not use them as constraints

                if corr is not None and corr < -0.6:
                    print(
                        f"⚠️ Strong negative correlation ignored: "
                        f"{ionA} vs {ionB} (r={corr:.2f})"
                    )

                # Use robust upper boundary based on the 95% quantile

                hi_robust = diff.quantile(0.95)

                # Only positive correlations are used as coupling constraints

                if corr is not None and corr >= 0.3:

                    if corr > 0.6:
                        lo = 0
                        hi = hi_robust

                    else:
                        lo = 0
                        hi = hi_robust * 1.1

                else:

                    lo = None
                    hi = None

                # Round up upper boundary to protect against rounding effects

                if hi is not None:
                    hi = float(np.ceil(hi))

                c[f"{key}_lo"] = lo
                c[f"{key}_hi"] = hi

        constraints[gid] = c

    return constraints


constraints_by_gid = lern_constraint_bands(df_typisch)


# Check whether a generated ion combination satisfies learned constraints

def ok_constraints(gid, combo_dict):
    """
    Checks whether an ion combination satisfies the learned constraints.
    """

    c = constraints_by_gid.get(gid, None)

    if c is None:
        return True

    for key, val in c.items():

        if key.endswith("_lo"):

            pair = key[:-3]
            ionA, ionB = pair.split("__")

            lo = val
            hi = c.get(f"{pair}_hi", None)

            if lo is None or hi is None:
                continue

            diff = abs(combo_dict[ionA] - combo_dict[ionB])

            if not (lo <= diff <= hi):
                return False

    return True


# Generate Cartesian product and filter by constraints

all_results = []

for _, row in minmax_typisch.iterrows():

    gid = row["Art"]
    gemeinde = row.get("Gemeindename", "")


    # Build ion ranges from min/max values

    ranges_loc = {}

    for ion in all_ions:

        min_val = row.get(f"Anteil_int_%_{ion}_min")
        max_val = row.get(f"Anteil_int_%_{ion}_max")

        if pd.isna(min_val) or pd.isna(max_val):

            print(f"❌ Skipping {gid} because of NaN in {ion}")
            ranges_loc = None
            break

        ranges_loc[ion] = range(int(min_val), int(max_val) + 1)


    if ranges_loc is None:
        continue


    print(
        f"👉 {gid}: range sizes:",
        {k: len(v) for k, v in ranges_loc.items()}
    )


    combos = itertools.product(*ranges_loc.values())


    # Calculate size of the raw search space

    total_raw = 1

    for r in ranges_loc.values():
        total_raw *= len(r)


    count_sum_ok = 0
    count_final = 0

    valid = []

    for combo in combos:

        # Keep only combinations where all ion percentages sum to 100

        if sum(combo) != 100:
            continue

        count_sum_ok += 1

        d = dict(zip(all_ions, combo))


        # Apply learned ion-pair constraints

        if not ok_constraints(gid, d):
            continue

        count_final += 1

        valid.append(combo)


    if len(valid) == 0:
        continue


    print(f"\n📍 {gid} ({gemeinde})")
    print(f"  🔢 Raw search space:     {total_raw:,}")
    print(f"  ➗ Sum = 100:            {count_sum_ok:,}")
    print(f"  🔗 After constraints:    {count_final:,}")

    if total_raw > 0:

        print(f"  📉 Reduction ratio:      {count_final / total_raw:.8f}")


    df_loc = pd.DataFrame(valid, columns=all_ions)

    df_loc["Art"] = gid
    df_loc["Gemeindename"] = gemeinde


    # Generate meta numbers

    def erzeuge_metazahl(row, spalten):
        return int("".join(f"{int(row[col]):02}" for col in spalten))


    df_loc["Metazahl_Kationen"] = df_loc.apply(
        lambda r: erzeuge_metazahl(r, ion_pairs_kat),
        axis=1
    )

    df_loc["Metazahl_Anionen"] = df_loc.apply(
        lambda r: erzeuge_metazahl(r, ion_pairs_ani),
        axis=1
    )


    all_results.append(df_loc)

# Collect all generated results

if all_results:
    df_cartesian = pd.concat(all_results, ignore_index=True)
else:
    df_cartesian = pd.DataFrame()


# Segment analysis

def split_metazahl(x):
    s = str(int(x)).zfill(8)
    return int(s[0:2]), int(s[2:4]), int(s[4:6]), int(s[6:8])


if not df_cartesian.empty:

    df_cartesian[["Kat1", "Kat2", "Kat3", "Kat4"]] = (
        df_cartesian["Metazahl_Kationen"]
        .apply(lambda x: pd.Series(split_metazahl(x)))
    )

    df_cartesian[["Ani1", "Ani2", "Ani3", "Ani4"]] = (
        df_cartesian["Metazahl_Anionen"]
        .apply(lambda x: pd.Series(split_metazahl(x)))
    )

    kation_labels = ["Ca", "Mg", "Na", "K"]
    anion_labels = ["HCO3", "SO4", "Cl", "NO3"]

    results = []

    for gid, g in df_cartesian.groupby("Art"):

        gemeinde = g["Gemeindename"].iloc[0]

        kat_means = {
            lab: g[f"Kat{i + 1}"].mean()
            for i, lab in enumerate(kation_labels)
        }

        ani_means = {
            lab: g[f"Ani{i + 1}"].mean()
            for i, lab in enumerate(anion_labels)
        }

        results.append({
            "Art": gid,
            "Gemeindename": gemeinde,
            "Max_Kation_Segment": max(kat_means, key=kat_means.get),
            "Max_Kation_Mean": round(max(kat_means.values()), 2),
            "Max_Anionen_Segment": max(ani_means, key=ani_means.get),
            "Max_Anionen_Mean": round(max(ani_means.values()), 2),
        })

    df_segstats = pd.DataFrame(results)

else:

    df_segstats = pd.DataFrame()


# Export Cartesian product results

with pd.ExcelWriter(output_file_cartesian, engine="openpyxl") as writer:

    if not df_cartesian.empty:

        df_cartesian.to_excel(
            writer,
            sheet_name="Meta_Kombinationen",
            index=False
        )

    df_segstats.to_excel(
        writer,
        sheet_name="Segment_Maxima",
        index=False
    )


print("✔️ Constrained Cartesian product file saved.")
print("📁", output_file_cartesian)


# Das ist log euclidean... passe noch an
#
# -*- coding: utf-8 -*-
import math
import pandas as pd
import plotly.graph_objects as go
from plotly.express.colors import qualitative
import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial import distance_matrix
from scipy.spatial.distance import euclidean
import streamlit.components.v1 as components
def log_euclid(a, b):
    return euclidean(np.log1p(a), np.log1p(b))



from pathlib import Path

input_file = output_file_cartesian
raw_file = output_file
plot_output = OUTPUT_DIR / "Metanumber_Plot_Ca_HCO3_Bands.html"
df = pd.read_excel(
    input_file,
    sheet_name="Meta_Kombinationen"
)

raw_df = pd.read_excel(
    raw_file,
    sheet_name="Typical_Data_5_95"
)



# ============================================================
# --- Hilfsfunktion zur Transformation mit frei wählbarer Basis ---
def custom_transform_optimal(x, base=math.e +4): #12.1415926535
    try:
        x_str = str(int(x)).zfill(8)
        a = int(x_str[0:2])
        b = int(x_str[2:4])
        c = int(x_str[4:6])
        d = int(x_str[6:8])
        return a * base**3 + b * base**2 + c * base + d
    except:
        return None

# --- Hover-Helfer ---
def pairs_to_percentages(x, labels):
    try:
        s = str(int(x)).zfill(8)
        vals = [int(s[i:i+2]) for i in range(0, 8, 2)]
        return dict(zip(labels, vals)), vals
    except:
        return dict(zip(labels, [None]*4)), [None]*4

def format_hover(row):
    k_labels = ["Ca", "Mg", "Na", "K"]
    a_labels = ["HCO₃", "SO₄", "Cl", "NO₃"]
    k_perc, _ = pairs_to_percentages(row["Metazahl_Kationen"], k_labels)
    a_perc, _ = pairs_to_percentages(row["Metazahl_Anionen"], a_labels)
    k_lines = " · ".join([f"{lbl}: {k_perc[lbl]}%" if k_perc[lbl] is not None else f"{lbl}: –" for lbl in k_labels])
    a_lines = " · ".join([f"{lbl}: {a_perc[lbl]}%" if a_perc[lbl] is not None else f"{lbl}: –" for lbl in a_labels])
    return (
        f"<b>Art:</b> {row['Art']}<br>"
        f"<b>Kationen</b> (aus {str(row['Metazahl_Kationen']).zfill(8)}):<br>{k_lines}<br>"
        f"<b>Anionen</b> (aus {str(row['Metazahl_Anionen']).zfill(8)}):<br>{a_lines}"
    )

try:
    # Excel einlesen
    df = pd.read_excel(input_file, sheet_name="Meta_Kombinationen")

    required_cols = ['Metazahl_Kationen', 'Metazahl_Anionen', 'Art']
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Spalte '{col}' fehlt in der Datei!")

    # Transformation anwenden
    df["Kationen_trans_raw"] = df["Metazahl_Kationen"].apply(custom_transform_optimal)
    df["Anionen_trans_raw"]  = df["Metazahl_Anionen"].apply(custom_transform_optimal)

    # Normierung 0–100
    df["Kationen_trans"] = df["Kationen_trans_raw"] / df["Kationen_trans_raw"].max() * 100
    df["Anionen_trans"]  = df["Anionen_trans_raw"]  / df["Anionen_trans_raw"].max()  * 100

    # Duplikate entfernen
    df = df.drop_duplicates(subset=["Kationen_trans", "Anionen_trans", "Art"])

    # Überlappungen zählen
    koord_counts = (
        df.groupby(["Kationen_trans", "Anionen_trans"])
        .size()
        .reset_index(name="region_count")
    )
    df = df.merge(koord_counts, on=["Kationen_trans", "Anionen_trans"], how="left")
    df["Symbol"] = df["region_count"].apply(lambda x: "star" if x > 1 else "circle")



    # Hover vorbereiten
    df["hover_text"] = df.apply(format_hover, axis=1)

    # ============================================================
    # MAHALANOBIS-DISTANZ (VOR DEM PLOT!)
    # ============================================================


    from scipy.spatial.distance import mahalanobis

    # --- Ionen definieren ---
    ion_cols = [
        "meq_L_Ca2+",
        "meq_L_Mg2+",
        "meq_L_Na+",
        "meq_L_K+",
        "meq_L_Cl-",
        "meq_L_SO4_2-",
        "meq_L_NO3-",
        "meq_L_HCO3-"
    ]

    # --- Kovarianzmatrix ---
    cov = np.cov(raw_df[ion_cols].values.T)
    cov += np.eye(cov.shape[0]) * 1e-6
    cov_inv = np.linalg.pinv(cov)

    # --- Gruppenmittelwerte ---
    group_means = raw_df.groupby("Art")[ion_cols].mean()
    group_means.index = group_means.index.astype(str).str.strip()

    # ============================================================
    # 🔬 LOG-TRANSFORMIERTE MAHALANOBIS (NEU)
    # ============================================================

    # --- Log-Transformation der Rohdaten ---
    X_log = np.log1p(raw_df[ion_cols])

    # --- Kovarianzmatrix im Log-Raum ---
    cov_log = np.cov(X_log.values.T)
    cov_log += np.eye(cov_log.shape[0]) * 1e-6
    cov_log_inv = np.linalg.pinv(cov_log)

    # --- Gruppenmittelwerte im Log-Raum ---
    group_means_log = raw_df.groupby("Art")[ion_cols].mean()
    group_means_log = np.log1p(group_means_log)
    group_means_log.index = group_means_log.index.astype(str).str.strip()

    # --- DEBUG: Vergleich Hallstatt vs Ossiach ---
    from scipy.spatial.distance import euclidean, mahalanobis

    h_name = next(
        g for g in group_means.index
        if str(g).strip().lower() == "lake hallstatt"
    )

    o_name = next(
        g for g in group_means.index
        if str(g).strip().lower() == "lake ossiach"
    )

    h = group_means.loc[h_name].values
    o = group_means.loc[o_name].values

    print("\n🔍 Vergleich Hallstatt vs Ossiach")
    print("Hallstatt:", h_name)
    print("Ossiach:", o_name)

    print("\nMittelwerte Differenz:")
    print(group_means.loc[h_name] - group_means.loc[o_name])

    print("\nDistanzen:")
    print("Euclidean:   ", euclidean(h, o))
    print("Mahalanobis (raw): ", mahalanobis(o, h, cov_inv))

    h_log = group_means_log.loc[h_name].values
    o_log = group_means_log.loc[o_name].values

    print("Mahalanobis (log): ", mahalanobis(o_log, h_log, cov_log_inv))

    # --- Referenz (Hallstatt) ---
    # --- Referenz (Hallstatt) ---
    ref_group = next(
        g for g in group_means.index
        if str(g).strip().lower() == "lake hallstatt"
    )
    print(f"\n✅ Referenz: {ref_group}")

        # LOG-Version für Plot verwenden
    ref_vector = group_means.loc[ref_group].values

    # Calculate Log-Euclidean distances to reference group
    mah_dict = {}

    for g in group_means.index:
        vec = group_means.loc[g].values
        mah_dict[str(g).strip().lower()] = log_euclid(vec, ref_vector)

        # Clean plot group names
    df["Group_clean"] = df["Art"].astype(str).str.strip().str.lower()

    # LogEuclid direkt über identische Gruppennamen zuordnen
    df["LogEuclid"] = df["Group_clean"].map(mah_dict)

    print("\nLogEuclid Check:")
    print(df["LogEuclid"].head())
    print("NaN Anzahl:", df["LogEuclid"].isna().sum())

    missing = df[df["LogEuclid"].isna()]["Group_clean"].unique()

    print("\n❌ NICHT GEMATCHT:")
    for m in missing[:20]:
        print(m)


    fig = go.Figure()

    mah_sorted = sorted(mah_dict.items(), key=lambda x: x[1])

    mah_text = "<span style='font-size:15px'><b>Log-Euclidan distance – reference: Hallstatt</b></span><br>"

    for g, d in mah_sorted[:5]:  # Top 5
        mah_text += f"{g.title()}: {d:.2f}<br>"

    fig.add_annotation(
        xref="paper", yref="paper",
        x=0.58, y=0.98,
        xanchor="left", yanchor="top",
        text=mah_text,
        showarrow=False,
        font=dict(size=16),  # Basisgröße für den Rest
        align="left",
        bgcolor="rgba(255,255,255,0.95)",
        bordercolor="black", borderwidth=1.5
    )

    print("\n📏 Mahalanobis-Distanzen relativ zu Hallstatt:\n")

    for g, d in sorted(mah_dict.items(), key=lambda x: x[1]):
        print(f"{g:25s}  →  {d:.3f}")



    print("\n📏 Log-Euclidean Distanzen relativ zu Hallstatt:\n")
    for g, d in sorted(mah_dict.items(), key=lambda x: x[1]):
        print(f"{g:25s}  →  {d:.3f}")



    # 🔥 HIER HINZUFÜGEN
    max_maha = df["LogEuclid"].max()
    # 🔥 DEBUG HIER EINBAUEN
    missing = df[df["LogEuclid"].isna()]["Group_clean"].unique()

    print("\n❌ NICHT GEMATCHT:")
    for m in missing[:20]:
        print(m)

    print("\nLogEuclid Check:")
    print(df["LogEuclid"].head())
    print("NaN Anzahl:", df["LogEuclid"].isna().sum())

    print("\nDEBUG MATCHING:")
    print(df["Group_clean"].unique()[:10])
    print(list(mah_dict.keys())[:10])

    print("\nLogEuclid Check:")
    print(df["LogEuclid"].head())
    print("NaN Anzahl:", df["LogEuclid"].isna().sum())

    # === Ca- und HCO3-Werte berechnen ===
    df["Ca_val"] = df["Metazahl_Kationen"].apply(
        lambda x: pairs_to_percentages(x, ["Ca", "Mg", "Na", "K"])[0]["Ca"]
    )

    df["HCO3_val"] = df["Metazahl_Anionen"].apply(
        lambda x: pairs_to_percentages(x, ["HCO₃", "SO₄", "Cl", "NO₃"])[0]["HCO₃"]
    )

    ca_max = df["Ca_val"].max()
    hco3_max = df["HCO3_val"].max()

    results_ca = []
    for ca_val in [2, 5, 10, 15, 20, 25, 30, 35, 40]:
        sub = df[df["Ca_val"] == ca_val]
        if sub.empty:
            continue
        y_min = sub["Kationen_trans_raw"].min() / df["Kationen_trans_raw"].max() * 100
        y_max = sub["Kationen_trans_raw"].max() / df["Kationen_trans_raw"].max() * 100
        results_ca.append(dict(Ca=ca_val, y_min=y_min, y_max=y_max))

        # Beispiel: Ca-Band
        fig.add_trace(go.Scatter(
            x=[0, 100, 100, 0],
            y=[y_min, y_min, y_max, y_max],
            fill="toself",
            fillpattern=dict(
                shape="/",  # Schraffur
                fgcolor="grey",
                size=6,
                solidity=0.08
            ),
            fillcolor="lightgrey",
            line=dict(width=0),
            opacity=0.3,
            name=f"Ca = {ca_val}%",
            showlegend=False,
            hoverinfo="skip"  # kein Hover
        ))

        fig.add_annotation(
            x=0,
            y=(y_min + y_max) / 2,
            text=f"<b>Ca = {ca_val}%</b>",
            showarrow=False,
            font=dict(size=12, color="grey"),
            xanchor="left",
            yanchor="middle"
        )

    # === HCO3-Referenzbänder für 20% und 40% ===
    results_hco3 = []
    for hco3_val in [5, 10, 15, 20, 25, 30, 35, 40, 45]:
        sub = df[df["HCO3_val"] == hco3_val]
        if sub.empty:
            continue
        x_min = sub["Anionen_trans_raw"].min() / df["Anionen_trans_raw"].max() * 100
        x_max = sub["Anionen_trans_raw"].max() / df["Anionen_trans_raw"].max() * 100
        results_hco3.append(dict(HCO3=hco3_val, x_min=x_min, x_max=x_max))

        fig.add_trace(go.Scatter(
            x=[x_min, x_max, x_max, x_min],
            y=[0, 0, 100, 100],
            fill="toself",
            fillpattern=dict(
                shape="\\",  # Schraffur andere Richtung
                fgcolor="blue",
                size=6,
                solidity=0.08
            ),
            fillcolor="lightblue",
            line=dict(width=0),
            opacity=0.2,
            name=f"HCO₃ = {hco3_val}%",
            showlegend=False,
            hoverinfo="skip"  # kein Hover
        ))

        # Beschriftung im Plot
        fig.add_annotation(
            x=(x_min + x_max) / 2, y=-3,
            text=f"<b>HCO₃ = {hco3_val}%</b>",
            showarrow=False,
            font=dict(size=12, color="blue"),
            xanchor="center",
            yanchor="bottom"
        )


    # === Theoretischer Balancepunkt berechnen ===
    x_theoretical = custom_transform_optimal(50000000)
    y_theoretical = custom_transform_optimal(50000000)

    x_theoretical_scaled = x_theoretical / df["Anionen_trans_raw"].max() * 100
    y_theoretical_scaled = y_theoretical / df["Kationen_trans_raw"].max() * 100

    print(f"📍 Theoretischer 50|50 Punkt: x={x_theoretical_scaled:.2f}, y={y_theoretical_scaled:.2f}")

    # Achsenreichweite so erweitern, dass der Punkt sichtbar ist (mit 5% Puffer)
    xmax = max(100, x_theoretical_scaled * 1.08)
    ymax = max(100, y_theoretical_scaled * 1.08)

    # === Layout ===
    # === Layout ===
    fig.update_layout(
        xaxis=dict(
            title=dict(text="", font=dict(size=20)),
            tickvals=[0, 100],
            ticktext=["", f"HCO₃ (≈ {hco3_max}%)"],
            tickfont=dict(size=14),
            showline=False,  # ❌ schwarze Achsenlinie ausschalten
            zeroline=False,
            range=[0, xmax]
        ),
        yaxis=dict(
            title=dict(text="", font=dict(size=20)),
            tickvals=[0, 100],
            ticktext=["", f"Ca (≈ {ca_max}%)"],
            tickfont=dict(size=14),
            tickangle=-90,
            showline=False,  # ❌ schwarze Achsenlinie ausschalten
            zeroline=False,
            range=[-3, ymax]
        ),
        title=dict(
            text="",
            font=dict(size=24), x=0.5, xanchor="center"
        ),
        legend=dict(
            font=dict(
                size=28,  # 🔼 größer
                color="black",
                family="Arial Black"  # 🔥 fett wie Labels
            ),

            itemsizing="trace",

            x=0.97,
            y=0.9,
            xanchor="left",
            yanchor="top",

            bgcolor="rgba(255,255,255,0.9)",  # 🔼 klarer
            bordercolor="black",
            borderwidth=2  # 🔥 kräftiger Rahmen
        ),
        hoverlabel=dict(font_size=20),
        margin=dict(l=0, r=80, t=60, b=40),
        plot_bgcolor="white"
    )
    fig.add_shape(
        type="rect",
        xref="x",
        yref="y",
        x0=0,
        x1=xmax,
        y0=0,
        y1=ymax,
        fillcolor="rgba(240,245,250,1)",
        line=dict(width=0),
        layer="below"
    )


    # X-Achse (HCO3)
    fig.add_annotation(
        x=100, y=0, ax=0, ay=0,  # statt x=hco3_max → x=100
        xref="x", yref="y", axref="x", ayref="y",
        showarrow=True, arrowhead=2, arrowsize=1.5,
        arrowwidth=1.5, arrowcolor="black", text=""
    )

    # Y-Achse (Ca)
    fig.add_annotation(
        x=0, y=100, ax=0, ay=0,  # statt y=ca_max → y=100
        xref="x", yref="y", axref="x", ayref="y",
        showarrow=True, arrowhead=2, arrowsize=1.5,
        arrowwidth=1.5, arrowcolor="black", text=""
    )

    # Vertikale Linie bei HCO₃ max
    fig.add_shape(
        type="line",
        x0=100, x1=100,
        y0=0, y1=100,
        xref="x", yref="y",
        line=dict(color="black", width=0.5, dash="dash")
    )

    # Horizontale Linie bei Ca max
    fig.add_shape(
        type="line",
        x0=0, x1=100,
        y0=100, y1=100,
        xref="x", yref="y",
        line=dict(color="black", width=0.5, dash="dash")
    )

        # ============================================================
    # 🎨 NICHTLINEARE COLORBAR (0–4 gestreckt)
    # ============================================================

    max_maha = df["LogEuclid"].max()
    if pd.isna(max_maha) or max_maha <= 0:
        max_maha = 1.0

    t = 1.2 / max_maha
    gamma = 0.5

    def stretch(x):
        return min((x ** gamma) * t, 1.0)

    custom_scale = [
        [0.0, "rgb(49,54,149)"],
        [stretch(0.01), "rgb(55,70,160)"],
        [stretch(0.02), "rgb(60,90,170)"],
        [stretch(0.03), "rgb(65,105,175)"],
        [stretch(0.05), "rgb(69,117,180)"],
        [stretch(0.07), "rgb(80,130,190)"],
        [stretch(0.10), "rgb(100,150,205)"],
        [stretch(0.15), "rgb(120,170,215)"],
        [stretch(0.20), "rgb(140,190,225)"],
        [stretch(0.25), "rgb(160,210,230)"],
        [stretch(0.30), "rgb(180,225,235)"],
        [stretch(0.35), "rgb(200,235,240)"],
        [stretch(0.40), "rgb(215,240,245)"],
        [stretch(0.45), "rgb(224,243,248)"],
        [min(t, 1.0), "rgb(255,255,191)"],
        [min(t + (1 - t) * 0.2, 1.0), "rgb(253,174,97)"],
        [min(t + (1 - t) * 0.5, 1.0), "rgb(244,109,67)"],
        [1.0, "rgb(165,0,38)"]
    ]

    custom_scale = sorted(custom_scale, key=lambda z: z[0])

    # ============================================================
    # 🎯 PUNKTE MIT LOG-EUCLIDEAN-FARBEN
    # ============================================================

    df["Group_clean"] = df["Art"].astype(str).str.strip().str.lower()
    df["LogEuclid"] = df["Group_clean"].map(mah_dict)

    print("Punkte gesamt:", len(df))
    print("LogEuclid gültig:", df["LogEuclid"].notna().sum())
    print("LogEuclid NaN:", df["LogEuclid"].isna().sum())

    if df["LogEuclid"].isna().any():
        print("Nicht gematchte Gruppen:")
        print(
            df.loc[df["LogEuclid"].isna(), "Art"]
            .drop_duplicates()
            .head(30)
        )

    df["LogEuclid"] = df["LogEuclid"].fillna(0)

    art_order = (
        df.groupby("Art")["LogEuclid"]
        .median()
        .sort_values(ascending=False)
        .index
    )

    for i, art in enumerate(art_order):

        sub = df[df["Art"] == art]

        if sub.empty:
            continue

        art_str = str(art).upper()

        if art_str.startswith("DA"):
            symbol_shape = "triangle-up"
            marker_size = 12
        elif art_str.startswith("GW"):
            symbol_shape = "square"
            marker_size = 10
        elif art_str.startswith("FW"):
            symbol_shape = "star"
            marker_size = 11
        else:
            symbol_shape = "circle"
            marker_size = 10

        fig.add_trace(go.Scatter(
            x=sub["Anionen_trans"],
            y=sub["Kationen_trans"],
            mode="markers",
            name=art,
            marker=dict(
                symbol=symbol_shape,
                size=marker_size,
                color=sub["LogEuclid"],
                colorscale=custom_scale,
                cmin=0,
                cmax=max_maha,
                showscale=(i == 0),
                colorbar=dict(
                    title=dict(
                        text="Log-Euclidean Distance<br>(to Hallstatt)",
                        font=dict(size=12, family="Arial Black", color="black")
                    ),
                    tickfont=dict(size=10),
                    tickvals=[0, 1, 2, 3, 4, round(max_maha, 1)],
                    ticktext=["0", "1", "2", "3", "4", f"{max_maha:.1f}"],
                    x=0.12,
                    y=0.5,
                    xanchor="right",
                    yanchor="middle",
                    len=1,
                    thickness=24
                ),
                line=dict(width=0.5, color="black")
            ),
            text=sub["hover_text"],
            hoverinfo="text"
        ))

              

        # Überlappungen (Ringe)
        overlaps = df[df["Symbol"] == "star"].copy()
        if not overlaps.empty:
            base_size = 6
            ring_width = 4

            grouped = overlaps.groupby(["Kationen_trans", "Anionen_trans"])

            for (y0, x0), g in grouped:
                arts = list(g["Art"])
                n = len(arts)

                # zentrales X
                fig.add_trace(go.Scatter(
                    x=[x0], y=[y0],
                    mode="markers",
                    marker=dict(
                        symbol="x",
                        size=8,
                        color="red",
                        line=dict(width=3, color="darkred")
                    ),
                    text=[f"Overlap with {n} groups"],
                    hoverinfo="text",
                    showlegend=False
                ))

                # konzentrische rote Halos
                for i, art in enumerate(arts):
                    row = g[g["Art"] == art].iloc[0]
                    size = base_size + i * ring_width

                    fig.add_trace(go.Scatter(
                        x=[x0], y=[y0],
                        mode="markers",
                        marker=dict(
                            symbol="circle",
                            size=size,
                            color="rgba(0,0,0,0)",  # transparent innen
                            line=dict(
                                width=5,
                                color="red"
                            ),
                        ),
                        text=[row["hover_text"]],
                        hoverinfo="text",
                        showlegend=False
                    ))

    # ============================================================
    # 📍 ZENTRALE PUNKTE DER SUBGRUPPEN
    # ============================================================
    group_centers = (
        df.groupby("Art")[["Anionen_trans", "Kationen_trans"]]
        .median()
    )

    from scipy.spatial.distance import pdist, squareform
    from scipy.stats import pearsonr, spearmanr

    center_dist = pd.DataFrame(
        squareform(pdist(group_centers.values, metric="euclidean")),
        index=group_centers.index,
        columns=group_centers.index
    )

    print("\n📏 Distanzmatrix der Plot-Zentren:")
    print(center_dist.round(2))

    # ============================================================
    # 🔗 Korrelation Plotdistanz vs LED
    # ============================================================

    common_groups = [g for g in center_dist.index if g in group_means.index]

    plot_vals = []
    led_vals = []

    for i in range(len(common_groups)):
        for j in range(i + 1, len(common_groups)):
            g1 = common_groups[i]
            g2 = common_groups[j]

            plot_vals.append(center_dist.loc[g1, g2])

            led = log_euclid(
                group_means.loc[g1].values,
                group_means.loc[g2].values
            )
            led_vals.append(led)

    pear_r, pear_p = pearsonr(plot_vals, led_vals)
    spear_r, spear_p = spearmanr(plot_vals, led_vals)

    print("\n🔗 Korrelation Plotdistanz vs LED")
    print(f"Pearson r  = {pear_r:.3f}  (p={pear_p:.4f})")
    print(f"Spearman ρ = {spear_r:.3f}  (p={spear_p:.4f})")

    # ============================================================
    # 🔷 CONVEX HULL PRO SUBGRUPPE
    # ============================================================

    for art in df["Art"].unique():

        sub = df[df["Art"] == art]

        # nur sinnvoll wenn genug Punkte
        if len(sub) < 3:
            continue

        points = sub[["Anionen_trans", "Kationen_trans"]].values

        try:
            hull = ConvexHull(points)

            hull_points = points[hull.vertices]

            # schließen der Linie
            hull_points = np.vstack([hull_points, hull_points[0]])

            fig.add_trace(go.Scatter(
                x=hull_points[:, 0],
                y=hull_points[:, 1],
                mode="lines",
                line=dict(
                    width=1.5,
                    color="rgba(0,0,0,0.8)"  # dünn + leicht transparent
                ),
                showlegend=False,
                hoverinfo="skip"
            ))

        except:
            pass  # falls numerische Probleme



    # === Referenzpunkt (Balancepunkt 50|50) hinzufügen ===
    # → Berechne reale transformierte Koordinaten für Metazahl 50000000
    x_theoretical = custom_transform_optimal(50000000)
    y_theoretical = custom_transform_optimal(50000000)

    # In dieselbe Skala wie die anderen Punkte (0–100 relativ zu Daten-Maximum)
    x_theoretical_scaled = x_theoretical / df["Anionen_trans_raw"].max() * 100
    y_theoretical_scaled = y_theoretical / df["Kationen_trans_raw"].max() * 100



    # === Diagonale Linie vom Ursprung (0,0) zum theoretischen Gleichgewichtspunkt ===
    # === Diagonale Linie vom Ursprung (0,0) zum theoretischen Gleichgewichtspunkt ===

    # Steigung der ursprünglichen Balance-Linie
    slope = y_theoretical_scaled / x_theoretical_scaled

    # Datenbereich bestimmen
    x_data_max = df["Anionen_trans"].max()
    y_data_max = df["Kationen_trans"].max()

    # Schnittpunkt der Linie mit dem Datenbereich berechnen
    y_at_xmax = slope * x_data_max

    if y_at_xmax <= y_data_max:
        x_end = x_data_max
        y_end = y_at_xmax
    else:
        y_end = y_data_max
        x_end = y_data_max / slope

    # Linie zeichnen (gekürzt auf Datenbereich)
    fig.add_shape(
        type="line",
        x0=0, y0=0,
        x1=x_end, y1=y_end,
        xref="x", yref="y",
        line=dict(color="grey", width=0.5, dash="dot"),
    )

    # Beschriftung mittig auf der gekürzten Linie
    fig.add_annotation(
        x=x_end * 0.5,
        y=y_end * 0.5,
        text="dotted line = Ca - HCO3 ~ 1:1 ",
        showarrow=False,
        font=dict(size=12, color="blue"),
        bgcolor="white",
        opacity=0.8
    )

    # === Overlap-Statistik & überlappende Gruppen (oben links) ===
    # === Overlap-Statistik & eindeutige überlappende Gruppen (oben links) ===
    from itertools import combinations

    # Overlap-Koordinaten zählen
    koord_counts = (
        df.groupby(["Kationen_trans", "Anionen_trans"])["Art"]
        .nunique()
        .reset_index(name="region_count")
    )
    overlap_coords_df = koord_counts[koord_counts["region_count"] > 1]
    n_overlap_coords = overlap_coords_df.shape[0]
    total_coords = koord_counts.shape[0]

    overlap_points = int(df["Symbol"].eq("star").sum())
    total_points = len(df)
    pct_overlap_points = (overlap_points / total_points * 100) if total_points else 0
    pct_overlap_coords = (n_overlap_coords / total_coords * 100) if total_coords else 0
    avg_arts_per_overlap = (
        float(overlap_coords_df["region_count"].mean()) if n_overlap_coords else 0.0
    )

    # Eindeutige Paare sammeln (jedes nur einmal)
    pair_set = set()
    for (_, _), g in df.groupby(["Kationen_trans", "Anionen_trans"]):
        arts_here = sorted(g["Art"].unique())
        if len(arts_here) > 1:
            for a, b in combinations(arts_here, 2):
                pair_set.add(f"{a} × {b}")

    # Text für Box
    if pair_set:
        overlap_text = (
                f"<span style='font-size:15px;'><b>Overlap statistic</b></span><br>"
                f"Points in overlaps: {overlap_points} / {total_points} ({pct_overlap_points:.1f}%)<br>"
                f"Coordinates with overlaps: {n_overlap_coords} / {total_coords} ({pct_overlap_coords:.1f}%)<br>"
                f"Ø Types per overlap coordinate: {avg_arts_per_overlap:.2f}<br>"
                f"<b>Overlapping groups:</b><br>"
                + "<br>".join(f"• {p}" for p in sorted(pair_set))
        )
    else:
        overlap_text = (
            f"<b>Overlap statistic</b><br>"
            f"No overlapping determined"
        )




    # Box oben links einfügen
    fig.add_annotation(
        xref="paper", yref="paper",
        x=0.28, y=0.98,  # 🔼 höher & zentriert
        xanchor="center", yanchor="top",
        text=overlap_text,
        showarrow=False,
        font=dict(size=16),
        align="left",
        bgcolor="rgba(255,255,255,0.95)",
        bordercolor="black",
        borderwidth=1.5,
        width=400  # 🔥 macht die Box breit!
    )


    def smart_label(name):

        name = name.replace("_", " ").strip()

        parts = name.split()

        # --- Lake ---
        if parts[0].lower() == "lake":
            if len(parts) > 1:
                return f"La {parts[1][:2].capitalize()}"
            return "La"

        # --- GW / DA / FW ---
        prefix = parts[0].upper()

        if prefix in ["GW", "DA", "FW"]:
            if len(parts) > 1:
                second = parts[1]

                # zusammengesetzte Namen kürzen
                second = second.replace("-", " ")
                subparts = second.split()

                if len(subparts) >= 2:
                    return f"{prefix} {subparts[0][:2].capitalize()}-{subparts[1][:2].capitalize()}"
                else:
                    return f"{prefix} {subparts[0][:3].capitalize()}"

            return prefix

        # fallback
        return name[:6]


    # ============================================================
    # ============================================================
    # 🏷️ LABEL COLLISION AVOIDANCE
    # ============================================================

    placed_labels = []


    def move_if_overlap(x, y, min_dx=6, min_dy=3):
        offsets = [
            (0, 0),
            (0, 5),
            (0, -5),
            (6, 0),
            (-6, 0),
            (6, 5),
            (-6, 5),
            (6, -5),
            (-6, -5),
        ]

        for dx, dy in offsets:
            new_x = x + dx
            new_y = y + dy

            overlap = False
            for px, py in placed_labels:
                if abs(new_x - px) < min_dx and abs(new_y - py) < min_dy:
                    overlap = True
                    break

            if not overlap:
                placed_labels.append((new_x, new_y))
                return new_x, new_y

        placed_labels.append((x, y))
        return x, y


    for art in df["Art"].unique():

        sub = df[df["Art"] == art]
        if sub.empty:
            continue

        x_center = sub["Anionen_trans"].median()
        y_center = sub["Kationen_trans"].median()


        # verschobene Position falls nötig
        x_lab, y_lab = move_if_overlap(x_center, y_center)

        label = smart_label(art)

        fig.add_annotation(
            x=x_lab,
            y=y_lab,
            text=label,
            showarrow=False,  # ← KEINE Pfeile
            font=dict(
                size=14,
                color="black",
                family="Arial Black"
            ),
            bgcolor="rgba(255,255,255,0.8)",
            bordercolor="black",
            borderwidth=1,
            xanchor="center",
            yanchor="middle"
        )


    print("Varianzen:")
    print(np.var(raw_df[ion_cols], axis=0))

    print("\nKorrelationsmatrix:")
    print(np.corrcoef(raw_df[ion_cols].values.T))
    # Export & Show
    # Feste Plotgröße wie im HTML/CMD-Output
       # Export & Show

    fig.update_layout(
        height=900,
        autosize=True,

        margin=dict(
            l=45,
            r=20,
            t=100,
            b=70
        ),

        xaxis=dict(
            domain=[0.01, 0.99],
            title=dict(text="", font=dict(size=18)),
            tickvals=[0, 100],
            ticktext=["", f"HCO₃ (≈ {hco3_max}%)"],
            tickfont=dict(size=18),
            showline=False,
            zeroline=False,
            range=[0, xmax]
        ),

        yaxis=dict(
            title=dict(text="", font=dict(size=18)),
            tickvals=[0, 100],
            ticktext=["", f"Ca (≈ {ca_max}%)"],
            tickfont=dict(size=18),
            tickangle=-90,
            showline=False,
            zeroline=False,
            range=[-3, ymax]
        ),

        legend=dict(
            x=1.02,
            y=0.98,
            xanchor="left",
            yanchor="top",
            font=dict(size=14, color="black", family="Arial"),
            bgcolor="rgba(255,255,255,0.95)",
            bordercolor="black",
            borderwidth=1
        ),

        hoverlabel=dict(font_size=16),
        plot_bgcolor="white"
    )

    html = fig.to_html(
        include_plotlyjs="cdn",
        full_html=False,
        config={"responsive": False}
    )

    components.html(
        html,
        height=750,
        scrolling=True
    )
    # Ergebnisse (Grenzen) auch ausgeben
    print("\nCa-Grenzen aus Daten:")
    for r in results_ca:
        print(f"Ca={r['Ca']}%  ->  y_min={r['y_min']:.2f}  y_max={r['y_max']:.2f}")

    print("\nHCO3-Grenzen aus Daten:")
    for r in results_hco3:
        print(f"HCO3={r['HCO3']}%  ->  x_min={r['x_min']:.2f}  x_max={r['x_max']:.2f}")



except Exception as e:
    print("❌ Fehler beim Plotten:", e)

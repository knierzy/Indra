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

from pathlib import Path

# Input / output paths for GitHub repository use

BASE_DIR = Path(__file__).resolve().parent

DATA_DIR = BASE_DIR / "data"
OUTPUT_DIR = BASE_DIR / "outputs"

OUTPUT_DIR.mkdir(exist_ok=True)

input_file = DATA_DIR / "compendium.xlsx"

preferred_sheet = "Sheet1"

output_file = OUTPUT_DIR / "compendium_processed.xlsx"

output_file_cartesian = (
    OUTPUT_DIR / "CartesianProduct_constraints.xlsx"
)
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
            probe.iloc[i].astype(str).tolist()
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
    'SAMPLING_DATE':       pick(cols, r'sampling|sampling[-\s]?date|date'),
    'ACID_NEUTRALIZING_CAPACITY': pick(cols, r'acid\s*neutralizing|anc|sbv'),

    'CALCIUM_mg_L':        pick(cols, r'calcium'),
    'MAGNESIUM_mg_L':      pick(cols, r'magnesium'),
    'SODIUM_mg_L':         pick(cols, r'sodium'),
    'POTASSIUM_mg_L':      pick(cols, r'potassium'),
    'NITRATE_mg_L':        pick(cols, r'nitrate'),
    'CHLORIDE_mg_L':       pick(cols, r'chloride'),
    'SULFATE_mg_L':        pick(cols, r'sulfate'),
    'BICARBONATE_mg_L':    pick(cols, r'bicarbonate|hydrogencarbonate|hco3'),

    'pH':                  pick(cols, r'\bph\b'),
}


# Create Art column from column B

df['Art'] = df.iloc[:, 1].astype(str).str.strip()
print("✅ Group column created.")


print("\n🔎 Group size by Art:")
print(df.groupby('Art').size().describe())



print("\n🔎 Group sizes for percentile filtering:")


print("\n🔎 Column mapping:")
for k, v in mapping.items():
    print(f"  {k:15s} → {v if (isinstance(v, str) and v in df.columns) else str(v)}")



# Use existing HCO3 values if available

if hco3_col:
    df['HCO3_mg_L_original'] = df[hco3_col].apply(to_num)
else:
    df['HCO3_mg_L_original'] = np.nan


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

df['Ca_mg_L'] = df[mapping['CALCIUM mg/l']].apply(to_num) if mapping['CALCIUM mg/l'] else np.nan
df['Mg_mg_L'] = df[mapping['MAGNESIUM mg/l']].apply(to_num) if mapping['MAGNESIUM mg/l'] else np.nan
df['Na_mg_L'] = df[mapping['NATRIUM mg/l']].apply(to_num) if mapping['NATRIUM mg/l'] else np.nan
df['K_mg_L'] = df[mapping['KALIUM mg/l']].apply(to_num) if mapping['KALIUM mg/l'] else np.nan
df['Cl_mg_L'] = df[mapping['CHLORID mg/l']].apply(to_num) if mapping['CHLORID mg/l'] else np.nan
df['SO4_mg_L'] = df[mapping['SULFAT mg/l']].apply(to_num) if mapping['SULFAT mg/l'] else np.nan

no3_raw = df[mapping['NITRAT-N mg/l']].apply(to_num) if mapping['NITRAT-N mg/l'] else np.nan


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

print("\n🔎 HCO3 sources:")
print("HCO3 original:", df['HCO3_mg_L_original'].notna().sum())
print("SBV:", df['SBV_mmol_L'].notna().sum())
print("HCO3 quick:", df['HCO3_mg_L_quick'].notna().sum())
print("HCO3 final:", df['HCO3_mg_L_final'].notna().sum())






n_before = len(df)

df = df.dropna(subset=required_ions).copy()

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

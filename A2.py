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
# === Datei & Zielpfad ===
from scipy.spatial.distance import euclidean

def log_euclid(a, b):
    return euclidean(np.log1p(a), np.log1p(b))



from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
OUTPUT_DIR = BASE_DIR / "outputs"

OUTPUT_DIR.mkdir(exist_ok=True)

input_file = DATA_DIR / "CartesianProduct_constraints.xlsx"
raw_file = DATA_DIR / "compendium_processed.xlsx"

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
def custom_transform_optimal(x, base=math.e +14): #12.1415926535
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
    df = pd.read_excel(input_file)

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

    h_name = [g for g in group_means.index if "hall" in g.lower()][0]
    o_name = [g for g in group_means.index if "oss" in g.lower()][0]

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
    ref_candidates = [g for g in group_means_log.index if "hall" in g.lower()]
    if not ref_candidates:
        raise ValueError("❌ Keine Hallstatt-Gruppe gefunden!")

    ref_group = ref_candidates[0]
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
    
    # ============================================================
    # Mapping Plotgruppen → Referenzgruppen
    # ============================================================

    def match_maha(name):

        name = str(name).strip().lower()

        mapping = {
            "da_altheim": "tgw_altheim",
            "da_bad schallerbach": "tgw_bad schallerbach",
            "da_buch-st. magdalena": "tgw_buch-st. magdalena",
            "da_großwilfersdorf": "tgw_großwilfersdorf",
            "da_rottenbach": "tgw_rottenbach",
            "da_senftenbach": "tgw_senftenbach",

            "gw_gaweinstal": "gaweinstal_pg31600452",
            "gw_groß-enzersdorf": "groß-enzersdorf_pg30800302",
            "gw_laa_an_der_thaya": "laa_pg31600422",
            "gw_mureck": "mureck_pg61511062",
            "gw_traiskirchen": "traiskirchen_pg30600152",

            "fw_tux": "kk72410012_tux",

            "lake constance": "bodensee",
            "lake fuschl": "fuschlsee",
            "lake hallstatt": "hallstätter see",
            "lake millstatt": "millstätter see",
            "lake neusiedl": "neusiedlersee",
            "lake ossiach": "ossiacher see",
            "lake wolfgang": "wolfgangsee"
        }

        # 1️⃣ direkte Zuordnung
        if name in mapping:
            return mah_dict.get(mapping[name], np.nan)

        # 2️⃣ exakter Match
        if name in mah_dict:
            return mah_dict[name]

        # 3️⃣ unscharfer Match
        for key in mah_dict.keys():

            key_norm = str(key).strip().lower()

            if name in key_norm:
                return mah_dict[key]

            if key_norm in name:
                return mah_dict[key]

        return np.nan

    fig = go.Figure()

    mah_sorted = sorted(mah_dict.items(), key=lambda x: x[1])

    mah_text = "<span style='font-size:24px'><b>Log-Euclidan distance – reference: Hallstatt</b></span><br>"

    for g, d in mah_sorted[:5]:  # Top 5
        mah_text += f"{g.title()}: {d:.2f}<br>"

    fig.add_annotation(
        xref="paper", yref="paper",
        x=0.55, y=1.05,
        xanchor="left", yanchor="top",
        text=mah_text,
        showarrow=False,
        font=dict(size=24),  # Basisgröße für den Rest
        align="left",
        bgcolor="rgba(255,255,255,0.95)",
        bordercolor="black", borderwidth=1.5
    )

    print("\n📏 Mahalanobis-Distanzen relativ zu Hallstatt:\n")

    for g, d in sorted(mah_dict.items(), key=lambda x: x[1]):
        print(f"{g:25s}  →  {d:.3f}")

    df["LogEuclid"] = df["Group_clean"].apply(match_maha)

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
                solidity=0.2
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
            font=dict(size=20, color="grey"),
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
                solidity=0.2
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
            font=dict(size=20, color="blue"),
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
            tickfont=dict(size=22),
            showline=False,  # ❌ schwarze Achsenlinie ausschalten
            zeroline=False,
            range=[0, xmax]
        ),
        yaxis=dict(
            title=dict(text="", font=dict(size=20)),
            tickvals=[0, 100],
            ticktext=["", f"Ca (≈ {ca_max}%)"],
            tickfont=dict(size=22),
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
    ),


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
        line=dict(color="black", width=1.5, dash="dash")
    )

    # Horizontale Linie bei Ca max
    fig.add_shape(
        type="line",
        x0=0, x1=100,
        y0=100, y1=100,
        xref="x", yref="y",
        line=dict(color="black", width=1.5, dash="dash")
    )

    # ============================================================
    # 🎨 NICHTLINEARE COLORBAR (0–4 gestreckt)
    # ============================================================

    t = 1.2 / max_maha if max_maha > 0 else 0.5
    gamma = 0.5  # 🔥 Stärke der Verzerrung (0.3 = sehr stark, 0.6 = moderat)


    def stretch(x):
        return (x ** gamma) * t


    custom_scale = [
        [0.0, "rgb(49,54,149)"],

        # 🔥 EXTREM fein 0–1
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

        [t, "rgb(255,255,191)"],

        # 🔽 stark komprimiert oben
        [t + (1 - t) * 0.2, "rgb(253,174,97)"],
        [t + (1 - t) * 0.5, "rgb(244,109,67)"],
        [1.0, "rgb(165,0,38)"]
    ]

    # ============================================================
    # 🎯 PUNKTE MIT MAHALANOBIS-FARBEN
    # ============================================================

    # 🔥 Median-Mahalanobis pro Art berechnen
    art_order = (
        df.groupby("Art")["LogEuclid"]
        .median()
        .sort_values(ascending=False)
        .index
    )

    # 🔥 danach plotten
    for i, art in enumerate(art_order):

        sub = df[df["Art"] == art]

        if sub.empty:
            continue

        art_str = str(art).upper()

        # Symbol-Logik
        if art_str.startswith("DA"):
            symbol_shape = "triangle-up"
            marker_size = 26
        elif art_str.startswith("GW"):
            symbol_shape = "square"
            marker_size = 22
        elif art_str.startswith("FW"):
            symbol_shape = "star"
            marker_size = 26
        else:
            symbol_shape = "circle"
            marker_size = 24

        fig.add_trace(go.Scatter(
            x=sub["Anionen_trans"],
            y=sub["Kationen_trans"],
            mode="markers",
            name=art,
            marker=dict(
                symbol=symbol_shape,
                size=marker_size,

                # 🔥 ORIGINALWERTE (kein sqrt!)
                color=sub["LogEuclid"],
                colorscale=custom_scale,

                cmin=0,
                cmax=max_maha,  # 🔥 wieder korrekt

                showscale=(i == 0),

                colorbar=dict(
                    title=dict(
                        text="Log-Euclidean Distance<br>(to Hallstatt)",
                        font=dict(
                            size=22,
                            family="Arial Black",
                            color="black"
                        )
                    ),

                    tickfont=dict(
                        size=22
                    ),

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
            base_size = 18
            ring_width = 8

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
                        size=12,
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
                    width=3.5,
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
        line=dict(color="grey", width=3.5, dash="dot"),
    )

    # Beschriftung mittig auf der gekürzten Linie
    fig.add_annotation(
        x=x_end * 0.5,
        y=y_end * 0.5,
        text="dotted line = Ca - HCO3 ~ 1:1 ",
        showarrow=False,
        font=dict(size=20, color="blue"),
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
                f"<span style='font-size:22px;'><b>Overlap statistic</b></span><br>"
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
        x=0.32, y=1.05,  # 🔼 höher & zentriert
        xanchor="center", yanchor="top",
        text=overlap_text,
        showarrow=False,
        font=dict(size=24),
        align="left",
        bgcolor="rgba(255,255,255,0.95)",
        bordercolor="black",
        borderwidth=1.5,
        width=600  # 🔥 macht die Box breit!
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
                size=25,
                color="black",
                family="Arial Black"
            ),
            bgcolor="rgba(255,255,255,0.8)",
            bordercolor="black",
            borderwidth=1,
            xanchor="center",
            yanchor="middle"
        )





    import numpy as np

    print("Varianzen:")
    print(np.var(raw_df[ion_cols], axis=0))

    print("\nKorrelationsmatrix:")
    print(np.corrcoef(raw_df[ion_cols].values.T))

    # Export & Show
    fig.write_html(plot_output)
    print(f"\n✅ Plot gespeichert unter:\n→ {plot_output}")
    fig.show()

    # Ergebnisse (Grenzen) auch ausgeben
    print("\nCa-Grenzen aus Daten:")
    for r in results_ca:
        print(f"Ca={r['Ca']}%  ->  y_min={r['y_min']:.2f}  y_max={r['y_max']:.2f}")

    print("\nHCO3-Grenzen aus Daten:")
    for r in results_hco3:
        print(f"HCO3={r['HCO3']}%  ->  x_min={r['x_min']:.2f}  x_max={r['x_max']:.2f}")



except Exception as e:
    print("❌ Fehler beim Plotten:", e)

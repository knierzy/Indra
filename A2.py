# This script generates an INDRA Projection (IP) plot for multivariate
# hydrochemical datasets based on transformed cation and anion metanumbers.
# Constrained compositional fields are transformed using an irrational base and 
# visualized in 2D. The script displays subgroup fields, theoretical Ca and HCO3
# reference bands, and compositional overlaps between subgroups. Log-Euclidean 
# distances are calculated from the original eight-dimensional
# hydrochemical compositions and used to assess compositional similarity 
# between subgroups and a selected reference group.

import math
import pandas as pd
import plotly.graph_objects as go
from plotly.express.colors import qualitative
import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial import distance_matrix

# === File and output paths ===
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
# --- Helper function for transformation using a freely selectable base ---
def custom_transform_optimal(x, base=math.e + 14):  # Example base: e + 14
    try:
        x_str = str(int(x)).zfill(8)
        a = int(x_str[0:2])
        b = int(x_str[2:4])
        c = int(x_str[4:6])
        d = int(x_str[6:8])
        return a * base**3 + b * base**2 + c * base + d
    except:
        return None

# --- Helper functions for hover information ---
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

    k_perc, _ = pairs_to_percentages(
        row["Metazahl_Kationen"], k_labels
    )
    a_perc, _ = pairs_to_percentages(
        row["Metazahl_Anionen"], a_labels
    )

    k_lines = " · ".join([
        f"{lbl}: {k_perc[lbl]}%"
        if k_perc[lbl] is not None
        else f"{lbl}: –"
        for lbl in k_labels
    ])

    a_lines = " · ".join([
        f"{lbl}: {a_perc[lbl]}%"
        if a_perc[lbl] is not None
        else f"{lbl}: –"
        for lbl in a_labels
    ])

    return (
        f"<b>Group:</b> {row['Art']}<br>"
        f"<b>Cations</b> (from {str(row['Metazahl_Kationen']).zfill(8)}):<br>"
        f"{k_lines}<br>"
        f"<b>Anions</b> (from {str(row['Metazahl_Anionen']).zfill(8)}):<br>"
        f"{a_lines}"
    )


try:
    # Read Excel input file
    df = pd.read_excel(input_file)

    required_cols = [
        'Metazahl_Kationen',
        'Metazahl_Anionen',
        'Art'
    ]

    for col in required_cols:
        if col not in df.columns:
            raise ValueError(
                f"Required column '{col}' is missing from the input file."
            )

    # Apply metanumber transformation
    df["Kationen_trans_raw"] = (
        df["Metazahl_Kationen"]
        .apply(custom_transform_optimal)
    )

    df["Anionen_trans_raw"] = (
        df["Metazahl_Anionen"]
        .apply(custom_transform_optimal)
    )

    # Normalize transformed coordinates to a 0–100 scale
    df["Kationen_trans"] = (
        df["Kationen_trans_raw"]
        / df["Kationen_trans_raw"].max()
        * 100
    )

    df["Anionen_trans"] = (
        df["Anionen_trans_raw"]
        / df["Anionen_trans_raw"].max()
        * 100
    )

    # Remove duplicate coordinates within the same subgroup
    df = df.drop_duplicates(
        subset=[
            "Kationen_trans",
            "Anionen_trans",
            "Art"
        ]
    )

    # Count coordinate overlaps
    koord_counts = (
        df.groupby([
            "Kationen_trans",
            "Anionen_trans"
        ])
        .size()
        .reset_index(name="region_count")
    )

    df = df.merge(
        koord_counts,
        on=[
            "Kationen_trans",
            "Anionen_trans"
        ],
        how="left"
    )

    df["Symbol"] = df["region_count"].apply(
        lambda x: "star" if x > 1 else "circle"
    )

    # Generate hover information
    df["hover_text"] = df.apply(
        format_hover,
        axis=1
    )
    # ============================================================
    # LOG-EUCLIDEAN DISTANCES
    # ============================================================

    # --- Define major-ion variables ---
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

    # --- Calculate subgroup mean compositions ---
    group_means = raw_df.groupby("Art")[ion_cols].mean()
    group_means.index = group_means.index.astype(str).str.strip()

    # --- Define reference subgroup ---
    ref_group = "Lake Hallstatt"

    if ref_group not in group_means.index:
        raise ValueError(
            "Reference subgroup 'Lake Hallstatt' not found."
        )

    print(f"\nReference subgroup: {ref_group}")

    ref_vector = group_means.loc[ref_group].values

    # --- Calculate Log-Euclidean distances to the reference subgroup ---
    led_dict = {}

    for g in group_means.index:
        vec = group_means.loc[g].values
        led_dict[str(g).strip().lower()] = log_euclid(
            vec,
            ref_vector
        )

    # --- Standardize plot subgroup names for matching ---
    df["Group_clean"] = (
        df["Art"]
        .astype(str)
        .str.strip()
        .str.lower()
    )

    # ============================================================
    # Map plot subgroup names to reference dataset subgroup names
    # ============================================================
    def match_led(name):

        name = str(name).strip().lower()

        # Exact match
        if name in led_dict:
            return led_dict[name]

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

        # Direct mapping
        if name in mapping:
            return led_dict.get(mapping[name], np.nan)

        # Fuzzy matching
        for key in led_dict.keys():

            key_norm = str(key).strip().lower()

            if name in key_norm:
                return led_dict[key]

            if key_norm in name:
                return led_dict[key]

        return np.nan

    fig = go.Figure()
    
    # Sort Log-Euclidean distances
    led_sorted = sorted(led_dict.items(), key=lambda x: x[1])

    led_text = (
        "<span style='font-size:24px'>"
        "<b>Log-Euclidean distance – reference: Hallstatt</b>"
        "</span><br>"
    )

    for g, d in led_sorted[:5]:  # Top 5
        led_text += f"{g.title()}: {d:.2f}<br>"

    fig.add_annotation(
        xref="paper",
        yref="paper",
        x=0.55,
        y=1.05,
        xanchor="left",
        yanchor="top",
        text=led_text,
        showarrow=False,
        font=dict(size=24),
        align="left",
        bgcolor="rgba(255,255,255,0.95)",
        bordercolor="black",
        borderwidth=1.5
    )

    # Assign Log-Euclidean distances to plot subgroups
    df["LogEuclid"] = df["Group_clean"].apply(match_led)

    # Print Log-Euclidean distances relative to the reference subgroup
    print("\nLog-Euclidean distances relative to Hallstatt:\n")

    for g, d in led_sorted:
        print(f"{g:25s}  →  {d:.3f}")

    # Maximum Log-Euclidean distance used for color scaling
    max_led = df["LogEuclid"].max()

    # Identify subgroup names that could not be matched
    missing = df[df["LogEuclid"].isna()]["Group_clean"].unique()

    if len(missing) > 0:
        print("\nUnmatched subgroup names:")
        for m in missing[:20]:
            print(m)

    # Check assigned Log-Euclidean distance values
    print("\nLog-Euclidean distance check:")
    print(df["LogEuclid"].head())
    print("Number of missing values:", df["LogEuclid"].isna().sum())

     # === Calculate Ca and HCO3 percentage values ===
    df["Ca_val"] = df["Metazahl_Kationen"].apply(
        lambda x: pairs_to_percentages(
            x,
            ["Ca", "Mg", "Na", "K"]
        )[0]["Ca"]
    )

    df["HCO3_val"] = df["Metazahl_Anionen"].apply(
        lambda x: pairs_to_percentages(
            x,
            ["HCO₃", "SO₄", "Cl", "NO₃"]
        )[0]["HCO₃"]
    )

    ca_max = df["Ca_val"].max()
    hco3_max = df["HCO3_val"].max()

    results_ca = []

    for ca_val in [2, 5, 10, 15, 20, 25, 30, 35, 40]:
        sub = df[df["Ca_val"] == ca_val]

        if sub.empty:
            continue

        y_min = (
            sub["Kationen_trans_raw"].min()
            / df["Kationen_trans_raw"].max()
            * 100
        )

        y_max = (
            sub["Kationen_trans_raw"].max()
            / df["Kationen_trans_raw"].max()
            * 100
        )

        results_ca.append(
            dict(
                Ca=ca_val,
                y_min=y_min,
                y_max=y_max
            )
        )

        # Add Ca reference band
        fig.add_trace(
            go.Scatter(
                x=[0, 100, 100, 0],
                y=[y_min, y_min, y_max, y_max],
                fill="toself",
                fillpattern=dict(
                    shape="/",  # Diagonal hatching
                    fgcolor="grey",
                    size=6,
                    solidity=0.2
                ),
                fillcolor="lightgrey",
                line=dict(width=0),
                opacity=0.3,
                name=f"Ca = {ca_val}%",
                showlegend=False,
                hoverinfo="skip"  # Disable hover information
            )
        )

        fig.add_annotation(
            x=0,
            y=(y_min + y_max) / 2,
            text=f"<b>Ca = {ca_val}%</b>",
            showarrow=False,
            font=dict(size=20, color="grey"),
            xanchor="left",
            yanchor="middle"
        )

    # === Generate HCO3 reference bands ===
    results_hco3 = []

    for hco3_val in [5, 10, 15, 20, 25, 30, 35, 40, 45]:
        sub = df[df["HCO3_val"] == hco3_val]

        if sub.empty:
            continue

        x_min = (
            sub["Anionen_trans_raw"].min()
            / df["Anionen_trans_raw"].max()
            * 100
        )

        x_max = (
            sub["Anionen_trans_raw"].max()
            / df["Anionen_trans_raw"].max()
            * 100
        )

        results_hco3.append(
            dict(
                HCO3=hco3_val,
                x_min=x_min,
                x_max=x_max
            )
        )

        # Add HCO3 reference band
        fig.add_trace(
            go.Scatter(
                x=[x_min, x_max, x_max, x_min],
                y=[0, 0, 100, 100],
                fill="toself",
                fillpattern=dict(
                    shape="\\",  # Opposite diagonal hatching
                    fgcolor="blue",
                    size=6,
                    solidity=0.2
                ),
                fillcolor="lightblue",
                line=dict(width=0),
                opacity=0.2,
                name=f"HCO₃ = {hco3_val}%",
                showlegend=False,
                hoverinfo="skip"  # Disable hover information
            )
        )

        # Add HCO3 reference-band label
        fig.add_annotation(
            x=(x_min + x_max) / 2,
            y=-3,
            text=f"<b>HCO₃ = {hco3_val}%</b>",
            showarrow=False,
            font=dict(size=20, color="blue"),
            xanchor="center",
            yanchor="bottom"
        )


    # === Calculate theoretical 50|50 balance point ===
    x_theoretical = custom_transform_optimal(50000000)
    y_theoretical = custom_transform_optimal(50000000)

    x_theoretical_scaled = (
        x_theoretical
        / df["Anionen_trans_raw"].max()
        * 100
    )

    y_theoretical_scaled = (
        y_theoretical
        / df["Kationen_trans_raw"].max()
        * 100
    )

    print(
        f"Theoretical 50|50 balance point: "
        f"x={x_theoretical_scaled:.2f}, "
        f"y={y_theoretical_scaled:.2f}"
    )

    # Extend axis ranges to include the theoretical balance point
    xmax = max(100, x_theoretical_scaled * 1.08)
    ymax = max(100, y_theoretical_scaled * 1.08)


    # === Plot layout ===
    fig.update_layout(
        xaxis=dict(
            title=dict(
                text="",
                font=dict(size=20)
            ),
            tickvals=[0, 100],
            ticktext=[
                "",
                f"HCO₃ (≈ {hco3_max}%)"
            ],
            tickfont=dict(size=22),
            showline=False,  # Hide default axis line
            zeroline=False,
            range=[0, xmax]
        ),

        yaxis=dict(
            title=dict(
                text="",
                font=dict(size=20)
            ),
            tickvals=[0, 100],
            ticktext=[
                "",
                f"Ca (≈ {ca_max}%)"
            ],
            tickfont=dict(size=22),
            tickangle=-90,
            showline=False,  # Hide default axis line
            zeroline=False,
            range=[-3, ymax]
        ),

        title=dict(
            text="",
            font=dict(size=24),
            x=0.5,
            xanchor="center"
        ),

        legend=dict(
            font=dict(
                size=28,  # Increase legend font size
                color="black",
                family="Arial Black"  # Match subgroup label style
            ),

            itemsizing="trace",

            x=0.97,
            y=0.9,
            xanchor="left",
            yanchor="top",

            bgcolor="rgba(255,255,255,0.9)",  # Semi-opaque background
            bordercolor="black",
            borderwidth=2  # Emphasize legend frame
        ),

        hoverlabel=dict(
            font_size=20
        ),

        margin=dict(
            l=0,
            r=80,
            t=60,
            b=40
        ),

        plot_bgcolor="white"
    )


    # Add background rectangle to the plotting area
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

      # X-axis arrow (HCO3)
    fig.add_annotation(
        x=100, y=0, ax=0, ay=0,  # Use normalized x-axis maximum
        xref="x", yref="y", axref="x", ayref="y",
        showarrow=True, arrowhead=2, arrowsize=1.5,
        arrowwidth=1.5, arrowcolor="black", text=""
    )

    # Y-axis arrow (Ca)
    fig.add_annotation(
        x=0, y=100, ax=0, ay=0,  # Use normalized y-axis maximum
        xref="x", yref="y", axref="x", ayref="y",
        showarrow=True, arrowhead=2, arrowsize=1.5,
        arrowwidth=1.5, arrowcolor="black", text=""
    )

    # Vertical reference line at normalized HCO3 maximum
    fig.add_shape(
        type="line",
        x0=100, x1=100,
        y0=0, y1=100,
        xref="x", yref="y",
        line=dict(color="black", width=1.5, dash="dash")
    )

    # Horizontal reference line at normalized Ca maximum
    fig.add_shape(
        type="line",
        x0=0, x1=100,
        y0=100, y1=100,
        xref="x", yref="y",
        line=dict(color="black", width=1.5, dash="dash")
    )

    # ============================================================
    # NONLINEAR COLOR SCALE FOR LOG-EUCLIDEAN DISTANCES
    # ============================================================

    t = 1.2 / max_led if max_led > 0 else 0.5

    # Controls nonlinear stretching of the lower distance range
    gamma = 0.5  # Lower values produce stronger stretching


    def stretch(x):
        return (x ** gamma) * t


    custom_scale = [
        [0.0, "rgb(49,54,149)"],

        # Fine color resolution for small Log-Euclidean distances
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

        # Compress the upper distance range
        [t + (1 - t) * 0.2, "rgb(253,174,97)"],
        [t + (1 - t) * 0.5, "rgb(244,109,67)"],
        [1.0, "rgb(165,0,38)"]
    ]
    # ============================================================
    # PLOT POINTS COLORED BY LOG-EUCLIDEAN DISTANCE
    # ============================================================

    # Order subgroups by median Log-Euclidean distance
    art_order = (
        df.groupby("Art")["LogEuclid"]
        .median()
        .sort_values(ascending=False)
        .index
    )

    # Plot subgroup points
    for i, art in enumerate(art_order):

        sub = df[df["Art"] == art]

        if sub.empty:
            continue

        art_str = str(art).upper()

        # Define marker symbols by subgroup type
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

        fig.add_trace(
            go.Scatter(
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
                    cmax=max_led,

                    # Display colorbar only once
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

                        tickvals=[
                            0, 1, 2, 3, 4,
                            round(max_led, 1)
                        ],
                        ticktext=[
                            "0", "1", "2", "3", "4",
                            f"{max_led:.1f}"
                        ],

                        x=0.12,
                        y=0.5,
                        xanchor="right",
                        yanchor="middle",

                        len=1,
                        thickness=24
                    ),

                    line=dict(
                        width=0.5,
                        color="black"
                    )
                ),
                text=sub["hover_text"],
                hoverinfo="text"
            )
        )

        # Identify overlapping coordinates
        overlaps = df[df["Symbol"] == "star"].copy()

        if not overlaps.empty:
            base_size = 18
            ring_width = 8

            grouped = overlaps.groupby(
                ["Kationen_trans", "Anionen_trans"]
            )

            for (y0, x0), g in grouped:
                arts = list(g["Art"])
                n = len(arts)

                # Mark overlap center with an X
                fig.add_trace(
                    go.Scatter(
                        x=[x0],
                        y=[y0],
                        mode="markers",
                        marker=dict(
                            symbol="x",
                            size=12,
                            color="red",
                            line=dict(
                                width=3,
                                color="darkred"
                            )
                        ),
                        text=[f"Overlap with {n} groups"],
                        hoverinfo="text",
                        showlegend=False
                    )
                )

                # Add concentric red rings for overlapping subgroups
                for i, art in enumerate(arts):
                    row = g[g["Art"] == art].iloc[0]
                    size = base_size + i * ring_width

                    fig.add_trace(
                        go.Scatter(
                            x=[x0],
                            y=[y0],
                            mode="markers",
                            marker=dict(
                                symbol="circle",
                                size=size,
                                color="rgba(0,0,0,0)",  # Transparent fill
                                line=dict(
                                    width=5,
                                    color="red"
                                )
                            ),
                            text=[row["hover_text"]],
                            hoverinfo="text",
                            showlegend=False
                        )
                    )

    # ============================================================
    # CALCULATE SUBGROUP CENTERS AND PAIRWISE PLOT DISTANCES
    # ============================================================

    group_centers = (
        df.groupby("Art")[["Anionen_trans", "Kationen_trans"]]
        .median()
    )

    from scipy.spatial.distance import pdist, squareform
    from scipy.stats import pearsonr, spearmanr

    center_dist = pd.DataFrame(
        squareform(
            pdist(
                group_centers.values,
                metric="euclidean"
            )
        ),
        index=group_centers.index,
        columns=group_centers.index
    )

    print("\nEuclidean distance matrix of subgroup centers:")
    print(center_dist.round(2))

    # ============================================================
    # CORRELATION BETWEEN PLOT DISTANCES AND LOG-EUCLIDEAN DISTANCES
    # ============================================================

    common_groups = [
        g for g in center_dist.index
        if g in group_means.index
    ]

    plot_vals = []
    led_vals = []

    for i in range(len(common_groups)):
        for j in range(i + 1, len(common_groups)):
            g1 = common_groups[i]
            g2 = common_groups[j]

            plot_vals.append(
                center_dist.loc[g1, g2]
            )

            led = log_euclid(
                group_means.loc[g1].values,
                group_means.loc[g2].values
            )

            led_vals.append(led)

    pear_r, pear_p = pearsonr(
        plot_vals,
        led_vals
    )

    spear_r, spear_p = spearmanr(
        plot_vals,
        led_vals
    )

    print(
        "\nCorrelation between plot distances "
        "and Log-Euclidean distances"
    )
    print(
        f"Pearson r  = {pear_r:.3f}  "
        f"(p={pear_p:.4f})"
    )
    print(
        f"Spearman ρ = {spear_r:.3f}  "
        f"(p={spear_p:.4f})"
    )

    # ============================================================
    # CONVEX HULL FOR EACH SUBGROUP
    # ============================================================

    for art in df["Art"].unique():

        sub = df[df["Art"] == art]

        # At least three points are required to construct a convex hull
        if len(sub) < 3:
            continue

        points = sub[
            ["Anionen_trans", "Kationen_trans"]
        ].values

        try:
            hull = ConvexHull(points)

            hull_points = points[
                hull.vertices
            ]

            # Close the hull polygon
            hull_points = np.vstack(
                [
                    hull_points,
                    hull_points[0]
                ]
            )

            fig.add_trace(
                go.Scatter(
                    x=hull_points[:, 0],
                    y=hull_points[:, 1],
                    mode="lines",
                    line=dict(
                        width=3.5,
                        color="rgba(0,0,0,0.8)"  # Slightly transparent outline
                    ),
                    showlegend=False,
                    hoverinfo="skip"
                )
            )

        except:
            pass  # Skip subgroups for which a valid hull cannot be computed


    # === Calculate theoretical 50|50 balance point ===
    # Calculate transformed coordinates for metanumber 50000000
    x_theoretical = custom_transform_optimal(50000000)
    y_theoretical = custom_transform_optimal(50000000)

    # Scale to the same 0–100 range as the plotted data
    x_theoretical_scaled = (
        x_theoretical
        / df["Anionen_trans_raw"].max()
        * 100
    )

    y_theoretical_scaled = (
        y_theoretical
        / df["Kationen_trans_raw"].max()
        * 100
    )


    # === Add diagonal Ca–HCO3 balance line ===

    # Calculate slope of the theoretical balance line
    slope = y_theoretical_scaled / x_theoretical_scaled

    # Determine maximum extent of the plotted data
    x_data_max = df["Anionen_trans"].max()
    y_data_max = df["Kationen_trans"].max()

    # Calculate intersection of the balance line with the plotted data range
    y_at_xmax = slope * x_data_max

    if y_at_xmax <= y_data_max:
        x_end = x_data_max
        y_end = y_at_xmax
    else:
        y_end = y_data_max
        x_end = y_data_max / slope

    # Draw balance line truncated to the plotted data range
    fig.add_shape(
        type="line",
        x0=0,
        y0=0,
        x1=x_end,
        y1=y_end,
        xref="x",
        yref="y",
        line=dict(
            color="grey",
            width=3.5,
            dash="dot"
        )
    )

    # Add label at the center of the balance line
    fig.add_annotation(
        x=x_end * 0.5,
        y=y_end * 0.5,
        text="Dotted line = Ca–HCO₃ ≈ 1:1",
        showarrow=False,
        font=dict(
            size=20,
            color="blue"
        ),
        bgcolor="white",
        opacity=0.8
    )


    # ============================================================
    # OVERLAP STATISTICS
    # ============================================================

    from itertools import combinations

    # Count coordinates occupied by more than one subgroup
    koord_counts = (
        df.groupby(
            ["Kationen_trans", "Anionen_trans"]
        )["Art"]
        .nunique()
        .reset_index(name="region_count")
    )

    overlap_coords_df = koord_counts[
        koord_counts["region_count"] > 1
    ]

    n_overlap_coords = overlap_coords_df.shape[0]
    total_coords = koord_counts.shape[0]

    overlap_points = int(
        df["Symbol"].eq("star").sum()
    )

    total_points = len(df)

    pct_overlap_points = (
        overlap_points / total_points * 100
        if total_points
        else 0
    )

    pct_overlap_coords = (
        n_overlap_coords / total_coords * 100
        if total_coords
        else 0
    )

    avg_arts_per_overlap = (
        float(
            overlap_coords_df["region_count"].mean()
        )
        if n_overlap_coords
        else 0.0
    )

    # Collect unique overlapping subgroup pairs
    pair_set = set()

    for (_, _), g in df.groupby(
        ["Kationen_trans", "Anionen_trans"]
    ):
        arts_here = sorted(
            g["Art"].unique()
        )

        if len(arts_here) > 1:
            for a, b in combinations(
                arts_here,
                2
            ):
                pair_set.add(
                    f"{a} × {b}"
                )

    # Generate overlap-statistics text box
    if pair_set:
        overlap_text = (
            f"<span style='font-size:22px;'>"
            f"<b>Overlap statistics</b>"
            f"</span><br>"
            f"Points in overlaps: "
            f"{overlap_points} / {total_points} "
            f"({pct_overlap_points:.1f}%)<br>"
            f"Coordinates with overlaps: "
            f"{n_overlap_coords} / {total_coords} "
            f"({pct_overlap_coords:.1f}%)<br>"
            f"Mean number of groups per overlap coordinate: "
            f"{avg_arts_per_overlap:.2f}<br>"
            f"<b>Overlapping groups:</b><br>"
            + "<br>".join(
                f"• {p}"
                for p in sorted(pair_set)
            )
        )
    else:
        overlap_text = (
            f"<b>Overlap statistics</b><br>"
            f"No overlaps detected"
        )

       # Add overlap-statistics box
    fig.add_annotation(
        xref="paper",
        yref="paper",
        x=0.32,
        y=1.05,  # Position above the main plotting area
        xanchor="center",
        yanchor="top",
        text=overlap_text,
        showarrow=False,
        font=dict(size=24),
        align="left",
        bgcolor="rgba(255,255,255,0.95)",
        bordercolor="black",
        borderwidth=1.5,
        width=600  # Set width of the annotation box
    )


    # Generate abbreviated subgroup labels
    def smart_label(name):

        name = name.replace("_", " ").strip()

        parts = name.split()

        # Lake subgroups
        if parts[0].lower() == "lake":
            if len(parts) > 1:
                return f"La {parts[1][:2].capitalize()}"
            return "La"

        # Groundwater (GW), deep aquifer (DA), and freshwater (FW) subgroups
        prefix = parts[0].upper()

        if prefix in ["GW", "DA", "FW"]:
            if len(parts) > 1:
                second = parts[1]

                # Abbreviate compound subgroup names
                second = second.replace("-", " ")
                subparts = second.split()

                if len(subparts) >= 2:
                    return (
                        f"{prefix} "
                        f"{subparts[0][:2].capitalize()}-"
                        f"{subparts[1][:2].capitalize()}"
                    )
                else:
                    return f"{prefix} {subparts[0][:3].capitalize()}"

            return prefix

        # Fallback for other subgroup names
        return name[:6]


    # ============================================================
    # LABEL COLLISION AVOIDANCE
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
                if (
                    abs(new_x - px) < min_dx
                    and abs(new_y - py) < min_dy
                ):
                    overlap = True
                    break

            if not overlap:
                placed_labels.append((new_x, new_y))
                return new_x, new_y

        placed_labels.append((x, y))
        return x, y


    # Calculate median position of each subgroup
    for art in df["Art"].unique():

        sub = df[df["Art"] == art]

        if sub.empty:
            continue

        x_center = sub["Anionen_trans"].median()
        y_center = sub["Kationen_trans"].median()


        
    # Shift label position if necessary to avoid overlap
    x_lab, y_lab = move_if_overlap(
        x_center,
        y_center
    )

    label = smart_label(art)

    fig.add_annotation(
        x=x_lab,
        y=y_lab,
        text=label,
        showarrow=False,  # No arrows
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


    from playwright.sync_api import sync_playwright

    # Print variance of major-ion variables
    print("Variances:")
    print(
        np.var(
            raw_df[ion_cols],
            axis=0
        )
    )

    # Print correlation matrix of major-ion variables
    print("\nCorrelation matrix:")
    print(
        np.corrcoef(
            raw_df[ion_cols].values.T
        )
    )


    # ============================================================
    # EXPORT HTML AND PNG OUTPUTS
    # ============================================================

    fig.write_html(
        plot_output,
        include_plotlyjs="cdn",
        full_html=True
    )

    print(
        f"\nHTML file saved to:\n"
        f"→ {plot_output}"
    )

    png_output = (
        OUTPUT_DIR
        / "Metanumber_Plot_Ca_HCO3_Bands.png"
    )

    html_path = (
        "file:///"
        + str(plot_output).replace("\\", "/")
    )

    with sync_playwright() as p:
        browser = p.chromium.launch(
            headless=True
        )

        page = browser.new_page(
            viewport={
                "width": 2260,
                "height": 1210
            },
            device_scale_factor=2
        )

        page.goto(
            html_path,
            wait_until="domcontentloaded",
            timeout=120000
        )

        page.wait_for_timeout(3000)

        page.screenshot(
            path=str(png_output),
            full_page=True
        )

        browser.close()

    print(
        f"PNG file saved to:\n"
        f"→ {png_output}"
    )

    fig.show()

    # Print calculated Ca reference-band limits
    print("\nCa reference-band limits:")
    for r in results_ca:
        print(
            f"Ca={r['Ca']}%  ->  "
            f"y_min={r['y_min']:.2f}  "
            f"y_max={r['y_max']:.2f}"
        )

    # Print calculated HCO3 reference-band limits
    print("\nHCO3 reference-band limits:")
    for r in results_hco3:
        print(
            f"HCO3={r['HCO3']}%  ->  "
            f"x_min={r['x_min']:.2f}  "
            f"x_max={r['x_max']:.2f}"
        )

except Exception as e:
    print("Error while generating the plot:", e)

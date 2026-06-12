"""
NMF Spatial Analysis Pipeline
==============================
primary vs. brain metastasis comparison in EAC spatial transcriptomics data.

Figures generated
-----------------
1  - Per-patient paired proportions (all neighborhoods)
2  - Survivor-stratified paired lines (per neighborhood)
3  - Mean ± SEM summary (LTS vs STS)
4  - Overlay paired lines
5  - Overlay paired lines + statistics
6  - Overlay mean ± SEM
7  - Neighborhood cellular composition heatmap
8  - Violin plots for significant neighborhoods
9  - Representative spatial maps (individual patient pairs + all-patient overview)
10 - Highlighted significant neighborhood maps
11 - NH5 cell-type enrichment heatmap (log2 observed / expected)
11B - NH5 enrichment restricted to Long-term survivor Brain Met

"""

# ── Imports ───────────────────────────────────────────────────────────────────

import os
import re
import warnings

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import stats

warnings.filterwarnings("ignore")


# ── Configuration ─────────────────────────────────────────────────────────────

DATA_PATH = os.path.expanduser(
    "~/cellcharter_scvi_embedding.h5ad"
)

MAPPING_FILE = os.path.expanduser(
    "~/mapping.csv"
)

OUTPUT_DIR = "./FINAL_SPATIAL_PUBPIPE"

PRIMARY_LABEL = "Primary EAC"
MET_LABEL = "EAC Brain Met"
GROUP_ORDER = [PRIMARY_LABEL, MET_LABEL]

NH_TARGET = "5"

# Figure style constants
LABEL_FS  = 7
TITLE_FS  = 7.5
TICK_FS   = 6.5
DOT_S     = 26
LINE_ALPHA = 0.55
LINE_LW   = 0.9

GROUP_COLORS = {"Long-term": "#2166ac", "Short-term": "#d6604d"}

SHORT_CELL_LABELS = {
    "CEACAM-high tumor epithelial cells":              "CEACAM+ Tumor",
    "Cycling Tumor Cells":                             "Cycling Tumor",
    "Mucin-producing tumor cells":                     "Mucin Tumor",
    "Inflamed primary tumor epithelial cells":         "Inflamed Tumor",
    "Complement immunosuppressive macrophages (TAMs)": "Complement TAMs",
    "Systemic inflammatory macrophage program (TAMs)": "Inflammatory TAMs",
    "CAFs (Cancer associated fibroblasts)":            "CAFs",
    "Pericyte-enriched endothelial cells":             "Pericyte ECs",
    "Cytotoxic T cells":                               "Cytotoxic T",
    "B lymphocytes":                                   "B cells",
    "Plasma Cells":                                    "Plasma",
}

PREFERRED_CELL_ORDER = [
    "CEACAM+ Tumor", "Cycling Tumor", "Mucin Tumor", "Inflamed Tumor",
    "Complement TAMs", "Inflammatory TAMs", "CAFs", "Pericyte ECs",
    "Cytotoxic T", "B cells", "Plasma",
]


# ── Utility helpers ───────────────────────────────────────────────────────────

def shorten_label(label: str, max_len: int = 45) -> str:
    return label if len(label) <= max_len else label[:max_len - 3] + "..."


def shorten_celltype(name: str) -> str:
    return SHORT_CELL_LABELS.get(str(name), str(name))


def fmt_p(p: float) -> str:
    if np.isnan(p):  return "P = n/a"
    if p < 0.001:    return "P < 0.001 ***"
    if p < 0.01:     return f"P = {p:.3f} **"
    if p < 0.05:     return f"P = {p:.3f} *"
    return f"P = {p:.2f} ns"


def wilcoxon_p(a, b) -> float:
    """Paired Wilcoxon signed-rank test; falls back to t-test for n < 4."""
    if len(a) < 4:
        return stats.ttest_rel(a, b).pvalue
    return stats.wilcoxon(a, b).pvalue


def style_ax(ax) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Primary", "Brain Met"], fontsize=TICK_FS)
    ax.set_xlim(-0.45, 1.45)
    ax.tick_params(axis="y", labelsize=TICK_FS)


def save_fig(fig, stem: str) -> None:
    fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")


def reorder_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Reorder heatmap columns by preferred biological order."""
    ordered = [c for c in PREFERRED_CELL_ORDER if c in df.columns]
    rest    = [c for c in df.columns if c not in ordered]
    return df[ordered + rest]


# ── Step 1: Load data ─────────────────────────────────────────────────────────

def load_adata(path: str):
    print("Loading AnnData object...")
    adata = sc.read_h5ad(path)
    print(f"AnnData shape: {adata.shape}")
    print(f"First 20 obs columns: {list(adata.obs.columns[:20])}")
    return adata


# ── Step 2: Identify sample column ───────────────────────────────────────────

def find_sample_column(adata) -> str:
    for col in adata.obs.columns:
        if any(k in col.lower() for k in ("sample", "patient", "pid")):
            print(f"Using sample column: {col}")
            return col
    raise ValueError("Could not find a sample column in adata.obs.")


# ── Step 3: Build per-sample metadata ────────────────────────────────────────

def _normalize_acc(s: str) -> str:
    if pd.isna(s):
        return ""
    return re.sub(r"[^\w\-]", "", str(s).lower())


def _match_sample(sample_name: str, mapping_df: pd.DataFrame):
    key = _normalize_acc(sample_name)
    m = mapping_df[mapping_df["Acc1_key"] == key]
    if not m.empty:
        return m.iloc[0]
    m = mapping_df[mapping_df["Acc1_key"].apply(lambda x: x in key or key in x)]
    if not m.empty:
        return m.iloc[0]
    print(f"⚠️  No mapping found for sample: {sample_name}")
    return None


def classify_survivor(pid: str) -> str:
    try:
        n = int(str(pid).split("-")[1])
    except (IndexError, ValueError):
        return "unclassified"
    if 1  <= n <= 14: return "Long-term"
    if 45 <= n <= 65: return "Short-term"
    return "unclassified"


def build_metadata(adata, sample_col: str) -> pd.DataFrame:
    mapping_df = pd.read_csv(MAPPING_FILE, dtype=str)
    mapping_df["Acc1_key"] = mapping_df["Acc1"].apply(_normalize_acc)

    rows = []
    for s in adata.obs[sample_col].unique():
        row = _match_sample(s, mapping_df)
        if row is None:
            rows.append({"sample": s, "base_pid": None, "Tumor Location": "unknown"})
        else:
            pid_match = re.search(r"(P-\d+)", row["Updated_ID-Dec2025"])
            rows.append({
                "sample": s,
                "base_pid": pid_match.group(1) if pid_match else None,
                "Tumor Location": row["Tumor Location"],
            })

    meta = pd.DataFrame(rows)
    meta["survivor_group"] = meta["base_pid"].apply(classify_survivor)
    return meta


# ── Step 4: Filter to paired patients ────────────────────────────────────────

def filter_paired_patients(adata, sample_col: str):
    sample_df = adata.obs[[sample_col, "base_pid", "Tumor Location"]].drop_duplicates()
    counts = (
        sample_df.groupby(["base_pid", "Tumor Location"])[sample_col]
        .nunique()
        .unstack(fill_value=0)
    )
    paired_pids = counts[
        (counts.get(PRIMARY_LABEL, 0) > 0) & (counts.get(MET_LABEL, 0) > 0)
    ].index
    print(f"Paired patients: {len(paired_pids)}")
    return adata[adata.obs["base_pid"].isin(paired_pids)].copy()


# ── Step 5: Assign semantic neighborhood labels ───────────────────────────────

def _map_to_broad_class(ct: str) -> str:
    ct = str(ct).lower()
    if "cancer associated fibroblasts" in ct or "caf" in ct:
        return "Fibroblast"
    if "ceacam" in ct or "tumor epithelial" in ct or "mucin" in ct or "cycling tumor" in ct:
        return "Tumor"
    if "cytotoxic t" in ct or "plasma" in ct:
        return "Adaptive-immune"
    if "macrophage" in ct:
        return "Myeloid"
    if "pericyte" in ct or "endothelial" in ct:
        return "Vascular"
    return "Other"


def _assign_tumor_subtype(tumor_props: pd.Series) -> str:
    if tumor_props.empty:
        return "Tumor-mixed"
    top = tumor_props.sort_values(ascending=False)
    name = top.index[0].lower()
    if "mucin" in name:
        return "Mucin-rich tumor" if top.iloc[0] > 50 else "Mucin-predominant tumor"
    if "cycling" in name:
        return "Cycling tumor–enriched"
    if "ceacam" in name:
        return "CEACAM-high tumor"
    return "Tumor-mixed"


def _assign_label(class_props: pd.Series, raw_props: pd.Series) -> str:
    if class_props.empty:
        return "Unclassified"
    top_class = class_props.idxmax()

    if top_class == "Tumor":
        tumor_only = raw_props[raw_props.index.map(lambda x: "tumor" in x.lower())]
        return _assign_tumor_subtype(tumor_only)
    if top_class == "Adaptive-immune":
        return "Immune-active" if class_props.get("Adaptive-immune", 0) > 40 else "Immune-mixed"
    if top_class == "Myeloid":
        return "Myeloid-rich" if class_props.get("Myeloid", 0) > 40 else "Myeloid-mixed"
    if top_class == "Fibroblast":
        return "Fibroblast-enriched"
    if top_class == "Vascular":
        return "Vascular-enriched"
    return "Mixed"


def build_neighborhood_composition(adata, neighborhoods) -> dict:
    composition = {}
    for nh in neighborhoods:
        cells      = adata.obs[adata.obs["spatial_cluster_k6"] == nh]
        raw_props  = cells["cell_type"].value_counts(normalize=True) * 100
        class_props = raw_props.groupby(raw_props.index.map(_map_to_broad_class)).sum()
        label       = f"NH{nh}: {_assign_label(class_props, raw_props)}"
        composition[nh] = {
            "label":       label,
            "class_props": class_props.to_dict(),
            "raw_props":   raw_props.nlargest(3).to_dict(),
        }
    return composition


# ── Step 6: Compute proportions ───────────────────────────────────────────────

def compute_proportions(adata, sample_col: str, neighborhoods) -> pd.DataFrame:
    rows = []
    for s in adata.obs[sample_col].unique():
        sdata = adata.obs[adata.obs[sample_col] == s]
        for tissue in GROUP_ORDER:
            tdata = sdata[sdata["Tumor Location"] == tissue]
            if len(tdata) == 0:
                continue
            total  = len(tdata)
            counts = tdata["spatial_cluster_k6"].value_counts()
            for nh in neighborhoods:
                rows.append({
                    "base_pid":       sdata["base_pid"].iloc[0],
                    "sample":         s,
                    "neighborhood":   nh,
                    "tissue":         tissue,
                    "proportion":     counts.get(nh, 0) / total,
                    "survivor_group": sdata["survivor_group"].iloc[0],
                })

    df = pd.DataFrame(rows)
    df_wide = df.pivot_table(
        index=["base_pid", "neighborhood", "survivor_group"],
        columns="tissue",
        values="proportion",
        fill_value=0,
    ).reset_index().rename(columns={
        PRIMARY_LABEL: "primary_proportion",
        MET_LABEL:     "brain_met_proportion",
    })
    return df, df_wide


# ── Step 7: Build spatial coordinate system ───────────────────────────────────

def detect_spatial_coords(adata) -> tuple:
    """Returns (use_obsm, x_col, y_col)."""
    if "spatial" in adata.obsm and adata.obsm["spatial"].shape[1] >= 2:
        print("✓ Using adata.obsm['spatial'] coordinates")
        return True, None, None

    preferred_pairs = [
        ("x_centroid", "y_centroid"), ("x_centroid_px", "y_centroid_px"),
        ("x", "y"), ("X", "Y"), ("x_coordinate", "y_coordinate"),
        ("xcoord", "ycoord"), ("spatial_x", "spatial_y"),
        ("center_x", "center_y"), ("Centroid X µm", "Centroid Y µm"),
        ("centroid_x", "centroid_y"),
    ]
    for xc, yc in preferred_pairs:
        if xc in adata.obs.columns and yc in adata.obs.columns:
            print(f"✓ Using coordinate columns: {xc}, {yc}")
            return False, xc, yc

    # Flexible fallback
    px, py = [], []
    for col in adata.obs.columns:
        c = col.lower()
        if ("x" in c and any(k in c for k in ("centroid", "coord", "spatial", "pixel", "center"))) or c == "x":
            px.append(col)
        if ("y" in c and any(k in c for k in ("centroid", "coord", "spatial", "pixel", "center"))) or c == "y":
            py.append(col)
    if px and py:
        print(f"✓ Auto-detected coordinates: {px[0]}, {py[0]}")
        return False, px[0], py[0]

    print("Available obs columns:", list(adata.obs.columns))
    print("Available obsm keys:",   list(adata.obsm.keys()))
    raise ValueError(
        "Could not detect spatial coordinates. "
        "Inspect printed columns/obsm keys and manually set x_col/y_col."
    )


def get_spatial_coords(subset, use_obsm: bool, x_col, y_col):
    if use_obsm:
        return subset.obsm["spatial"][:, 0], subset.obsm["spatial"][:, 1]
    return subset.obs[x_col].values, subset.obs[y_col].values


# ── Step 8: Build patient → sample lookup ────────────────────────────────────

def build_patient_sample_map(adata, sample_col: str) -> dict:
    sample_lookup = (
        adata.obs[[sample_col, "base_pid", "Tumor Location", "survivor_group"]]
        .drop_duplicates()
        .query("`Tumor Location` in @GROUP_ORDER")
    )

    patient_map = {}
    for pid in sorted(sample_lookup["base_pid"].dropna().unique()):
        sub        = sample_lookup[sample_lookup["base_pid"] == pid]
        tissue_map = {}
        for tissue in GROUP_ORDER:
            ts = sub[sub["Tumor Location"] == tissue]
            if len(ts) > 0:
                tissue_map[tissue] = ts[sample_col].iloc[0]
        if all(t in tissue_map for t in GROUP_ORDER):
            patient_map[pid] = {
                "survivor_group": sub["survivor_group"].iloc[0],
                "samples": tissue_map,
            }

    print(f"Fully paired mapped patients: {len(patient_map)}")
    return patient_map


# ── Spatial plot helpers ───────────────────────────────────────────────────────

def _make_nh_palette(adata, nh_colors: dict) -> tuple:
    adata.obs["spatial_cluster_k6"] = (
        adata.obs["spatial_cluster_k6"].astype(str).astype("category")
    )
    ordered_nh = sorted(adata.obs["spatial_cluster_k6"].cat.categories)
    nh_palette = {
        nh: nh_colors[int(nh)] if str(nh).isdigit() and int(nh) in nh_colors else "gray"
        for nh in ordered_nh
    }
    scanpy_palette = [nh_palette[nh] for nh in ordered_nh]
    return nh_palette, scanpy_palette


def plot_spatial_neighborhoods(ax, adata, sample_id, tissue_type, sample_col,
                                scanpy_nh_palette, title_prefix="") -> None:
    subset = adata[
        (adata.obs[sample_col] == sample_id) &
        (adata.obs["Tumor Location"] == tissue_type)
    ].copy()

    if subset.n_obs == 0:
        ax.axis("off")
        ax.set_title(f"No data: {sample_id} — {tissue_type}", fontsize=7)
        return

    subset.obs["spatial_cluster_k6"] = (
        subset.obs["spatial_cluster_k6"].astype(str).astype("category")
    )
    try:
        sc.pl.spatial(
            subset, color="spatial_cluster_k6", spot_size=30, frameon=False,
            palette=scanpy_nh_palette, legend_loc=None, ax=ax, show=False,
        )
    except Exception as e:
        ax.axis("off")
        ax.set_title(f"Plot failed: {sample_id}\n{str(e)[:50]}", fontsize=6)
        return

    ax.set_title(f"{title_prefix}{sample_id} — {tissue_type}",
                 fontsize=8, fontweight="bold")


def plot_highlighted_neighborhood(ax, adata, sample_id, tissue_type, target_nh,
                                   sample_col, nh_palette,
                                   neighborhood_composition, title_prefix="") -> None:
    subset = adata[
        (adata.obs[sample_col] == sample_id) &
        (adata.obs["Tumor Location"] == tissue_type)
    ].copy()

    if subset.n_obs == 0:
        ax.axis("off")
        ax.set_title(f"No data: {sample_id} — {tissue_type}", fontsize=7)
        return

    target_nh = str(target_nh)
    subset.obs["highlight_nh"] = pd.Categorical(
        np.where(subset.obs["spatial_cluster_k6"].astype(str) == target_nh,
                 f"NH{target_nh}", "Other"),
        categories=["Other", f"NH{target_nh}"], ordered=True,
    )
    highlight_palette = ["lightgray", nh_palette.get(target_nh, "red")]

    try:
        sc.pl.spatial(
            subset, color="highlight_nh", spot_size=30, frameon=False,
            palette=highlight_palette, legend_loc=None, ax=ax, show=False,
        )
    except Exception as e:
        ax.axis("off")
        ax.set_title(f"Plot failed: {sample_id}\n{str(e)[:50]}", fontsize=6)
        return

    nh_label = neighborhood_composition.get(
        int(target_nh), {"label": f"NH{target_nh}"}
    )["label"]
    ax.set_title(f"{title_prefix}{sample_id} — {tissue_type}\n{shorten_label(nh_label, 42)}",
                 fontsize=7, fontweight="bold")


def add_neighborhood_legend(fig, neighborhoods, neighborhood_composition,
                             nh_colors, ncol: int = 2) -> None:
    handles = [
        mpatches.Patch(
            facecolor=nh_colors.get(nh, "gray"), edgecolor="black", linewidth=0.3,
            label=shorten_label(neighborhood_composition[nh]["label"], 55),
        )
        for nh in neighborhoods
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.01),
               fontsize=6, ncol=ncol, frameon=False,
               title="Spatial Neighborhoods", title_fontsize=8)


# ── Schema validation & sample resolution (v3.2 pipeline) ────────────────────

def validate_schema(adata, df_wide) -> None:
    missing = {"base_pid", "survivor_group"} - set(df_wide.columns)
    if missing:
        raise ValueError(f"df_wide missing required columns: {missing}")
    if "sample" not in adata.obs.columns and "sample_id" not in adata.obs.columns:
        raise ValueError(
            "adata.obs must contain 'sample' or 'sample_id'. "
            "Map spatial samples before running this pipeline."
        )


def resolve_sample_column(adata) -> str:
    return "sample_id" if "sample_id" in adata.obs.columns else "sample"


def enforce_sample_mapping(df_wide: pd.DataFrame,
                            mapping_df: pd.DataFrame = None) -> pd.DataFrame:
    """Ensure df_wide has a 'sample_id' column.

    v3.3 fix: recognises the 'sample' column produced by the earlier pipeline
    stage and renames it so downstream schema-locked functions work correctly.
    """
    if "sample_id" in df_wide.columns:
        return df_wide
    if "sample" in df_wide.columns:
        df_wide = df_wide.copy()
        df_wide["sample_id"] = df_wide["sample"]
        print("✓ Mapped 'sample' column → 'sample_id' for schema compatibility")
        return df_wide
    if mapping_df is None:
        raise ValueError("No sample_id found. Provide mapping_df with base_pid → sample_id.")
    return df_wide.merge(mapping_df, on="base_pid", how="left")


def get_examples(df_wide: pd.DataFrame) -> tuple:
    if "sample_id" not in df_wide.columns:
        raise ValueError("sample_id missing — schema not resolved")
    lts = sorted(df_wide[df_wide["survivor_group"] == "Long-term"]["sample_id"].dropna().unique())
    sts = sorted(df_wide[df_wide["survivor_group"] == "Short-term"]["sample_id"].dropna().unique())
    return lts, sts


# ── NH5 enrichment helpers ────────────────────────────────────────────────────

def nh5_enrichment(adata, nh_target: str = NH_TARGET) -> pd.DataFrame:
    """Cohort-wide NH enrichment (observed/expected cell-type proportions)."""
    nh_cells  = adata.obs[adata.obs["spatial_cluster_k6"] == nh_target]
    overall   = adata.obs["cell_type"].value_counts(normalize=True)
    nh_ct     = nh_cells["cell_type"].value_counts(normalize=True)
    rows = [
        [ct, nh_ct.get(ct, 0), overall.get(ct, 0),
         nh_ct.get(ct, 0) / overall.get(ct, 1e-9)]
        for ct in overall.index
    ]
    return pd.DataFrame(rows, columns=["cell_type", "observed", "expected", "ratio"])


def _add_log2_ratio(df: pd.DataFrame, pseudocount: float = 1e-9) -> pd.DataFrame:
    df = df.copy()
    df["log2_ratio"] = np.log2(
        (df["observed"] + pseudocount) / (df["expected"] + pseudocount)
    )
    return df.sort_values("log2_ratio", ascending=False)


def save_nh5_with_log2(nh5_df: pd.DataFrame, output_dir: str,
                        nh_target: str = NH_TARGET) -> pd.DataFrame:
    nh5_df = _add_log2_ratio(nh5_df)
    out = os.path.join(output_dir, f"NH{nh_target}_enrichment_log2.csv")
    nh5_df.to_csv(out, index=False)
    print(f"✓ NH{nh_target} log2 enrichment CSV saved: {out}")
    return nh5_df


def nh5_lts_brainmet_enrichment(adata, sample_col: str,
                                  nh_target: str = NH_TARGET,
                                  survivor_group: str = "Long-term",
                                  tissue_label: str = MET_LABEL) -> pd.DataFrame:
    """NH enrichment restricted to Long-term survivor Brain Met tissue."""
    nh_target = str(nh_target)
    subset = adata[
        (adata.obs["survivor_group"] == survivor_group) &
        (adata.obs["Tumor Location"] == tissue_label)
    ].copy()
    if subset.n_obs == 0:
        raise ValueError(f"No cells for {survivor_group} + {tissue_label}")

    nh_subset = subset[subset.obs["spatial_cluster_k6"].astype(str) == nh_target].copy()
    if nh_subset.n_obs == 0:
        raise ValueError(f"No NH{nh_target} cells in {survivor_group} {tissue_label}")

    print(f"✓ Total LTS Brain Met cells: {subset.n_obs}")
    print(f"✓ NH{nh_target} LTS Brain Met cells: {nh_subset.n_obs}")

    overall = subset.obs["cell_type"].value_counts(normalize=True)
    nh_ct   = nh_subset.obs["cell_type"].value_counts(normalize=True)

    rows = [
        [ct, nh_ct.get(ct, 0), overall.get(ct, 0),
         nh_ct.get(ct, 0) / overall.get(ct, 1e-9) if overall.get(ct, 0) > 0 else np.nan,
         nh_subset.obs[sample_col].nunique()]
        for ct in overall.index
    ]
    df = pd.DataFrame(rows, columns=["cell_type", "observed", "expected", "ratio", "num_samples"])
    return _add_log2_ratio(df)


def nh5_lts_brainmet_enrichment_with_observed_counts(
    adata, output_dir: str, sample_col: str,
    nh_target: str = NH_TARGET,
    survivor_group: str = "Long-term",
    tissue_label: str = MET_LABEL,
):
    """Extended enrichment table that also saves raw observed/expected cell counts."""
    nh_target = str(nh_target)

    expected_subset = adata[adata.obs["Tumor Location"] == tissue_label].copy()
    expected_subset = expected_subset[
        expected_subset.obs["cell_type"] != "Inflamed primary tumor epithelial cells"
    ].copy()

    observed_subset = adata[
        (adata.obs["survivor_group"] == survivor_group) &
        (adata.obs["Tumor Location"] == tissue_label)
    ].copy()
    observed_subset = observed_subset[
        observed_subset.obs["cell_type"] != "Inflamed primary tumor epithelial cells"
    ].copy()

    if observed_subset.n_obs == 0:
        raise ValueError(f"No cells for {survivor_group} + {tissue_label}")

    nh_subset = observed_subset[
        observed_subset.obs["spatial_cluster_k6"].astype(str) == nh_target
    ].copy()
    if nh_subset.n_obs == 0:
        raise ValueError(f"No NH{nh_target} cells found")

    print(f"✓ Expected pool (all Brain Met): {expected_subset.n_obs}")
    print(f"✓ Observed pool (NH{nh_target} LTS BM): {nh_subset.n_obs}")

    overall = expected_subset.obs["cell_type"].value_counts(normalize=True)
    nh_ct   = nh_subset.obs["cell_type"].value_counts(normalize=True)

    expected_sample_counts = pd.crosstab(
        expected_subset.obs[sample_col], expected_subset.obs["cell_type"]
    )
    expected_total_counts = expected_subset.obs["cell_type"].value_counts()

    observed_sample_counts = pd.crosstab(
        nh_subset.obs[sample_col], nh_subset.obs["cell_type"]
    )
    observed_total_counts = nh_subset.obs["cell_type"].value_counts()

    rows = [
        [ct, nh_ct.get(ct, 0), overall.get(ct, 0),
         nh_ct.get(ct, 0) / overall.get(ct, 1e-9) if overall.get(ct, 0) > 0 else np.nan,
         observed_total_counts.get(ct, 0),
         expected_total_counts.get(ct, 0),
         nh_subset.obs[sample_col].nunique()]
        for ct in overall.index
    ]
    out_df = pd.DataFrame(rows, columns=[
        "cell_type", "observed", "expected", "ratio",
        "observed_total_count", "expected_total_count", "num_samples",
    ])
    out_df = _add_log2_ratio(out_df)

    # Save tables
    stems = {
        f"NH{nh_target}_LTS_BM_enrichment_with_observed_counts.csv": out_df,
        f"NH{nh_target}_EXPECTED_sample_by_celltype_counts.csv":       expected_sample_counts,
        f"NH{nh_target}_OBSERVED_sample_by_celltype_counts.csv":       observed_sample_counts,
        f"NH{nh_target}_EXPECTED_total_celltype_counts.csv":           expected_total_counts,
        f"NH{nh_target}_OBSERVED_total_celltype_counts.csv":           observed_total_counts,
    }
    for fname, data in stems.items():
        path = os.path.join(output_dir, fname)
        data.to_csv(path, index=(data is not out_df))
        print(f"  Saved: {path}")
    print("✓ Observed + expected count tables saved")

    return (out_df, expected_sample_counts, observed_sample_counts,
            expected_total_counts, observed_total_counts)


def print_lts_bm_nh_samples(adata, sample_col: str,
                              nh_target: str = NH_TARGET,
                              tissue_label: str = MET_LABEL) -> None:
    subset = adata[
        (adata.obs["survivor_group"] == "Long-term") &
        (adata.obs["Tumor Location"] == tissue_label) &
        (adata.obs["spatial_cluster_k6"].astype(str) == str(nh_target))
    ]
    print(f"\nLTS Brain Met samples containing NH{nh_target}:")
    print(subset.obs[sample_col].value_counts().sort_values(ascending=False))


# ── Figure generation ─────────────────────────────────────────────────────────

def fig1_per_patient(df_wide, neighborhoods, neighborhood_composition,
                     nh_colors) -> None:
    print("Generating Figure 1 (per-patient proportions)...")
    samples = sorted(df_wide["base_pid"].unique())
    ncols, nrows = 5, int(np.ceil(len(samples) / 5))

    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(ncols * 3.2, nrows * 3.6),
                              constrained_layout=True)
    axes = axes.flatten()

    for idx, pid in enumerate(samples):
        ax  = axes[idx]
        sub = df_wide[df_wide["base_pid"] == pid].set_index("neighborhood")
        group = sub["survivor_group"].iloc[0]
        gc    = GROUP_COLORS.get(group, "gray")

        pv_all, mv_all = [], []
        for nh in neighborhoods:
            if nh not in sub.index:
                continue
            pv, mv = sub.loc[nh, "primary_proportion"], sub.loc[nh, "brain_met_proportion"]
            pv_all.append(pv); mv_all.append(mv)
            lc = "#d62728" if pv > mv else nh_colors[nh]
            ax.plot([0, 1], [pv, mv], color=lc, alpha=LINE_ALPHA, linewidth=LINE_LW, zorder=1)

        ax.scatter(np.zeros(len(pv_all)), pv_all, color="white",
                   edgecolors="black", s=DOT_S, zorder=3, linewidths=0.7)
        ax.scatter(np.ones(len(mv_all)), mv_all, color="white",
                   edgecolors="black", s=DOT_S, zorder=3, linewidths=0.7)

        p = wilcoxon_p(pv_all, mv_all) if len(pv_all) >= 4 else np.nan
        tag = "[Lo]" if group == "Long-term" else "[Sh]" if group == "Short-term" else ""
        ax.set_title(f"{pid} {tag}\n{fmt_p(p)}", fontsize=TITLE_FS, color=gc, pad=3)
        ax.set_ylabel("Proportion", fontsize=LABEL_FS)
        style_ax(ax)

    for ax in axes[len(samples):]:
        ax.set_visible(False)

    legend_handles = [
        mpatches.Patch(facecolor=nh_colors[nh],
                       label=shorten_label(neighborhood_composition[nh]["label"]))
        for nh in neighborhoods
    ]
    fig.legend(handles=legend_handles, fontsize=5.5, loc="lower right",
               ncol=2, frameon=False, title="Neighborhood", title_fontsize=7)
    fig.suptitle(
        "Per-sample spatial neighborhood proportions — Primary vs Brain Met\n"
        "[Lo] = Long-term | [Sh] = Short-term | red = decrease in met",
        fontsize=10, fontweight="bold",
    )
    save_fig(fig, "fig1_pub")
    print("✓ Figure 1 saved.")


def fig2_stratified(df_wide, neighborhoods, neighborhood_composition) -> None:
    print("Generating Figure 2 (survivor-stratified)...")
    ncols = 3
    nrows = int(np.ceil(len(neighborhoods) / ncols))

    fig, axes = plt.subplots(nrows, ncols * 2,
                              figsize=(ncols * 2 * 2.8, nrows * 3.8),
                              constrained_layout=True)

    for idx, nh in enumerate(neighborhoods):
        row, col_base = idx // ncols, (idx % ncols) * 2
        for g_idx, group in enumerate(["Long-term", "Short-term"]):
            ax  = axes[row, col_base + g_idx]
            sub = df_wide[(df_wide["neighborhood"] == nh) &
                          (df_wide["survivor_group"] == group)]
            gc  = GROUP_COLORS[group]
            if sub.empty:
                ax.set_visible(False); continue
            pv  = sub["primary_proportion"].values
            mv  = sub["brain_met_proportion"].values
            for i in range(len(pv)):
                lc = "#d62728" if pv[i] > mv[i] else gc
                ax.plot([0, 1], [pv[i], mv[i]], color=lc, alpha=0.6, linewidth=1.0)
            ax.scatter(np.zeros(len(pv)), pv, color="white", edgecolors="black", s=DOT_S)
            ax.scatter(np.ones(len(mv)), mv,  color="white", edgecolors="black", s=DOT_S)
            p = wilcoxon_p(pv, mv) if len(pv) >= 4 else np.nan
            ax.set_title(
                f"{shorten_label(neighborhood_composition[nh]['label'])}\n"
                f"{group} (n={len(pv)})\n{fmt_p(p)}",
                fontsize=TITLE_FS, color=gc, fontweight="bold",
            )
            ax.set_facecolor((*plt.matplotlib.colors.to_rgb(gc), 0.04))
            ax.set_ylabel("Proportion", fontsize=LABEL_FS)
            style_ax(ax)

    fig.suptitle("Survivor-stratified paired lines — Primary vs Brain Met",
                 fontsize=10, fontweight="bold")
    save_fig(fig, "fig2_pub")
    print("✓ Figure 2 saved.")


def _mean_sem_grid(df_wide, neighborhoods, neighborhood_composition,
                   ncols: int = 3, markers: dict = None) -> plt.Figure:
    """Shared logic for Figures 3 and 6."""
    if markers is None:
        markers = {"Long-term": "o", "Short-term": "o"}
    nrows = int(np.ceil(len(neighborhoods) / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(ncols * 3.2, nrows * 3.6),
                              constrained_layout=True)
    axes = np.array(axes).reshape(-1)

    for idx, nh in enumerate(neighborhoods):
        ax       = axes[idx]
        y_offset = 0
        for group in ["Long-term", "Short-term"]:
            sub = df_wide[(df_wide["neighborhood"] == nh) &
                          (df_wide["survivor_group"] == group)]
            if len(sub) < 2:
                continue
            pv    = sub["primary_proportion"].values
            mv    = sub["brain_met_proportion"].values
            means = [np.mean(pv), np.mean(mv)]
            sems  = [stats.sem(pv), stats.sem(mv)]
            gc    = GROUP_COLORS[group]
            p     = wilcoxon_p(pv, mv)

            ax.fill_between([0, 1],
                            [means[0] - sems[0], means[1] - sems[1]],
                            [means[0] + sems[0], means[1] + sems[1]],
                            color=gc, alpha=0.15)
            ax.plot([0, 1], means, color=gc, linewidth=2.2,
                    label=f"{group} (n={len(pv)})")
            for xi, m, se in zip([0, 1], means, sems):
                ax.errorbar(xi, m, yerr=se, fmt=markers[group], color=gc,
                            markerfacecolor="white", markeredgewidth=1.5,
                            markersize=7, capsize=4)
            delta = means[1] - means[0]
            sign  = "+" if delta >= 0 else "−"
            ax.annotate(f"{group[:2]} Δ={sign}{abs(delta):.3f} {fmt_p(p)}",
                        xy=(0.5, max(means) + max(sems) + 0.02 + y_offset),
                        ha="center", fontsize=5.6, color=gc)
            y_offset += max(sems) + 0.03

        ax.set_title(shorten_label(neighborhood_composition[nh]["label"]), fontsize=TITLE_FS)
        ax.set_ylabel("Mean proportion", fontsize=LABEL_FS)
        ax.legend(fontsize=5.5, frameon=False)
        style_ax(ax)

    for ax in axes[len(neighborhoods):]:
        ax.set_visible(False)
    return fig


def fig3_mean_sem(df_wide, neighborhoods, neighborhood_composition) -> None:
    print("Generating Figure 3 (mean ± SEM)...")
    fig = _mean_sem_grid(df_wide, neighborhoods, neighborhood_composition)
    fig.suptitle("Cohort-level Mean ± SEM — LTS vs STS\nPrimary vs Brain Met",
                 fontsize=10, fontweight="bold")
    save_fig(fig, "fig3_pub")
    print("✓ Figure 3 saved.")


def _overlay_lines_grid(df_wide, neighborhoods, neighborhood_composition,
                         ncols: int = 3,
                         markers: dict = None,
                         show_stats: bool = False) -> plt.Figure:
    """Shared logic for Figures 4 and 5."""
    if markers is None:
        markers = {"Long-term": "o", "Short-term": "s"}
    nrows = int(np.ceil(len(neighborhoods) / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(ncols * 3.2, nrows * 3.6),
                              constrained_layout=True)
    axes = np.array(axes).reshape(-1)

    for idx, nh in enumerate(neighborhoods):
        ax = axes[idx]
        for group in ["Long-term", "Short-term"]:
            sub = df_wide[(df_wide["neighborhood"] == nh) &
                          (df_wide["survivor_group"] == group)]
            if sub.empty:
                continue
            pv = sub["primary_proportion"].values
            mv = sub["brain_met_proportion"].values
            for i in range(len(pv)):
                lc = "#d62728" if pv[i] > mv[i] else GROUP_COLORS[group]
                ax.plot([0, 1], [pv[i], mv[i]], color=lc,
                        alpha=LINE_ALPHA, linewidth=LINE_LW, zorder=1)
            ax.scatter(np.zeros(len(pv)), pv,
                       marker=markers[group], edgecolors=GROUP_COLORS[group],
                       facecolors="white", s=DOT_S, linewidths=0.9, zorder=3)
            ax.scatter(np.ones(len(mv)), mv,
                       marker=markers[group], edgecolors=GROUP_COLORS[group],
                       facecolors="white", s=DOT_S, linewidths=0.9, zorder=3)
            if show_stats:
                p = wilcoxon_p(pv, mv)
                y_pos = 1.02 if group == "Long-term" else 0.92
                ax.set_facecolor((*plt.matplotlib.colors.to_rgb(GROUP_COLORS[group]), 0.04))
                ax.text(0.5, y_pos, f"{group[:2]}: {fmt_p(p)}",
                        transform=ax.transAxes, ha="center",
                        fontsize=6, color=GROUP_COLORS[group])

        title_kw = {"fontweight": "bold"} if show_stats else {}
        ax.set_title(shorten_label(neighborhood_composition[nh]["label"]),
                     fontsize=TITLE_FS, **title_kw)
        ax.set_ylabel("Proportion", fontsize=LABEL_FS)
        style_ax(ax)

    for ax in axes[len(neighborhoods):]:
        ax.set_visible(False)
    return fig


def fig4_overlay(df_wide, neighborhoods, neighborhood_composition) -> None:
    print("Generating Figure 4 (overlay lines)...")
    fig = _overlay_lines_grid(df_wide, neighborhoods, neighborhood_composition)
    fig.suptitle("Overlay Paired Lines — LTS vs STS\n"
                 "Primary vs Brain Met (circles=LTS, squares=STS, red=decrease)",
                 fontsize=10, fontweight="bold")
    save_fig(fig, "fig4_overlay_pub")
    print("✓ Figure 4 saved.")


def fig5_overlay_stats(df_wide, neighborhoods, neighborhood_composition) -> None:
    print("Generating Figure 5 (overlay lines + stats)...")
    fig = _overlay_lines_grid(df_wide, neighborhoods, neighborhood_composition,
                               show_stats=True)
    fig.suptitle("Overlay Paired Lines with Statistics — LTS vs STS",
                 fontsize=10, fontweight="bold")
    save_fig(fig, "fig5_overlay_stats_pub")
    print("✓ Figure 5 saved.")


def fig6_overlay_mean_sem(df_wide, neighborhoods, neighborhood_composition) -> None:
    print("Generating Figure 6 (overlay mean ± SEM)...")
    fig = _mean_sem_grid(df_wide, neighborhoods, neighborhood_composition,
                          markers={"Long-term": "o", "Short-term": "s"})
    fig.suptitle("Overlay Mean ± SEM — LTS vs STS (Δ = Met − Primary)",
                 fontsize=10, fontweight="bold")
    save_fig(fig, "fig6_overlay_mean_sem_pub")
    print("✓ Figure 6 saved.")


def fig7_composition_heatmap(adata, neighborhoods, neighborhood_composition) -> None:
    print("Generating Figure 7 (cellular composition heatmap)...")
    records = []
    for nh in neighborhoods:
        cells = adata.obs[adata.obs["spatial_cluster_k6"] == nh]
        props = cells["cell_type"].value_counts(normalize=True) * 100
        props.name = nh
        records.append(props)
    comp_df = pd.DataFrame(records).fillna(0).loc[neighborhoods]

    fig = plt.figure(figsize=(max(12, len(comp_df.columns) * 0.45), 8))
    sns.heatmap(comp_df, cmap="viridis", annot=False, linewidths=0.3,
                cbar_kws={"label": "% of neighborhood"},
                yticklabels=[shorten_label(neighborhood_composition[nh]["label"], 55)
                             for nh in neighborhoods])
    plt.title("Cellular Composition of Spatial Neighborhoods",
              fontsize=12, fontweight="bold")
    plt.xlabel("Cell Type", fontsize=LABEL_FS)
    plt.ylabel("Neighborhood", fontsize=LABEL_FS)
    plt.xticks(rotation=90, fontsize=6)
    plt.yticks(fontsize=7)
    plt.tight_layout()
    save_fig(fig, "fig7_neighborhood_celltype_heatmap")
    print("✓ Figure 7 saved.")


def select_significant_neighborhoods(df_wide, neighborhoods,
                                      neighborhood_composition) -> list:
    results = []
    for nh in neighborhoods:
        sub = df_wide[df_wide["neighborhood"] == nh]
        if len(sub) < 3:
            continue
        p     = wilcoxon_p(sub["primary_proportion"].values,
                           sub["brain_met_proportion"].values)
        delta = sub["brain_met_proportion"].mean() - sub["primary_proportion"].mean()
        results.append({"neighborhood": nh,
                        "label":        neighborhood_composition[nh]["label"],
                        "p_value":      p,
                        "delta":        delta})
    sig_df = pd.DataFrame(results).sort_values("p_value")
    sig_nh = sig_df[sig_df["p_value"] < 0.05]["neighborhood"].tolist()
    if not sig_nh:
        sig_nh = sig_df.head(min(2, len(sig_df)))["neighborhood"].tolist()
    print(f"Key neighborhoods selected: {sig_nh}")
    return sig_nh


def fig8_violin(df, sig_neighborhoods, neighborhood_composition) -> None:
    print("Generating Figure 8 (violin plots)...")
    if not sig_neighborhoods:
        return

    violin_palette = {
        "Long-term":   GROUP_COLORS.get("Long-term",  "#2166ac"),
        "Short-term":  GROUP_COLORS.get("Short-term", "#d6604d"),
        "unclassified": "#9e9e9e",
    }

    fig, axes = plt.subplots(len(sig_neighborhoods), 1,
                              figsize=(8, 4.5 * len(sig_neighborhoods)),
                              constrained_layout=True)
    if len(sig_neighborhoods) == 1:
        axes = [axes]

    for ax, nh in zip(axes, sig_neighborhoods):
        sub = df[df["neighborhood"] == nh].copy()
        hue_order = [g for g in ["Long-term", "Short-term", "unclassified"]
                     if g in sub["survivor_group"].unique()]
        sns.violinplot(data=sub, x="tissue", y="proportion",
                       hue="survivor_group", hue_order=hue_order,
                       palette=violin_palette, split=False,
                       inner="box", cut=0, ax=ax)
        ax.set_title(shorten_label(neighborhood_composition[nh]["label"]),
                     fontsize=TITLE_FS, fontweight="bold")
        ax.set_ylabel("Neighborhood proportion", fontsize=LABEL_FS)
        ax.set_xlabel("")
        ax.tick_params(axis="x", labelsize=TICK_FS)
        ax.tick_params(axis="y", labelsize=TICK_FS)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, title="Survivor Group",
                  fontsize=6, title_fontsize=7, frameon=False)
        ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle("Distribution of Key Spatial Neighborhoods\n"
                 "Primary vs Brain Met by Survivor Group",
                 fontsize=11, fontweight="bold")
    save_fig(fig, "fig8_violin_key_neighborhoods")
    print("✓ Figure 8 saved.")


def fig9_individual_pairs(all_patients, patient_sample_map, adata, sample_col,
                           neighborhoods, neighborhood_composition,
                           nh_colors, scanpy_nh_palette) -> None:
    print("Generating Figure 9 individual patient-pair files...")
    for survivor_type, pid in all_patients:
        fig, axes = plt.subplots(1, 2, figsize=(14, 7), constrained_layout=True)
        for col_idx, tissue in enumerate(GROUP_ORDER):
            real_id = patient_sample_map[pid]["samples"].get(tissue)
            if real_id is None:
                axes[col_idx].axis("off"); continue
            tag = "Lo" if survivor_type == "Long-term" else "Sh"
            plot_spatial_neighborhoods(axes[col_idx], adata, real_id, tissue, sample_col,
                                        scanpy_nh_palette, title_prefix=f"[{tag}] {pid}\n")
        add_neighborhood_legend(fig, neighborhoods, neighborhood_composition, nh_colors, ncol=3)
        fig.suptitle(f"Spatial Neighborhood Architecture — {pid} ({survivor_type})",
                     fontsize=13, fontweight="bold")
        safe_pid = str(pid).replace("/", "_")
        save_fig(fig, f"fig9_pair_{safe_pid}_spatial")
        plt.close(fig)
    print("✓ Individual patient-pair files saved.")


def fig9_overview(all_patients, patient_sample_map, adata, sample_col,
                   scanpy_nh_palette) -> None:
    print("Generating Figure 9 all-patient overview...")
    nrows = len(all_patients)
    fig, axes = plt.subplots(nrows, 2, figsize=(16, 4 * nrows), constrained_layout=True)
    if nrows == 1:
        axes = np.array([axes])

    for row_idx, (survivor_type, pid) in enumerate(all_patients):
        for col_idx, tissue in enumerate(GROUP_ORDER):
            ax      = axes[row_idx, col_idx]
            real_id = patient_sample_map[pid]["samples"].get(tissue)
            if real_id is None:
                ax.axis("off"); continue
            tag = "Lo" if survivor_type == "Long-term" else "Sh"
            plot_spatial_neighborhoods(ax, adata, real_id, tissue, sample_col,
                                        scanpy_nh_palette, title_prefix=f"[{tag}] {pid}\n")

    for col_idx, tissue in enumerate(GROUP_ORDER):
        axes[0, col_idx].set_title(tissue, fontsize=11, fontweight="bold")

    fig.suptitle("Spatial Neighborhood Architecture — All Paired Patients",
                 fontsize=14, fontweight="bold")
    save_fig(fig, "fig9_all_patients_true_sampleIDs")
    print("✓ Figure 9 overview saved.")


def fig10_highlighted(all_patients, patient_sample_map, adata, sample_col,
                       sig_neighborhoods, nh_colors, nh_palette,
                       neighborhood_composition) -> None:
    print("Generating Figure 10 (highlighted neighborhoods)...")
    for nh in sig_neighborhoods:
        nrows = len(all_patients)
        fig, axes = plt.subplots(nrows, 2, figsize=(16, 4 * nrows), constrained_layout=True)
        if nrows == 1:
            axes = np.array([axes])

        for row_idx, (survivor_type, pid) in enumerate(all_patients):
            for col_idx, tissue in enumerate(GROUP_ORDER):
                ax      = axes[row_idx, col_idx]
                real_id = patient_sample_map[pid]["samples"].get(tissue)
                if real_id is None:
                    ax.axis("off"); continue
                tag = "Lo" if survivor_type == "Long-term" else "Sh"
                plot_highlighted_neighborhood(
                    ax, adata, real_id, tissue, nh, sample_col,
                    nh_palette, neighborhood_composition,
                    title_prefix=f"[{tag}] {pid}\n",
                )

        for col_idx, tissue in enumerate(GROUP_ORDER):
            axes[0, col_idx].set_title(tissue, fontsize=11, fontweight="bold")

        highlight_handles = [
            mpatches.Patch(color="lightgray", label="Other neighborhoods"),
            mpatches.Patch(color=nh_colors.get(nh, "red"),
                           label=shorten_label(neighborhood_composition[nh]["label"], 55)),
        ]
        fig.legend(handles=highlight_handles, loc="lower center",
                   bbox_to_anchor=(0.5, -0.01), ncol=2, frameon=False, fontsize=7)
        fig.suptitle(
            f"Spatial Distribution of {shorten_label(neighborhood_composition[nh]['label'])}\n"
            "All Paired Patients",
            fontsize=14, fontweight="bold",
        )
        save_fig(fig, f"fig10_true_sampleIDs_nh{nh}")
        plt.close(fig)
    print("✓ Figure 10 saved.")


def fig11_enrichment_heatmap(nh5_df: pd.DataFrame, output_dir: str,
                               nh_target: str = NH_TARGET) -> None:
    print(f"Generating Figure 11 (NH{nh_target} enrichment heatmap)...")
    required = {"cell_type", "observed", "expected", "ratio"}
    missing  = required - set(nh5_df.columns)
    if missing:
        raise ValueError(f"NH dataframe missing required columns: {missing}")

    nh5_plot = _add_log2_ratio(nh5_df)
    heatmap_df = nh5_plot.set_index("cell_type")[["log2_ratio"]]
    fig_height = max(8, len(heatmap_df) * 0.38)

    plt.figure(figsize=(7, fig_height))
    sns.heatmap(heatmap_df, cmap="bwr", center=0, annot=True, fmt=".2f",
                linewidths=0.3, linecolor="lightgray",
                cbar_kws={"label": "log2(Observed / Expected)"})
    plt.title(f"NH{nh_target} Cell-Type Enrichment\nlog2(Observed / Expected)",
              fontsize=13, fontweight="bold", pad=12)
    plt.ylabel("Cell Type", fontsize=10)
    plt.xticks(fontsize=9); plt.yticks(fontsize=8)
    plt.tight_layout()

    for ext in ("png", "pdf"):
        path = os.path.join(output_dir, f"FIG11_NH{nh_target}_log2_enrichment_heatmap.{ext}")
        plt.savefig(path, dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close()
    print(f"✓ Figure 11 saved to {output_dir}")


def fig11b_lts_bm_heatmap(nh_df: pd.DataFrame, output_dir: str,
                            nh_target: str = NH_TARGET) -> None:
    print(f"Generating Figure 11B (LTS Brain Met NH{nh_target} heatmap)...")
    heatmap_df = nh_df.set_index("cell_type")[["log2_ratio"]]
    fig_height = max(8, len(heatmap_df) * 0.38)

    plt.figure(figsize=(7.5, fig_height))
    sns.heatmap(heatmap_df, cmap="bwr", center=0, annot=True, fmt=".2f",
                linewidths=0.3, linecolor="lightgray",
                cbar_kws={"label": "log2(Observed / Expected)"})
    plt.title(f"LTS Brain Met Only — NH{nh_target} Cell-Type Enrichment\n"
              "log2(Observed / Expected)",
              fontsize=13, fontweight="bold", pad=12)
    plt.ylabel("Cell Type")
    plt.xticks(fontsize=9); plt.yticks(fontsize=8)
    plt.tight_layout()

    stems = {
        "png": os.path.join(output_dir, f"FIG11B_LTS_BrainMet_NH{nh_target}_log2_heatmap.png"),
        "pdf": os.path.join(output_dir, f"FIG11B_LTS_BrainMet_NH{nh_target}_log2_heatmap.pdf"),
    }
    plt.savefig(stems["png"], dpi=300, bbox_inches="tight")
    plt.savefig(stems["pdf"], bbox_inches="tight")
    plt.close()

    csv_path = os.path.join(output_dir, f"FIG11B_LTS_BrainMet_NH{nh_target}_log2_values.csv")
    nh_df.to_csv(csv_path, index=False)
    print(f"✓ Figure 11B saved to {output_dir}")


def _build_nh_heatmap_df(adata, nh_target: str) -> tuple[pd.DataFrame, dict]:
    """
    Per-sample NH enrichment table for LTS Brain Met.
    Returns (heatmap_df, count_map).
    """
    subset = adata[
        (adata.obs["survivor_group"] == "Long-term") &
        (adata.obs["Tumor Location"] == MET_LABEL)
    ].copy()
    if subset.n_obs == 0:
        return None, {}

    sample_col_local = resolve_sample_column(adata)
    global_expected  = subset.obs["cell_type"].value_counts(normalize=True)
    rows = []

    for sample_id in sorted(subset.obs[sample_col_local].dropna().unique()):
        s_data = subset[subset.obs[sample_col_local] == sample_id]
        nh_cells = s_data.obs[s_data.obs["spatial_cluster_k6"] == nh_target]
        if len(nh_cells) == 0:
            continue
        observed = nh_cells["cell_type"].value_counts(normalize=True)
        ratios   = {
            ct: np.log2((observed.get(ct, 0) + 1e-6) / (global_expected[ct] + 1e-6))
            for ct in global_expected.index
        }
        ratios["sample_id"] = sample_id
        rows.append(ratios)

    if not rows:
        return None, {}

    df = pd.DataFrame(rows).set_index("sample_id")
    df.columns = [shorten_celltype(c) for c in df.columns]
    df = reorder_columns(df)
    df["__sort__"] = df.abs().mean(axis=1)
    df = df.sort_values("__sort__", ascending=False).drop(columns="__sort__")

    count_map = (
        subset.obs[subset.obs["spatial_cluster_k6"] == nh_target]["cell_type"]
        .value_counts()
        .to_dict()
    )
    return df, count_map


def fig_nh_horizontal_heatmap(adata, output_dir: str,
                                nh_target: str = NH_TARGET) -> None:
    """Per-sample horizontal heatmap: NH enrichment in LTS Brain Met."""
    print(f"\n[NH{nh_target} LTS BM HORIZONTAL HEATMAP]")
    heatmap_df, _ = _build_nh_heatmap_df(adata, nh_target)
    if heatmap_df is None:
        print("⚠️  No data — skipped."); return

    plt.figure(figsize=(max(16, len(heatmap_df.columns) * 1.15),
                         max(5, len(heatmap_df.index) * 0.9)))
    ax = sns.heatmap(heatmap_df, cmap="RdBu_r", center=0,
                     linewidths=0.5, linecolor="lightgray",
                     cbar_kws={"label": "log2(observed / expected)"})

    plt.title(f"NH{nh_target} Cell-Type Enrichment — Long-Term Survivor Brain Metastases",
              fontsize=20, fontweight="bold", pad=18)
    plt.xlabel("Cell Type",           fontsize=15, fontweight="bold")
    plt.ylabel("LTS Brain Met Sample", fontsize=15, fontweight="bold")
    plt.xticks(rotation=45, ha="right", fontsize=13, fontweight="bold")
    plt.yticks(rotation=0, fontsize=13, fontweight="bold")
    ax.collections[0].colorbar.set_label("log2(observed / expected)",
                                          fontsize=14, fontweight="bold")
    plt.tight_layout()

    stem = os.path.join(output_dir, f"NH{nh_target}_LTS_BM_horizontal_heatmap_log2_custom")
    plt.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{stem}.pdf", bbox_inches="tight")
    heatmap_df.to_csv(f"{stem}.csv")
    plt.close()
    print(f"✓ Saved: {stem}.*")


def fig_nh_aggregated_heatmap(adata, output_dir: str,
                               nh_target: str = NH_TARGET,
                               cell_types_to_include: list = None) -> None:
    """Aggregated single-row horizontal heatmap for NH enrichment in LTS Brain Met."""
    print(f"\n[NH{nh_target} LTS BM AGGREGATED HORIZONTAL HEATMAP]")

    subset = adata[
        (adata.obs["survivor_group"] == "Long-term") &
        (adata.obs["Tumor Location"] == MET_LABEL)
    ].copy()
    if subset.n_obs == 0:
        print("⚠️  No data — skipped."); return

    expected = subset.obs["cell_type"].value_counts(normalize=True)
    nh_cells = subset.obs[subset.obs["spatial_cluster_k6"] == nh_target]
    if len(nh_cells) == 0:
        print("⚠️  No NH cells found — skipped."); return

    observed     = nh_cells["cell_type"].value_counts(normalize=True)
    enrich_dict  = {
        ct: np.log2((observed.get(ct, 0) + 1e-6) / (expected[ct] + 1e-6))
        for ct in expected.index
    }

    heatmap_df = pd.DataFrame([enrich_dict], index=[f"LTS Brain Met NH{nh_target}"])
    heatmap_df.columns = [shorten_celltype(c) for c in heatmap_df.columns]
    heatmap_df = reorder_columns(heatmap_df)

    if cell_types_to_include is not None:
        keep = [shorten_celltype(c) for c in cell_types_to_include
                if shorten_celltype(c) in heatmap_df.columns]
        heatmap_df = heatmap_df[keep]

    obs_counts = {shorten_celltype(ct): nh_cells["cell_type"].value_counts().get(ct, 0)
                  for ct in nh_cells["cell_type"].unique()}

    plt.figure(figsize=(max(16, len(heatmap_df.columns) * 1.2), 4.8))
    ax = sns.heatmap(heatmap_df, cmap="RdBu_r", center=0, linewidths=0.6,
                     linecolor="lightgray", annot=True, fmt=".2f",
                     annot_kws={"fontsize": 11, "fontweight": "bold"},
                     cbar_kws={"label": "log2(observed / expected)"})

    for i, ct in enumerate(heatmap_df.columns):
        ax.text(i + 0.5, -0.32, f"n={obs_counts.get(ct, 0)}",
                ha="center", va="center", fontsize=9, fontweight="bold", rotation=90)

    plt.title(f"Aggregated NH{nh_target} Cell-Type Enrichment — LTS Brain Metastases",
              fontsize=18, fontweight="bold", pad=40)
    plt.xlabel("Cell Type", fontsize=14, fontweight="bold")
    plt.xticks(rotation=45, ha="right", fontsize=12, fontweight="bold")
    plt.yticks(rotation=0, fontsize=12, fontweight="bold")
    ax.collections[0].colorbar.set_label("log2(observed / expected)",
                                          fontsize=13, fontweight="bold")
    plt.tight_layout()

    stem = os.path.join(output_dir, f"NH{nh_target}_LTS_BM_aggregated_heatmap_observed_counts")
    plt.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{stem}.pdf", bbox_inches="tight")
    heatmap_df.to_csv(f"{stem}.csv")
    plt.close()
    print(f"✓ Saved: {stem}.*")


# ── v3.2 schema-locked spatial example plot ───────────────────────────────────

def fig9_v32(adata, df_wide, output_dir: str) -> None:
    lts, sts = get_examples(df_wide)
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    sample_col_local = resolve_sample_column(adata)

    def _plot(ax, sample_id, title):
        sub = adata[adata.obs[sample_col_local] == sample_id].copy()
        if len(sub) == 0:
            ax.axis("off"); ax.set_title(f"{title}\n(NO DATA)"); return
        sc.pl.spatial(sub, color="spatial_cluster_k6", spot_size=25,
                      alpha=0.9, ax=ax, show=False, frameon=False)
        ax.set_title(title, fontsize=8, fontweight="bold")

    for i in range(3):
        if i < len(lts): _plot(axes[0, i], lts[i], f"LTS {lts[i]}")
        if i < len(sts): _plot(axes[1, i], sts[i], f"STS {sts[i]}")

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "FIG9_spatial_v32.png"), dpi=300)
    plt.close()
    print("✓ FIG9 v3.2 saved.")


# ── Main pipeline ─────────────────────────────────────────────────────────────

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ── Load & prepare ────────────────────────────────────────────────────────
    adata      = load_adata(DATA_PATH)
    sample_col = find_sample_column(adata)

    meta_df = build_metadata(adata, sample_col)
    adata.obs = adata.obs.merge(meta_df, left_on=sample_col,
                                right_on="sample", how="left")

    adata = filter_paired_patients(adata, sample_col)

    neighborhoods         = sorted(adata.obs["spatial_cluster_k6"].unique())
    neighborhood_composition = build_neighborhood_composition(adata, neighborhoods)

    df, df_wide = compute_proportions(adata, sample_col, neighborhoods)

    # ── Spatial coordinate system ─────────────────────────────────────────────
    use_obsm, x_col, y_col = detect_spatial_coords(adata)

    # ── Neighborhood color palette ────────────────────────────────────────────
    nh_colors = dict(zip(neighborhoods,
                         plt.cm.tab10(np.linspace(0, 1, len(neighborhoods)))))
    nh_palette, scanpy_nh_palette = _make_nh_palette(adata, nh_colors)

    # ── Patient sample map ────────────────────────────────────────────────────
    patient_sample_map = build_patient_sample_map(adata, sample_col)
    lts_patients = sorted(pid for pid, info in patient_sample_map.items()
                           if info["survivor_group"] == "Long-term")
    sts_patients = sorted(pid for pid, info in patient_sample_map.items()
                           if info["survivor_group"] == "Short-term")
    all_patients  = ([(  "Long-term", pid) for pid in lts_patients] +
                     [("Short-term", pid) for pid in sts_patients])

    print(f"LTS patients ({len(lts_patients)}): {lts_patients}")
    print(f"STS patients ({len(sts_patients)}): {sts_patients}")

    # ── Figures 1–8 ───────────────────────────────────────────────────────────
    fig1_per_patient(df_wide, neighborhoods, neighborhood_composition, nh_colors)
    fig2_stratified(df_wide, neighborhoods, neighborhood_composition)
    fig3_mean_sem(df_wide, neighborhoods, neighborhood_composition)
    fig4_overlay(df_wide, neighborhoods, neighborhood_composition)
    fig5_overlay_stats(df_wide, neighborhoods, neighborhood_composition)
    fig6_overlay_mean_sem(df_wide, neighborhoods, neighborhood_composition)
    fig7_composition_heatmap(adata, neighborhoods, neighborhood_composition)

    sig_neighborhoods = select_significant_neighborhoods(df_wide, neighborhoods,
                                                          neighborhood_composition)
    fig8_violin(df, sig_neighborhoods, neighborhood_composition)

    # ── Figures 9–10 (spatial) ────────────────────────────────────────────────
    fig9_individual_pairs(all_patients, patient_sample_map, adata, sample_col,
                           neighborhoods, neighborhood_composition,
                           nh_colors, scanpy_nh_palette)
    fig9_overview(all_patients, patient_sample_map, adata, sample_col,
                  scanpy_nh_palette)
    fig10_highlighted(all_patients, patient_sample_map, adata, sample_col,
                       sig_neighborhoods, nh_colors, nh_palette, neighborhood_composition)

    # ── v3.2 schema-locked components ────────────────────────────────────────
    mapping_rows = [
        {"base_pid": pid, "sample_id": sid}
        for pid, info in patient_sample_map.items()
        for sid in info["samples"].values()
    ]
    mapping_df = pd.DataFrame(mapping_rows).drop_duplicates()
    validate_schema(adata, df_wide)
    df_wide = enforce_sample_mapping(df_wide, mapping_df)
    fig9_v32(adata, df_wide, OUTPUT_DIR)

    # ── NH enrichment figures ─────────────────────────────────────────────────
    print("\n[NH5 ENRICHMENT — GLOBAL]")
    nh5_df = nh5_enrichment(adata)
    nh5_df = save_nh5_with_log2(nh5_df, OUTPUT_DIR)
    fig11_enrichment_heatmap(nh5_df, OUTPUT_DIR)

    print("\n[NH5 ENRICHMENT — LTS BRAIN MET]")
    nh5_lts_df = nh5_lts_brainmet_enrichment(adata, sample_col)
    fig11b_lts_bm_heatmap(nh5_lts_df, OUTPUT_DIR)
    print_lts_bm_nh_samples(adata, sample_col)

    print("\n[NH5 ENRICHMENT — OBSERVED COUNTS]")
    (nh5_lts_counts_df, *_) = nh5_lts_brainmet_enrichment_with_observed_counts(
        adata, OUTPUT_DIR, sample_col
    )

    fig_nh_horizontal_heatmap(adata, OUTPUT_DIR)
    fig_nh_aggregated_heatmap(adata, OUTPUT_DIR)

    selected_cell_types = [
        "Cytotoxic T cells",
        "Plasma Cells",
        "CEACAM-high tumor epithelial cells",
        "Cycling Tumor Cells",
    ]
    fig_nh_aggregated_heatmap(adata, OUTPUT_DIR,
                               cell_types_to_include=selected_cell_types)

    # ── Summary ───────────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("✓ ALL FIGURES GENERATED SUCCESSFULLY (v3.3)")
    print("=" * 70)
    figure_summary = [
        "1   Per-patient paired proportions",
        "2   Survivor-stratified paired lines",
        "3   Mean ± SEM",
        "4   Overlay paired lines",
        "5   Overlay paired lines + statistics",
        "6   Overlay mean ± SEM",
        "7   Neighborhood cellular composition heatmap",
        "8   Violin plots for significant neighborhoods",
        "9   Representative spatial maps (individual pairs + overview)",
        "10  Highlighted significant neighborhood maps",
        "11  NH5 enrichment heatmap (global)",
        "11B NH5 enrichment heatmap (LTS Brain Met only)",
    ]
    for line in figure_summary:
        print(f"  {line}")
    print(f"\nAll outputs saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()

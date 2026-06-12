"""
──────────────────────────────────────────────────────────────────────────────
Cohort-level paired plots of gene expression comparing Primary EAC vs EAC Brain Met,
with LTS and STS OVERLAID on the same axes for direct comparison.

For each gene, creates ONE panel with:
  • LTS lines in dark blue
  • STS lines in dark red
  • Separate dots for each survivor group
  • Two sets of statistics (one for LTS, one for STS)

Generates:
  Figure 1 – Overlay LTS/STS paired lines (1 panel per gene)
  Figure 2 – Overlay LTS/STS mean ± SEM (1 panel per gene)

Usage
──────
  import scanpy as sc
  adata = sc.read_h5ad("data.h5ad")
  main(adata)
"""

import os, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib import patches as mpatches

# ─────────────────────────────────────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────────────────────────────────────
genes_of_interest = ['CDK1', 'CDK4', 'TUBB', 'MKI67']



MALIGNANT_LABELS = [
    # --- Malignant epithelial ---
    'CEACAM-high tumor epithelial cells',   # softened azure blue
    'Cycling Tumor Cells',                 # softened amber
    'Mucin-producing tumor cells',         # softened coral
    'Inflamed primary tumor epithelial cells'
]


GROUP_ORDER = ["Primary EAC", "EAC Brain Met"]
PALETTE = {
    "Primary EAC":   "black",   # blue
    "EAC Brain Met": "black",   # orange
}

# Survivor-specific line colors (for overlay)
SURVIVOR_LINE_COLORS = {
    "Long-term":  "#2166AC",   # dark blue
    "Short-term": "#292929",   # dark red
}

# Survivor-specific marker styles
SURVIVOR_MARKERS = {
    "Long-term":  {"marker": "o", "facecolor": "white"},
    "Short-term": {"marker": "s", "facecolor": "white"},  # squares for STS
}

output_dir = "results"

mapping_file = "sample_metadata.csv"
IGNORE_PREFIX_N = 0

# Style constants
LABEL_FS   = 14
TITLE_FS   = 12
TICK_FS    = 12
DOT_S      = 26
LINE_ALPHA = 0.55
LINE_LW    = 0.9


# ═════════════════════════════════════════════════════════════════════════════
# 1. SURVIVOR CLASSIFICATION
# ═════════════════════════════════════════════════════════════════════════════
def classify_survivor(pid):
    """
    Assign samples to survival groups based on study-specific IDs.
    """
    try:
        n = int(str(pid).split("-")[1])
    except (IndexError, ValueError):
        return "unclassified"
    if 1 <= n <= 14:
        return "Long-term"
    if 45 <= n <= 65:
        return "Short-term"
    return "unclassified"


# ═════════════════════════════════════════════════════════════════════════════
# 2. SAMPLE → PATIENT MAPPING
# ═════════════════════════════════════════════════════════════════════════════
def normalize_acc1(s, ignore_n=IGNORE_PREFIX_N):
    if pd.isna(s):
        return ""
    s = str(s).strip().lower()
    s = re.sub(r'[^\w\-]', '', s)
    return s[ignore_n:]


def match_sample_acc1(sample_name, mapping_df):
    sample_key = normalize_acc1(sample_name)
    matches = mapping_df[mapping_df["Acc1_key"] == sample_key]
    if not matches.empty:
        return matches.iloc[0]
    matches = mapping_df[
        mapping_df["Acc1_key"].apply(
            lambda ak: ak in sample_key or sample_key in ak
        )
    ]
    if not matches.empty:
        return matches.iloc[0]
    print("Warning: sample not found in metadata mapping.")
    return None


def build_sample_metadata(adata, mapping_file, sample_col="sample"):
    mapping_df = pd.read_csv(mapping_file, dtype=str)
    mapping_df = mapping_df.loc[:, ~mapping_df.columns.duplicated()]
    mapping_df["Acc1_key"] = mapping_df["Acc1"].apply(normalize_acc1)
    mapping_df["mapped_ID_raw"] = mapping_df["Patient_ID"]

    meta_rows = []
    for sample in adata.obs[sample_col].unique():
        row = match_sample_acc1(sample, mapping_df)
        if row is None:
            meta_rows.append({
                "sample": sample, "mapped_ID_raw": sample,
                "base_pid": None, "Tumor Location": "unknown",
                "Sample ID": "unknown",
            })
        else:
            mapped_id = row["mapped_ID_raw"]
            pid_match = re.search(r"(P-\d+)", mapped_id)
            meta_rows.append({
                "sample": sample,
                "mapped_ID_raw": mapped_id,
                "base_pid": pid_match.group(1) if pid_match else None,
                "Tumor Location": row["Tumor Location"],
                "Sample ID": row["Sample ID"],
            })
    
    meta_df = pd.DataFrame(meta_rows)
    meta_df["survivor_group"] = meta_df["base_pid"].apply(classify_survivor)
    return meta_df


# ═════════════════════════════════════════════════════════════════════════════
# 3. DATA EXTRACTION
# ═════════════════════════════════════════════════════════════════════════════
def extract_malignant_expression(adata, genes_of_interest, sample_meta,
                                  malignant_labels=MALIGNANT_LABELS,
                                  sample_col="sample"):
    """Returns DataFrame with per-sample mean expression + survivor_group."""
    missing = [g for g in genes_of_interest if g not in adata.var_names]
    if missing:
        raise ValueError(f"Genes missing in adata: {missing}")

    mask = adata.obs['cell_type'].astype(str).isin(malignant_labels)
    adata_mal = adata[mask]
    print(f"Retained {mask.sum()} malignant cells ({mask.mean()*100:.1f}% of total)")

    gene_indices = {g: adata.var_names.get_loc(g) for g in genes_of_interest}
    sample_to_meta = sample_meta.set_index("sample")

    rows = []
    for sample in adata_mal.obs[sample_col].unique():
        if sample not in sample_to_meta.index:
            print(f"⚠️  Skipping unmapped sample: {sample}")
            continue
        meta = sample_to_meta.loc[sample]
        if meta["Tumor Location"] not in GROUP_ORDER:
            continue

        adata_s = adata_mal[adata_mal.obs[sample_col] == sample]
        row = {
            "sample": sample,
            "base_pid": meta["base_pid"],
            "Tumor Location": meta["Tumor Location"],
            "survivor_group": meta["survivor_group"],
        }
        for gene in genes_of_interest:
            idx = gene_indices[gene]
            expr = adata_s.layers['lognormal'][:, idx].toarray().flatten()
            row[f"{gene}_mean"] = float(np.mean(expr))
        rows.append(row)

    return pd.DataFrame(rows)


# ═════════════════════════════════════════════════════════════════════════════
# 4. FILTER TO PAIRED PATIENTS
# ═════════════════════════════════════════════════════════════════════════════
def filter_paired(sample_df, group_order=GROUP_ORDER):
    """Keep only samples from patients with ≥1 sample in each tumor location."""
    counts = (
        sample_df.groupby(["base_pid", "Tumor Location"])["sample"]
        .nunique()
        .unstack(fill_value=0)
    )
    paired_pids = counts[
        (counts.get(group_order[0], pd.Series(0, index=counts.index)) > 0) &
        (counts.get(group_order[1], pd.Series(0, index=counts.index)) > 0)
    ].index
    print(f"Paired patients retained: {len(paired_pids)}")
    return sample_df[sample_df["base_pid"].isin(paired_pids)].copy()


# ═════════════════════════════════════════════════════════════════════════════
# 5. PLOT HELPERS
# ═════════════════════════════════════════════════════════════════════════════
def style_ax(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(GROUP_ORDER, fontsize=TICK_FS)
    ax.set_xlim(-0.55, 1.65)
    ax.tick_params(axis="y", labelsize=TICK_FS)
    for tick, grp in zip(ax.get_xticklabels(), GROUP_ORDER):
        tick.set_color(PALETTE.get(grp, "black"))
        tick.set_fontweight("bold")


def wilcoxon_p_per_patient(sample_df, col, group_order=GROUP_ORDER):
    """Wilcoxon signed-rank using per-patient means."""
    patient_means = (
        sample_df.groupby(["base_pid", "Tumor Location"])[col]
        .mean()
        .unstack()
        .dropna(subset=group_order)
    )
    a = patient_means[group_order[0]].values
    b = patient_means[group_order[1]].values
    if len(a) < 4 or np.all(a == b):
        _, p = stats.ttest_rel(a, b)
    else:
        _, p = stats.wilcoxon(a, b, zero_method="wilcox", alternative="two-sided")
    return p


def fmt_p(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"


# ═════════════════════════════════════════════════════════════════════════════
# 6. FIGURE 1 – Overlay LTS/STS paired lines
# ═════════════════════════════════════════════════════════════════════════════
def plot_fig1_overlay(paired_df, genes_of_interest, plot_dir):
    """
    One panel per gene, with LTS and STS overlaid.
    Different line colors and marker shapes distinguish the groups.
    """
    survivor_groups = ["Long-term", "Short-term"]
    n_genes = len(genes_of_interest)
    ncols = min(2, n_genes)
    nrows = int(np.ceil(n_genes / ncols))
    
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(ncols * 3.8, nrows * 4.2),
        constrained_layout=True,
    )
    axes_flat = np.array(axes).flatten()
    
    pval_summary = {}
    
    for idx, gene in enumerate(genes_of_interest):
        ax = axes_flat[idx]
        col = f"{gene}_mean"
        
        # Plot each survivor group with different styling
        for survivor in survivor_groups:
            df_surv = paired_df[paired_df["survivor_group"] == survivor]
            
            if df_surv.empty:
                continue
            
            primary_df = df_surv[df_surv["Tumor Location"] == GROUP_ORDER[0]]
            met_df = df_surv[df_surv["Tumor Location"] == GROUP_ORDER[1]]
            
            line_color = SURVIVOR_LINE_COLORS[survivor]
            marker_style = SURVIVOR_MARKERS[survivor]
            
            # Paired lines with survivor-specific color
            for pid in df_surv["base_pid"].unique():
                p_vals = primary_df.loc[primary_df["base_pid"] == pid, col].values
                m_vals = met_df.loc[met_df["base_pid"] == pid, col].values
                for pv in p_vals:
                    for mv in m_vals:
                        # Use survivor color, but make decreasing lines slightly darker
                        if pv > mv:
                            plot_color = "#8B0000" if survivor == "Short-term" else "#000080"
                        else:
                            plot_color = line_color
                        
                        ax.plot(
                            [0, 1], [pv, mv],
                            color=plot_color, alpha=LINE_ALPHA,
                            linewidth=LINE_LW, zorder=1,
                        )
            
            # Sample dots with survivor-specific markers
        # Primary samples (x=0)
            ax.scatter(
                np.zeros(len(primary_df)), primary_df[col].values,
                color=marker_style["facecolor"],
                edgecolors=line_color,
                marker=marker_style["marker"],
                s=DOT_S, zorder=3, linewidths=1.2,
            )
            # Met samples (x=1)
            ax.scatter(
                np.ones(len(met_df)), met_df[col].values,
                color=marker_style["facecolor"],
                edgecolors=line_color,
                marker=marker_style["marker"],
                s=DOT_S, zorder=3, linewidths=1.2,
            )
            
            # Calculate p-value for this survivor group
            if len(df_surv["base_pid"].unique()) >= 2:
                p = wilcoxon_p_per_patient(df_surv, col)
                pval_summary[f"{gene}_{survivor}"] = p
            else:
                pval_summary[f"{gene}_{survivor}"] = np.nan
        
        # Add legend for survivor groups
        legend_elements = [
            mpatches.Patch(color=SURVIVOR_LINE_COLORS["Long-term"],
                          label=f'LTS ({fmt_p(pval_summary.get(f"{gene}_Long-term", np.nan))})'),
            mpatches.Patch(color=SURVIVOR_LINE_COLORS["Short-term"],
                          label=f'STS ({fmt_p(pval_summary.get(f"{gene}_Short-term", np.nan))})')
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=12,
                 frameon=True, fancybox=False, edgecolor='gray')
        
        # Count total samples by location (across both survivor groups)
        all_primary = paired_df[paired_df["Tumor Location"] == GROUP_ORDER[0]]
        all_met = paired_df[paired_df["Tumor Location"] == GROUP_ORDER[1]]
        
        for x_pos, grp_df in zip([0, 1], [all_primary, all_met]):
            # Count by survivor group
            lts_n = grp_df[grp_df["survivor_group"] == "Long-term"]["sample"].nunique()
            sts_n = grp_df[grp_df["survivor_group"] == "Short-term"]["sample"].nunique()
            total_patients = grp_df["base_pid"].nunique()
            
        
        ax.set_title(f"{gene}", fontsize=TITLE_FS, pad=3, fontweight='bold')
        ax.set_ylabel("Mean log-norm expression", fontsize=LABEL_FS)
        style_ax(ax)
        ax.set_ylim(0, 3.5)
    
    for ax in axes_flat[n_genes:]:
        ax.set_visible(False)
    
    fig.suptitle(
        "Gene Expression: Primary vs Brain Met – LTS & STS Overlaid\n"
        "(Malignant cells; circles=LTS, squares=STS; darker lines=decrease)",
        fontsize=10.5, fontweight="bold"
    )
    
    stem = "fig1_overlay_LTS_STS_paired_lines"
    fig.savefig(os.path.join(plot_dir, f"{stem}.png"), bbox_inches="tight", dpi=200)
    fig.savefig(os.path.join(plot_dir, f"{stem}.pdf"), bbox_inches="tight")
    plt.close(fig)
    print("Figure 1 (overlay) saved.")
    return pval_summary


# ═════════════════════════════════════════════════════════════════════════════
# 7. FIGURE 2 – Overlay LTS/STS mean ± SEM
# ═════════════════════════════════════════════════════════════════════════════
def plot_fig2_overlay(paired_df, genes_of_interest, pval_summary, plot_dir):
    """Mean ± SEM with LTS and STS overlaid on same axes."""
    survivor_groups = ["Long-term", "Short-term"]
    n_genes = len(genes_of_interest)
    ncols = min(2, n_genes)
    nrows = int(np.ceil(n_genes / ncols))
    
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(ncols * 3.6, nrows * 3.2),
        constrained_layout=True,
    )
    axes_flat = np.array(axes).flatten()
    
    for idx, gene in enumerate(genes_of_interest):
        ax = axes_flat[idx]
        col = f"{gene}_mean"
        
        # Plot each survivor group
        for survivor in survivor_groups:
            df_surv = paired_df[paired_df["survivor_group"] == survivor]
            
            if df_surv.empty:
                continue
            
            primary_vals = df_surv.loc[
                df_surv["Tumor Location"] == GROUP_ORDER[0], col
            ].values
            met_vals = df_surv.loc[
                df_surv["Tumor Location"] == GROUP_ORDER[1], col
            ].values
                
            
            means = [np.mean(primary_vals), np.mean(met_vals)]
            sems = [stats.sem(primary_vals), stats.sem(met_vals)]
        
            color = SURVIVOR_LINE_COLORS[survivor]
            marker = SURVIVOR_MARKERS[survivor]["marker"]
                
                # SEM band
            ax.fill_between(
                [0, 1],
                [means[0] - sems[0], means[1] - sems[1]],
                [means[0] + sems[0], means[1] + sems[1]],
                color=color, alpha=0.15, zorder=1,
            )
            
            # Mean line
            ax.plot([0, 1], means, color=color, linewidth=2.2,
                   label=f"{survivor}", zorder=2)
            
            # Error bars
            for xi, m, se in zip([0, 1], means, sems):
                ax.errorbar(
                    xi, m, yerr=se,
                    fmt=marker, color=color,
                    markerfacecolor="white", markeredgewidth=1.5,
                    markersize=8, capsize=5, capthick=1.5,
                    linewidth=1.5, zorder=3,
                )
            

        legend_elements = [
            mpatches.Patch(color=SURVIVOR_LINE_COLORS["Long-term"],
                          label=f'LTS {fmt_p(pval_summary.get(f"{gene}_Long-term", np.nan))}'),
            mpatches.Patch(color=SURVIVOR_LINE_COLORS["Short-term"],
                          label=f'STS {fmt_p(pval_summary.get(f"{gene}_Short-term", np.nan))}')
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=10,
                 frameon=False, fancybox=False)
        
        # Count annotations
        all_primary = paired_df[paired_df["Tumor Location"] == GROUP_ORDER[0]]
        all_met = paired_df[paired_df["Tumor Location"] == GROUP_ORDER[1]]
        
        for x_pos, grp_df in zip([0, 1], [all_primary, all_met]):
            lts_n = grp_df[grp_df["survivor_group"] == "Long-term"]["sample"].nunique()
            sts_n = grp_df[grp_df["survivor_group"] == "Short-term"]["sample"].nunique()
            total_patients = grp_df["base_pid"].nunique()
            
        
        ax.set_title(f"{gene}", fontsize=TITLE_FS, pad=3, fontweight='bold')
        ax.set_ylabel("Mean log-norm expression", fontsize=LABEL_FS)
        style_ax(ax)
        ax.set_ylim(0, 3.5)
    
    for ax in axes_flat[n_genes:]:
        ax.set_visible(False)
    
    fig.suptitle(
        "Mean ± SEM: Primary vs Brain Met – LTS & STS Overlaid\n"
        "(Malignant cells; SEM bands and error bars by survivor group)",
        fontsize=10.5, fontweight="bold"
    )
    
    stem = "fig2_overlay_LTS_STS_mean_sem"
    fig.savefig(os.path.join(plot_dir, f"{stem}.png"), bbox_inches="tight", dpi=200)
    fig.savefig(os.path.join(plot_dir, f"{stem}.pdf"), bbox_inches="tight")
    plt.close(fig)
    print("Figure 2 (overlay) saved.")


# ═════════════════════════════════════════════════════════════════════════════
# 8. P-VALUE SUMMARY
# ═════════════════════════════════════════════════════════════════════════════
def print_pval_summary(pval_summary):
    print("\n── P-values by Gene and Survivor Group ────────────────────────────")
    print(f"{'Gene_SurvivorGroup':<25} {'P-value':>9}  Sig")
    print("─" * 45)
    for key, p in pval_summary.items():
        if np.isnan(p):
            print(f"{key:<25} {'N/A':>9}  N/A")
        else:
            sig = ("***" if p < 0.001 else
                   "**"  if p < 0.01  else
                   "*"   if p < 0.05  else "ns")
            print(f"{key:<25} {p:>9.4f}  {sig}")


# ═════════════════════════════════════════════════════════════════════════════
# 9. MAIN
# ═════════════════════════════════════════════════════════════════════════════
def main(adata):
    plot_dir = os.path.join(output_dir, "cohort_gene_expression_by_survivor_overlay")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)

    print("Building sample metadata (with survivor groups)…")
    sample_meta = build_sample_metadata(adata, mapping_file)

    print("Extracting malignant cell expression (per sample)…")
    sample_df = extract_malignant_expression(adata, genes_of_interest, sample_meta)

    print("\nSurvivor group distribution:")
    print(sample_df.groupby(["Tumor Location", "survivor_group"])["sample"].nunique())

    print("\nFiltering to paired patients (within each survivor group)…")
    # Filter separately for each survivor group
    lts_df = sample_df[sample_df["survivor_group"] == "Long-term"]
    sts_df = sample_df[sample_df["survivor_group"] == "Short-term"]
    
    print("\n--- Long-term survivors ---")
    lts_paired = filter_paired(lts_df)
    
    print("\n--- Short-term survivors ---")
    sts_paired = filter_paired(sts_df)
    
    # Combine back
    paired_df = pd.concat([lts_paired, sts_paired], ignore_index=True)
    
    print(f"\nTotal paired samples: {len(paired_df)}")
    print(paired_df.groupby(["Tumor Location", "survivor_group"])["sample"].nunique())

    print("\nPlotting Figure 1 (overlay paired lines)…")
    pval_summary = plot_fig1_overlay(paired_df, genes_of_interest, plot_dir)

    print("Plotting Figure 2 (overlay mean ± SEM)…")
    plot_fig2_overlay(paired_df, genes_of_interest, pval_summary, plot_dir)

    print_pval_summary(pval_summary)
    print("\nDone. Outputs in:", plot_dir)

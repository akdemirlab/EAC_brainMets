#!/usr/bin/env python3
"""
spatial_niche_plot.py
=====================
Visualise a spatial neighbourhood ("niche") around a focal cell, replicating
the panel style from the reference figure:

  • Cell boundaries coloured by cell type (only specified types are shown)
  • Unspecified cell types are rendered as faint grey outlines (or hidden)
  • If a gene is supplied, transcripts for that gene are overlaid as small dots
    ONLY within the footprint of the shown cell types
  • A circular region-of-interest is drawn around the focal cell

Usage
-----
  # Centre on a specific cell ID
  python spatial_niche_plot.py \\
      --cell_boundaries  cell_boundaries.csv.gz \\
      --annotations      annotations.csv \\
      --transcripts      transcripts.zarr.zip \\
      --focal_cell_id    12345 \\
      --radius_um        120 \\
      --cell_types       "Malignant" "myCAF" "TAM" \\
      --gene             ERBB2 \\
      --outfile          niche_plot.pdf

  # Centre on x/y coordinates (nearest cell centroid is used as focus)
  python spatial_niche_plot.py \\
      --cell_boundaries  cell_boundaries.csv.gz \\
      --annotations      annotations.csv \\
      --x_um  1234.5 --y_um  678.9 \\
      --radius_um        120 \\
      --cell_types       "Malignant" "myCAF" "TAM" \\
      --outfile          niche_plot.pdf

Arguments
---------
--cell_boundaries   CSV/CSV.GZ with columns: cell_id, vertex_x, vertex_y
--annotations       CSV/TSV with columns: cell_id, cell_type  (+ optional extra cols)
--transcripts       transcripts.zarr or transcripts.zarr.zip
--focal_cell_id     The cell ID around which the niche is centred.
                    Mutually exclusive with --x_um / --y_um.
--x_um              X coordinate (in µm) of the niche centre.
                    The nearest cell centroid is found and used as the focal cell.
                    Must be supplied together with --y_um.
--y_um              Y coordinate (in µm) of the niche centre (see --x_um).
--radius_um         Radius of the circular ROI in microns (default: 120)
--cell_types        One or more cell_type labels to colour and show.
                    Cells NOT in this list are drawn as faint grey outlines
                    (use --hide_other to suppress them entirely).
--hide_other        Flag: if set, cells outside --cell_types are not drawn at all
--gene              Gene name to overlay as transcript dots (optional)
--min_qv            Minimum Phred QV for transcripts (default: 20)
--palette           Comma-separated HEX colours aligned with --cell_types
                    (auto-generated if omitted)
--figsize           Width,Height in inches (default: 5,5)
--dpi               Output DPI (default: 300)
--outfile           Output file path (default: niche_plot.pdf)
--px_scale          µm per coordinate unit (default: auto-detect)
--transcript_dot_size  Dot size for transcript scatter (default: 2)
--cell_alpha        Alpha for filled cell polygons (default: 0.75)
--other_alpha       Alpha for grey 'other' cell outlines (default: 0.25)
--show_circle       Draw the ROI circle boundary (default: True)
--title             Custom plot title (optional)

Dependencies
------------
  numpy, pandas, matplotlib, zarr, shapely
"""

import argparse
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle
from matplotlib.collections import PatchCollection, PolyCollection
import matplotlib.gridspec as gridspec

# ── optional shapely ───────────────────────────────────────────────────────────
try:
    from shapely.geometry import Point, Polygon, MultiPolygon
    from shapely.ops import unary_union
    HAS_SHAPELY = True
except ImportError:
    HAS_SHAPELY = False
    warnings.warn(
        "shapely not installed — transcript masking to cell footprints disabled. "
        "Run: pip install shapely",
        stacklevel=2,
    )

# ── optional zarr ──────────────────────────────────────────────────────────────
try:
    import zarr
    HAS_ZARR = True
except ImportError:
    HAS_ZARR = False
    warnings.warn("zarr not installed — transcripts cannot be loaded.", stacklevel=2)


# ══════════════════════════════════════════════════════════════════════════════
# CONSTANTS
# ══════════════════════════════════════════════════════════════════════════════

_NON_GENE_PREFIXES = (
    "negcontrolprobe_", "negcontrolcodeword_", "blank_",
    "antisense_", "unassigned", "deprecated", "intergenic_", "falsecodev",
)

# Colourblind-friendly default palette (tab10 extended)
_DEFAULT_COLOURS = [
    "#4E3B8B",  # dark purple
    "#C85FAB",  # pink/magenta
    "#2CA02C",  # green
    "#98DF8A",  # light green
    "#FF7F0E",  # orange
    "#D62728",  # red
    "#1F77B4",  # blue
    "#9467BD",  # purple
    "#8C564B",  # brown
    "#E377C2",  # pink
    "#7F7F7F",  # grey
    "#BCBD22",  # yellow-green
    "#17BECF",  # teal
]

_OTHER_CELL_COLOUR = "#CCCCCC"   # light grey for unspecified cells
_OTHER_CELL_LINEWIDTH = 0.4
_CIRCLE_COLOUR = "black"


_DEFAULT_GENE_COLOURS = [
    "black",  # deep sky blue
]

# ── Cell-type colour themes ────────────────────────────────────────────────────

# Theme 1: default (tab10-inspired, high contrast)
_THEME1_COLOURS = [
    "#4E3B8B", "#C85FAB", "#2CA02C", "#98DF8A",
    "#FF7F0E", "#D62728", "#1F77B4", "#9467BD",
    "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
]





# Theme 2: biology-informed palette
_THEME2_PALETTE: Dict[str, str] = {
    # Malignant epithelial
    "CEACAM-high tumor epithelial cells":       "#6BA4F8",
    "Cycling Tumor Cells":                      "#F4BA63",
    "Mucin-producing tumor cells":              "#E59973",
    "Inflamed primary tumor epithelial cells":  "#FF6259",
    # Tumour microenvironment
    "Complement immunosuppressive macrophages (TAMs)": "lightgray",
    "Systemic inflammatory macrophage program (TAMs)": "lightgray",
    "CAFs (Cancer associated fibroblasts)":     "lightgray",
    "Pericyte-enriched endothelial cells":      "lightgray",
    # Immune
    "Cytotoxic T cells":                        "#874284",
    "Plasma Cells":                             "#56B356",
}



# ══════════════════════════════════════════════════════════════════════════════
# I/O helpers
# ══════════════════════════════════════════════════════════════════════════════

def load_cell_boundaries(path: Path) -> pd.DataFrame:
    """
    Load cell boundary vertices.

    Accepts CSV or CSV.GZ.  Required columns: cell_id, vertex_x, vertex_y.
    Xenium also names them x_global_px / y_global_px — we handle both.
    """
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip().str.lower()

    # Normalise column names
    rename = {}
    for c in df.columns:
        if c in ("cell_id", "cellid", "cell id"):
            rename[c] = "cell_id"
        elif c in ("vertex_x", "x_global_px", "x", "x_um"):
            rename[c] = "vertex_x"
        elif c in ("vertex_y", "y_global_px", "y", "y_um"):
            rename[c] = "vertex_y"
    df = df.rename(columns=rename)

    required = {"cell_id", "vertex_x", "vertex_y"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"cell_boundaries file missing columns: {missing}. "
            f"Found: {list(df.columns)}"
        )
    df["cell_id"] = df["cell_id"].astype(str)
    return df[["cell_id", "vertex_x", "vertex_y"]]


def load_annotations(path: Path) -> pd.DataFrame:
    """Load cell-type annotation file (CSV or TSV, two-column minimum)."""
    sep = "\t" if str(path).endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep)
    df.columns = df.columns.str.strip().str.lower()

    rename = {}
    for c in df.columns:
        if c in ("cell_id", "cellid", "id"):
            rename[c] = "cell_id"
        if c in ("cell_type", "celltype", "type", "annotation", "label", "cluster"):
            rename[c] = "cell_type"
    df = df.rename(columns=rename)

    if "cell_id" not in df.columns or "cell_type" not in df.columns:
        raise ValueError(
            f"Annotation file must have 'cell_id' and 'cell_type' columns. "
            f"Found: {list(df.columns)}"
        )
    df["cell_id"] = df["cell_id"].astype(str)
    return df[["cell_id", "cell_type"]]


def detect_px_scale(boundaries_df: pd.DataFrame) -> float:
    """
    Heuristic: if coordinates exceed 10 000 assume pixel space (0.2125 µm/px),
    otherwise assume already in µm (scale = 1.0).
    """
    max_coord = max(boundaries_df["vertex_x"].max(), boundaries_df["vertex_y"].max())
    return 0.2125 if max_coord > 10_000 else 1.0


def build_cell_polygons(
    boundaries_df: pd.DataFrame,
    px_scale: float,
) -> Dict[str, np.ndarray]:
    """
    Group boundary rows by cell_id and return a dict:
      { cell_id: np.ndarray of shape (N, 2) in µm }
    """
    polys = {}
    for cid, grp in boundaries_df.groupby("cell_id"):
        pts = grp[["vertex_x", "vertex_y"]].values * px_scale
        polys[str(cid)] = pts
    return polys


# ══════════════════════════════════════════════════════════════════════════════
# Transcript loading — supports CSV(.gz), Parquet, and zarr transcript tables
# ══════════════════════════════════════════════════════════════════════════════

# Xenium transcript CSV/parquet column name candidates
_TX_X_COLS  = ["x_location", "x", "x_global_px", "x_px", "x_um"]
_TX_Y_COLS  = ["y_location", "y", "y_global_px", "y_px", "y_um"]
_TX_GENE_COLS = ["feature_name", "gene_name", "gene", "name"]
_TX_QV_COLS   = ["qv", "quality_value", "quality", "qv_score"]


def _load_transcript_table(path: Path) -> pd.DataFrame:
    """
    Load a transcript table from CSV(.gz), Parquet, or a zarr that contains
    a flat transcript table (keys: x_location, y_location, feature_name, qv).

    Returns a DataFrame with at minimum columns: x, y, gene  (and qv if present).
    """
    p = str(path)
    suffix = path.suffix.lower()
    stem_suffix = path.stem.lower()  # e.g. ".csv" from "transcripts.csv.gz"

    # ── parquet ────────────────────────────────────────────────────────────
    if suffix == ".parquet":
        try:
            df = pd.read_parquet(path)
        except Exception as e:
            raise RuntimeError(f"Cannot read parquet {path}: {e}") from e
        return _normalise_tx_df(df)

    # ── CSV / CSV.GZ ───────────────────────────────────────────────────────
    if suffix in (".csv", ".gz") or stem_suffix in (".csv",):
        try:
            df = pd.read_csv(path)
        except Exception as e:
            raise RuntimeError(f"Cannot read CSV {path}: {e}") from e
        return _normalise_tx_df(df)

    # ── zarr ───────────────────────────────────────────────────────────────
    if suffix in (".zarr", ".zip") or "zarr" in p:
        if not HAS_ZARR:
            raise RuntimeError("zarr not installed — cannot read zarr transcript store.")
        try:
            if p.endswith(".zip"):
                store = zarr.storage.ZipStore(p, mode="r")
            else:
                store = zarr.storage.DirectoryStore(p)
            z = zarr.open(store, mode="r")
        except Exception as e:
            raise RuntimeError(f"Cannot open zarr {path}: {e}") from e

        # Print keys to help diagnose unexpected structures
        top_keys = list(z.keys())
        print(f"  zarr top-level keys: {top_keys}")

        def _get(candidates):
            for k in candidates:
                if k in z:
                    return z[k][:]
            return None

        x_raw  = _get(_TX_X_COLS)
        y_raw  = _get(_TX_Y_COLS)
        g_raw  = _get(_TX_GENE_COLS)
        qv_raw = _get(_TX_QV_COLS)

        if x_raw is None or y_raw is None:
            # Maybe it's a cell×gene matrix zarr — not a transcript table
            raise RuntimeError(
                f"This zarr does not appear to be a transcript-coordinates store.\n"
                f"  Top-level keys: {top_keys}\n"
                f"  Expected keys like: x_location, y_location, feature_name\n"
                f"  Pass --transcripts pointing to transcripts.csv.gz or "
                f"transcripts.parquet instead."
            )

        d: Dict[str, np.ndarray] = {
            "x": np.asarray(x_raw, dtype=float),
            "y": np.asarray(y_raw, dtype=float),
        }
        if g_raw is not None:
            d["gene"] = np.array([
                g.decode() if isinstance(g, bytes) else str(g) for g in g_raw
            ])
        if qv_raw is not None:
            d["qv"] = np.asarray(qv_raw, dtype=float)

        return pd.DataFrame(d)

    raise RuntimeError(
        f"Unrecognised transcript file format: {path}\n"
        f"Supported: .csv, .csv.gz, .parquet, .zarr, .zarr.zip"
    )


def _normalise_tx_df(df: pd.DataFrame) -> pd.DataFrame:
    """Rename columns to canonical x / y / gene / qv."""
    df = df.copy()
    df.columns = df.columns.str.strip()

    rename = {}
    cols_lower = {c.lower(): c for c in df.columns}

    for cand in _TX_X_COLS:
        if cand in cols_lower:
            rename[cols_lower[cand]] = "x"
            break
    for cand in _TX_Y_COLS:
        if cand in cols_lower:
            rename[cols_lower[cand]] = "y"
            break
    for cand in _TX_GENE_COLS:
        if cand in cols_lower:
            rename[cols_lower[cand]] = "gene"
            break
    for cand in _TX_QV_COLS:
        if cand in cols_lower:
            rename[cols_lower[cand]] = "qv"
            break

    df = df.rename(columns=rename)

    missing = [c for c in ("x", "y", "gene") if c not in df.columns]
    if missing:
        raise RuntimeError(
            f"Transcript file missing required columns after normalisation: {missing}\n"
            f"  Found columns: {list(df.columns)}"
        )
    return df


def load_transcripts_for_gene(
    transcript_path: Path,
    gene: str,
    roi_xmin: float,
    roi_xmax: float,
    roi_ymin: float,
    roi_ymax: float,
    px_scale: float,
    min_qv: float = 20.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load x/y coordinates (in µm) of transcripts for *gene* within the bounding
    box.  Accepts CSV(.gz), Parquet, or zarr transcript tables.

    Returns (x_um, y_um) arrays.
    """
    df = _load_transcript_table(transcript_path)

    # Apply pixel scale
    df["x"] = df["x"] * px_scale
    df["y"] = df["y"] * px_scale

    # Quality filter
    if "qv" in df.columns:
        df = df[df["qv"] >= min_qv]

    # Bounding box pre-filter
    df = df[
        (df["x"] >= roi_xmin) & (df["x"] <= roi_xmax) &
        (df["y"] >= roi_ymin) & (df["y"] <= roi_ymax)
    ]

    # Gene filter (case-insensitive)
    df = df[df["gene"].str.upper() == gene.upper()]

    return df["x"].to_numpy(), df["y"].to_numpy()



# ══════════════════════════════════════════════════════════════════════════════
# Geometry helpers
# ══════════════════════════════════════════════════════════════════════════════

def find_nearest_cell(
    cell_polygons: Dict[str, np.ndarray],
    x_um: float,
    y_um: float,
) -> Tuple[str, float]:
    """
    Return (cell_id, distance_um) for the cell whose centroid is closest
    to the given (x_um, y_um) coordinate.
    """
    best_id, best_dist = None, np.inf
    for cid, pts in cell_polygons.items():
        cx, cy = pts[:, 0].mean(), pts[:, 1].mean()
        d = np.hypot(cx - x_um, cy - y_um)
        if d < best_dist:
            best_dist, best_id = d, cid
    return best_id, best_dist


def cell_centroid(poly_pts: np.ndarray) -> Tuple[float, float]:
    """Return (cx, cy) as the mean of boundary vertices."""
    return float(poly_pts[:, 0].mean()), float(poly_pts[:, 1].mean())


def points_in_circle(
    xs: np.ndarray, ys: np.ndarray,
    cx: float, cy: float, r: float,
) -> np.ndarray:
    """Boolean mask: which (xs, ys) fall within circle of radius r."""
    return (xs - cx) ** 2 + (ys - cy) ** 2 <= r ** 2


def build_union_polygon(polys_pts: List[np.ndarray]):
    """
    Build the shapely union of a list of polygon vertex arrays.
    Returns a shapely geometry or None if shapely is unavailable.
    """
    if not HAS_SHAPELY:
        return None
    shapes = []
    for pts in polys_pts:
        if len(pts) >= 3:
            try:
                shapes.append(Polygon(pts))
            except Exception:
                pass
    if not shapes:
        return None
    return unary_union(shapes)


def mask_points_to_union(
    xs: np.ndarray, ys: np.ndarray, union_geom
) -> np.ndarray:
    """
    Return boolean mask for points that fall inside *union_geom*.
    Falls back to all-True if shapely is unavailable or geometry is None.
    """
    if union_geom is None or not HAS_SHAPELY:
        return np.ones(len(xs), dtype=bool)
    mask = np.array([
        union_geom.contains(Point(x, y))
        for x, y in zip(xs, ys)
    ], dtype=bool)
    return mask


def _points_in_polygon_np(
    xs: np.ndarray, ys: np.ndarray, poly_pts: np.ndarray
) -> np.ndarray:
    """
    Vectorised ray-casting point-in-polygon for a single polygon.
    Used as a shapely-free fallback. poly_pts shape (M, 2).
    Returns a boolean mask of length len(xs).
    """
    n = len(poly_pts)
    if n < 3 or len(xs) == 0:
        return np.zeros(len(xs), dtype=bool)
    px = poly_pts[:, 0]
    py = poly_pts[:, 1]
    inside = np.zeros(len(xs), dtype=bool)
    j = n - 1
    for i in range(n):
        yi, yj = py[i], py[j]
        xi, xj = px[i], px[j]
        cond = ((yi > ys) != (yj > ys))
        # Avoid divide-by-zero when yj == yi (cond is False there so result ignored)
        denom = np.where(yj != yi, yj - yi, 1.0)
        x_intersect = (xj - xi) * (ys - yi) / denom + xi
        inside ^= cond & (xs < x_intersect)
        j = i
    return inside


def compute_per_cell_transcript_counts(
    cells: List[Tuple[str, np.ndarray]],
    tx: np.ndarray,
    ty: np.ndarray,
) -> np.ndarray:
    """
    Count transcripts falling inside each cell polygon.

    Parameters
    ----------
    cells : list of (cell_id, polygon_vertices) for the cells of interest.
            polygon_vertices is an (M, 2) array in the same units as tx/ty.
    tx, ty : transcript coordinates (already filtered to the relevant window if
             desired by the caller).

    Returns
    -------
    np.ndarray of length len(cells) — per-cell transcript counts (int).
    Order matches the input cell list.
    """
    counts = np.zeros(len(cells), dtype=int)
    if len(tx) == 0 or not cells:
        return counts

    tx = np.asarray(tx, dtype=float)
    ty = np.asarray(ty, dtype=float)

    for i, (_cid, pts) in enumerate(cells):
        if len(pts) < 3:
            continue
        # Cheap bbox prefilter
        xmin, ymin = pts[:, 0].min(), pts[:, 1].min()
        xmax, ymax = pts[:, 0].max(), pts[:, 1].max()
        cand = (tx >= xmin) & (tx <= xmax) & (ty >= ymin) & (ty <= ymax)
        if not cand.any():
            continue
        cand_x = tx[cand]
        cand_y = ty[cand]

        if HAS_SHAPELY:
            try:
                poly = Polygon(pts)
                inside = np.array([
                    poly.contains(Point(x, y)) for x, y in zip(cand_x, cand_y)
                ], dtype=bool)
            except Exception:
                inside = _points_in_polygon_np(cand_x, cand_y, pts)
        else:
            inside = _points_in_polygon_np(cand_x, cand_y, pts)

        counts[i] = int(inside.sum())

    return counts


def summarise_transcripts_per_cell(
    cells: List[Tuple[str, np.ndarray]],
    transcripts_per_gene: Dict[str, Tuple[np.ndarray, np.ndarray]],
    genes: List[str],
    cx: float,
    cy: float,
    radius_um: float,
) -> Dict[str, Dict[str, float]]:
    """
    For each gene, compute per-cell transcript count statistics across the
    given list of cells.

    Transcripts are first masked to the niche circle (cx, cy, radius_um) so the
    statistics reflect exactly what is visible in the inset. Each transcript is
    then attributed to at most one cell — a transcript that falls inside any
    shown cell polygon counts toward that cell's tally.

    Returns
    -------
    { gene: { 'mean': ..., 'median': ..., 'max': ...,
              'n_cells': int, 'n_transcripts_in_cells': int,
              'n_transcripts_in_circle': int } }
    """
    stats: Dict[str, Dict[str, float]] = {}
    n_cells = len(cells)
    for g in genes:
        xy = transcripts_per_gene.get(g)
        if xy is None:
            continue
        tx, ty = xy[0], xy[1]
        if len(tx) == 0 or n_cells == 0:
            stats[g] = {
                "mean": 0.0, "median": 0.0, "max": 0.0,
                "n_cells": n_cells,
                "n_transcripts_in_cells": 0,
                "n_transcripts_in_circle": 0,
            }
            continue

        circ_mask = points_in_circle(np.asarray(tx), np.asarray(ty),
                                     cx, cy, radius_um)
        tx_c = np.asarray(tx)[circ_mask]
        ty_c = np.asarray(ty)[circ_mask]

        counts = compute_per_cell_transcript_counts(cells, tx_c, ty_c)
        stats[g] = {
            "mean": float(counts.mean()) if n_cells > 0 else 0.0,
            "median": float(np.median(counts)) if n_cells > 0 else 0.0,
            "max": int(counts.max()) if n_cells > 0 else 0,
            "n_cells": n_cells,
            "n_transcripts_in_cells": int(counts.sum()),
            "n_transcripts_in_circle": int(circ_mask.sum()),
        }
    return stats


def format_transcript_stats(stats: Dict[str, Dict[str, float]]) -> str:
    """Format per-gene stats as a multi-line label for placing on a plot."""
    lines = []
    for g, s in stats.items():
        lines.append(
            f"{g}: mean={s['mean']:.1f}  med={s['median']:.1f}  "
            f"max={s['max']:.0f}  (n={s['n_cells']} cells)"
        )
    return "\n".join(lines)


# ══════════════════════════════════════════════════════════════════════════════
# Colour helpers
# ══════════════════════════════════════════════════════════════════════════════

def build_colour_map(
    cell_types: List[str],
    palette: Optional[List[str]] = None,
    theme: int = 1,
) -> Dict[str, str]:
    """
    Map each cell type to a hex colour.

    Priority:
      1. Explicit --palette (overrides everything)
      2. theme=2 → _THEME2_PALETTE lookup, fallback to _THEME1_COLOURS for unknowns
      3. theme=1 (default) → _THEME1_COLOURS in order
    """
    if palette and len(palette) >= len(cell_types):
        return dict(zip(cell_types, palette[: len(cell_types)]))

    if palette:
        warnings.warn(
            f"--palette has {len(palette)} entries but {len(cell_types)} "
            "cell types — using theme colours for the remainder."
        )

    result = {}
    theme1_idx = len(palette) if palette else 0  # offset if partial palette given

    for i, ct in enumerate(cell_types):
        if palette and i < len(palette):
            result[ct] = palette[i]
        elif theme == 2:
            # Use biology palette if the type is known, else fall back to theme1
            result[ct] = _THEME2_PALETTE.get(
                ct,
                _THEME1_COLOURS[theme1_idx % len(_THEME1_COLOURS)]
            )
            if ct not in _THEME2_PALETTE:
                theme1_idx += 1
        else:
            result[ct] = _THEME1_COLOURS[theme1_idx % len(_THEME1_COLOURS)]
            theme1_idx += 1

    return result


# ══════════════════════════════════════════════════════════════════════════════
# Core plotting
# ══════════════════════════════════════════════════════════════════════════════

def plot_niche(
    cell_polygons: Dict[str, np.ndarray],
    annotations: pd.DataFrame,
    focal_cell_id: str,
    radius_um: float,
    cell_types: List[str],
    colour_map: Dict[str, str],
    hide_other: bool = False,
    genes: Optional[List[str]] = None,
    transcripts_per_gene: Optional[Dict[str, Tuple[np.ndarray, np.ndarray]]] = None,
    gene_colours: Optional[Dict[str, str]] = None,
    figsize: Tuple[float, float] = (5.0, 5.0),
    dpi: int = 300,
    cell_alpha: float = 0.75,
    other_alpha: float = 0.25,
    show_circle: bool = True,
    transcript_dot_size: float = 2.0,
    cell_boundary_linewidth: float = 0.8,
    scalebar_fontsize: float = 9.0,
    title: Optional[str] = None,
) -> plt.Figure:
    """
    Create the niche spatial plot.

    Parameters
    ----------
    cell_polygons         : dict cell_id → (N,2) µm vertex array
    annotations           : DataFrame with columns cell_id, cell_type
    focal_cell_id         : string cell ID to centre the view on
    radius_um             : radius of circular ROI
    cell_types            : ordered list of cell types to colour
    colour_map            : { cell_type: hex_colour }
    hide_other            : if True, cells outside cell_types are hidden entirely
    genes                 : list of gene names to overlay (None = no overlay)
    transcripts_per_gene  : { gene: (x_um, y_um) } pre-loaded arrays
    gene_colours          : { gene: hex_colour } for transcript dots
    figsize / dpi         : figure dimensions
    cell_alpha            : alpha for filled cell bodies
    other_alpha           : alpha for grey 'other' cell outlines
    show_circle           : draw the ROI boundary circle
    transcript_dot_size   : scatter dot size
    cell_boundary_linewidth : thickness of the cell boundary line drawn when
                              transcripts are being overlaid (cells are rendered
                              hollow with this line width in the cell-type colour).
                              Ignored when no transcripts are plotted.
    title                 : optional axes title

    Returns
    -------
    matplotlib Figure
    """
    ann_map = dict(zip(annotations["cell_id"], annotations["cell_type"]))

    # ── focal cell centroid ────────────────────────────────────────────────
    if focal_cell_id not in cell_polygons:
        raise KeyError(
            f"Focal cell '{focal_cell_id}' not found in cell_boundaries. "
            f"Check your --focal_cell_id argument."
        )
    focal_pts = cell_polygons[focal_cell_id]
    cx, cy = cell_centroid(focal_pts)

    # ── collect cells whose centroid lies within radius ────────────────────
    cells_in_roi: List[Tuple[str, str, np.ndarray]] = []
    for cid, pts in cell_polygons.items():
        pcx, pcy = cell_centroid(pts)
        if (pcx - cx) ** 2 + (pcy - cy) ** 2 <= radius_um ** 2:
            ctype = ann_map.get(cid, "Unknown")
            cells_in_roi.append((cid, ctype, pts))

    if not cells_in_roi:
        warnings.warn("No cells found within the specified radius.")

    shown_types_set = set(cell_types)
    shown_polys: Dict[str, List[np.ndarray]] = {ct: [] for ct in cell_types}
    shown_cells: List[Tuple[str, np.ndarray]] = []  # for transcript stats
    other_polys: List[np.ndarray] = []

    for cid, ctype, pts in cells_in_roi:
        if ctype in shown_types_set:
            shown_polys[ctype].append(pts)
            shown_cells.append((cid, pts))
        else:
            other_polys.append(pts)

    # ── build union geometry for transcript masking ────────────────────────
    all_shown_pts = [pts for ct_list in shown_polys.values() for pts in ct_list]
    union_geom = build_union_polygon(all_shown_pts) if HAS_SHAPELY else None

    # ── figure setup ──────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.set_aspect("equal")
    ax.set_facecolor("white")
    ax.axis("off")

    # ── draw 'other' cells first (bottom layer) ───────────────────────────
    if not hide_other and other_polys:
        other_collection = PolyCollection(
            other_polys,
            facecolors=_OTHER_CELL_COLOUR,
            edgecolors=_OTHER_CELL_COLOUR,
            linewidths=_OTHER_CELL_LINEWIDTH,
            alpha=other_alpha,
            zorder=1,
        )
        ax.add_collection(other_collection)

    # When transcripts are being overlaid we render cells as hollow outlines so
    # transcript dots are clearly visible; otherwise the cells are filled.
    transcripts_active = bool(
        genes and transcripts_per_gene
        and any(
            transcripts_per_gene.get(g) is not None and len(transcripts_per_gene[g][0]) > 0
            for g in genes
        )
    )

    # ── draw shown cell types (coloured) ──────────────────────────────────
    for ct in cell_types:
        polys = shown_polys.get(ct, [])
        if not polys:
            continue
        fc = colour_map[ct]
        if transcripts_active:
            # Hollow cells: boundary in the cell-type colour, no fill
            coll = PolyCollection(
                polys,
                facecolors="none",
                edgecolors=fc,
                linewidths=cell_boundary_linewidth,
                zorder=2,
            )
        else:
            ec = mcolors.to_rgba(fc)
            ec = (*ec[:3], min(1.0, ec[3] * 1.5))
            coll = PolyCollection(
                polys,
                facecolors=fc,
                edgecolors=ec,
                linewidths=0.5,
                alpha=cell_alpha,
                zorder=2,
            )
        ax.add_collection(coll)

    # ── transcript overlays (one scatter per gene) ────────────────────────
    if genes and transcripts_per_gene:
        for gene in genes:
            xy = transcripts_per_gene.get(gene)
            if xy is None:
                continue
            tx, ty = xy
            if len(tx) == 0:
                continue

            # Circle mask
            circ_mask = points_in_circle(tx, ty, cx, cy, radius_um)
            tx, ty = tx[circ_mask], ty[circ_mask]

            # Cell-footprint mask
            if HAS_SHAPELY and union_geom is not None:
                tx, ty = tx[mask_points_to_union(tx, ty, union_geom)], \
                         ty[mask_points_to_union(tx, ty, union_geom)]
            elif not HAS_SHAPELY:
                warnings.warn(
                    "shapely not available — transcripts shown in full ROI, "
                    "not masked to cell boundaries."
                )

            dot_colour = (gene_colours or {}).get(gene, "#FFD700")
            ax.scatter(
                tx, ty,
                s=transcript_dot_size,
                c=dot_colour,
                linewidths=0.3,
                edgecolors="none",
                alpha=0.95,
                zorder=5,
                rasterized=True,
                label=gene,
            )

    # ── per-cell transcript statistics (over shown cells in this niche) ───
    if genes and transcripts_per_gene:
        stats = summarise_transcripts_per_cell(
            cells=shown_cells,
            transcripts_per_gene=transcripts_per_gene,
            genes=list(genes),
            cx=cx, cy=cy, radius_um=radius_um,
        )
        if stats:
            print(f"  Niche {focal_cell_id} — per-cell transcript stats "
                  f"({len(shown_cells)} shown cells):")
            for g, s in stats.items():
                print(f"    {g}: mean={s['mean']:.2f}  "
                      f"median={s['median']:.2f}  max={s['max']:.0f}  "
                      f"({s['n_transcripts_in_cells']} in cells / "
                      f"{s['n_transcripts_in_circle']} in circle)")

    # ── gene legend (bottom-left, small) ─────────────────────────────────
    if genes and transcripts_per_gene and gene_colours:
        gene_handles = [
            mpatches.Patch(facecolor=gene_colours[g], edgecolor="none", label=g)
            for g in genes if g in transcripts_per_gene
        ]
        if gene_handles:
            ax.legend(
                handles=gene_handles,
                loc="lower left",
                frameon=False,
                fontsize=7,
                handlelength=1.0,
                borderpad=0.3,
            )

    # ── ROI circle ────────────────────────────────────────────────────────
    if show_circle:
        circle = Circle(
            (cx, cy), radius_um,
            fill=False,
            edgecolor=_CIRCLE_COLOUR,
            linewidth=1.2,
            zorder=6,
        )
        ax.add_patch(circle)

    # ── axis limits ───────────────────────────────────────────────────────
    pad = radius_um * 0.05
    ax.set_xlim(cx - radius_um - pad, cx + radius_um + pad)
    ax.set_ylim(cy - radius_um - pad, cy + radius_um + pad)

    # ── scale bar ─────────────────────────────────────────────────────────
    sb_len = round(radius_um * 0.4 / 10) * 10
    sb_len = max(sb_len, 5)
    sb_x0 = cx - radius_um + radius_um * 0.07
    sb_y0 = cy - radius_um + radius_um * 0.07
    ax.plot(
        [sb_x0, sb_x0 + sb_len], [sb_y0, sb_y0],
        color="black", linewidth=1.5, solid_capstyle="butt", zorder=7,
    )
    ax.text(
        sb_x0 + sb_len / 2, sb_y0 + radius_um * 0.03,
        f"{sb_len} µm",
        ha="center", va="bottom", fontsize=scalebar_fontsize,
        color="black", zorder=7,
    )

    # ── title ─────────────────────────────────────────────────────────────
    if title:
        ax.set_title(title, fontsize=9, pad=4)
    else:
        gene_str = f" · {', '.join(genes)}" if genes else ""
        ax.set_title(
            f"Cell {focal_cell_id}{gene_str}  |  r = {radius_um} µm",
            fontsize=8,
            pad=4,
        )

    fig.tight_layout(pad=0.5)
    return fig


# ══════════════════════════════════════════════════════════════════════════════
# Composite layout: overview rectangle + per-focal-cell circular insets
# ══════════════════════════════════════════════════════════════════════════════

def _draw_cells_on_ax(
    ax: plt.Axes,
    cell_polygons: Dict[str, np.ndarray],
    ann_map: Dict[str, str],
    shown_types_set: set,
    colour_map: Dict[str, str],
    cell_types: List[str],
    xmin: float, xmax: float,
    ymin: float, ymax: float,
    hide_other: bool,
    cell_alpha: float,
    other_alpha: float,
    hollow_cells: bool = False,
    cell_boundary_linewidth: float = 0.8,
) -> Dict[str, List[np.ndarray]]:
    """
    Draw all cells whose centroid falls within [xmin,xmax] x [ymin,ymax] onto ax.
    Returns shown_polys dict for transcript masking.

    If hollow_cells is True, shown cells are drawn as boundary-only polygons
    in the cell-type colour using cell_boundary_linewidth, so transcript dots
    overlaid on top stand out more clearly.
    """
    shown_polys: Dict[str, List[np.ndarray]] = {ct: [] for ct in cell_types}
    other_polys: List[np.ndarray] = []

    for cid, pts in cell_polygons.items():
        pcx, pcy = pts[:, 0].mean(), pts[:, 1].mean()
        if not (xmin <= pcx <= xmax and ymin <= pcy <= ymax):
            continue
        ctype = ann_map.get(cid, "Unknown")
        if ctype in shown_types_set:
            shown_polys[ctype].append(pts)
        else:
            other_polys.append(pts)

    if not hide_other and other_polys:
        ax.add_collection(PolyCollection(
            other_polys,
            facecolors=_OTHER_CELL_COLOUR,
            edgecolors=_OTHER_CELL_COLOUR,
            linewidths=_OTHER_CELL_LINEWIDTH,
            alpha=other_alpha,
            zorder=1,
        ))

    for ct in cell_types:
        polys = shown_polys.get(ct, [])
        if not polys:
            continue
        fc = colour_map[ct]
        if hollow_cells:
            ax.add_collection(PolyCollection(
                polys,
                facecolors="none",
                edgecolors=fc,
                linewidths=cell_boundary_linewidth,
                zorder=2,
            ))
        else:
            ec = mcolors.to_rgba(fc)
            ec = (*ec[:3], min(1.0, ec[3] * 1.5))
            ax.add_collection(PolyCollection(
                polys,
                facecolors=fc,
                edgecolors=ec,
                linewidths=0.4,
                alpha=cell_alpha,
                zorder=2,
            ))

    return shown_polys


def _draw_transcripts_on_ax(
    ax: plt.Axes,
    genes: List[str],
    transcripts_per_gene: Dict[str, Tuple[np.ndarray, np.ndarray]],
    gene_colours: Dict[str, str],
    union_geom,
    cx: float, cy: float, radius_um: float,
    transcript_dot_size: float,
    clip_circle: bool = True,
) -> None:
    """Overlay transcript dots, masked to cell union and optionally a circle."""
    for gene in genes:
        xy = transcripts_per_gene.get(gene)
        if xy is None:
            continue
        tx, ty = xy[0].copy(), xy[1].copy()
        if len(tx) == 0:
            continue

        if clip_circle:
            circ_mask = points_in_circle(tx, ty, cx, cy, radius_um)
            tx, ty = tx[circ_mask], ty[circ_mask]

        if HAS_SHAPELY and union_geom is not None:
            m = mask_points_to_union(tx, ty, union_geom)
            tx, ty = tx[m], ty[m]

        dot_colour = gene_colours.get(gene, "#FFD700")
        ax.scatter(tx, ty, s=transcript_dot_size, c=dot_colour,
                   linewidths=0, alpha=0.95, zorder=5, rasterized=True)


def _add_scalebar(
    ax: plt.Axes,
    x0: float, y0: float,
    length_um: float,
    dy: float,
    fontsize: float = 6,
) -> None:
    ax.plot([x0, x0 + length_um], [y0, y0],
            color="black", linewidth=1.5, solid_capstyle="butt", zorder=10)
    ax.text(x0 + length_um / 2, y0 + dy,
            f"{int(length_um)} µm",
            ha="center", va="bottom", fontsize=fontsize, color="black", zorder=10)


def plot_overview_with_niches(
    cell_polygons: Dict[str, np.ndarray],
    annotations: pd.DataFrame,
    focal_cell_ids: List[str],
    niche_radius_um: float,
    overview_size_um: float,
    cell_types: List[str],
    colour_map: Dict[str, str],
    hide_other: bool = False,
    genes: Optional[List[str]] = None,
    transcripts_per_gene: Optional[Dict[str, Tuple[np.ndarray, np.ndarray]]] = None,
    gene_colours: Optional[Dict[str, str]] = None,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: int = 300,
    cell_alpha: float = 0.75,
    other_alpha: float = 0.25,
    transcript_dot_size: float = 2.0,
    cell_boundary_linewidth: float = 0.8,
    scalebar_fontsize: float = 9.0,
    title: Optional[str] = None,
    show_connectors: bool = True,
) -> plt.Figure:
    """
    Composite figure matching the reference panel style (B + C):

      Left panel  — large rectangular overview (~overview_size_um wide) showing
                    all cells in the region, with circular ROI outlines drawn
                    around each focal cell.
      Right panels — one circular inset per focal cell (niche_radius_um radius),
                     stacked vertically, with transcript dots if genes supplied.

    Parameters
    ----------
    focal_cell_ids   : one or more cell IDs; single ID → still uses composite layout
    niche_radius_um  : radius of each inset circle (200–500 µm recommended)
    overview_size_um : half-width of the rectangular overview (1500–2000 µm recommended)
    """
    ann_map = dict(zip(annotations["cell_id"], annotations["cell_type"]))
    shown_types_set = set(cell_types)
    n_niches = len(focal_cell_ids)

    # When transcripts are being overlaid we render cells as hollow outlines so
    # transcript dots are clearly visible; otherwise the cells are filled.
    transcripts_active = bool(
        genes and transcripts_per_gene
        and any(
            transcripts_per_gene.get(g) is not None and len(transcripts_per_gene[g][0]) > 0
            for g in genes
        )
    )

    # ── Compute focal centroids ────────────────────────────────────────────
    focal_centroids: List[Tuple[float, float]] = []
    for fid in focal_cell_ids:
        if fid not in cell_polygons:
            raise KeyError(f"Focal cell '{fid}' not found in cell_boundaries.")
        pts = cell_polygons[fid]
        focal_centroids.append(cell_centroid(pts))

    # ── Overview bounding box: centred on the mean of all focal cells ──────
    mean_cx = float(np.mean([c[0] for c in focal_centroids]))
    mean_cy = float(np.mean([c[1] for c in focal_centroids]))
    ov_half = overview_size_um / 2.0
    ov_xmin, ov_xmax = mean_cx - ov_half, mean_cx + ov_half
    ov_ymin, ov_ymax = mean_cy - ov_half, mean_cy + ov_half

    # ── Load transcripts over full overview bbox once ─────────────────────
    # (transcripts_per_gene already pre-loaded by caller over a wider region)

    # ── Figure layout ─────────────────────────────────────────────────────
    # Left col: overview (square aspect)
    # Right col: n_niches stacked insets
    if figsize is None:
        fw = 12.0
        fh = max(5.0, 3.5 * n_niches)
        figsize = (fw, fh)

    fig = plt.figure(figsize=figsize, dpi=dpi)
    fig.patch.set_facecolor("white")

    # GridSpec: 1 row × 2 cols, left col ~2× width of right col
    gs = gridspec.GridSpec(
        n_niches, 2,
        figure=fig,
        width_ratios=[2, 1],
        hspace=0.08,
        wspace=0.06,
        left=0.03, right=0.97, top=0.93, bottom=0.05,
    )

    # Overview ax spans all rows in left column
    ax_ov = fig.add_subplot(gs[:, 0])
    ax_ov.set_aspect("equal")
    ax_ov.set_facecolor("white")
    ax_ov.axis("off")

    # ── Draw overview cells ────────────────────────────────────────────────
    # The overview (left panel) always uses filled cell polygons coloured by
    # annotation, regardless of whether transcripts are being plotted. Hollow
    # outlines are reserved for the per-niche insets (right panel) where
    # transcript dots are actually overlaid.
    shown_polys_ov = _draw_cells_on_ax(
        ax_ov, cell_polygons, ann_map, shown_types_set, colour_map, cell_types,
        ov_xmin, ov_xmax, ov_ymin, ov_ymax,
        hide_other, cell_alpha, other_alpha,
        hollow_cells=False,
    )

    ax_ov.set_xlim(ov_xmin, ov_xmax)
    ax_ov.set_ylim(ov_ymin, ov_ymax)

    # Draw circle outlines for each focal cell on overview
    for i, (fid, (fcx, fcy)) in enumerate(zip(focal_cell_ids, focal_centroids)):
        circ = Circle((fcx, fcy), niche_radius_um,
                      fill=False, edgecolor="black", linewidth=1.2, zorder=6)
        ax_ov.add_patch(circ)

    # Overview scale bar
    sb_len_ov = round(overview_size_um * 0.15 / 100) * 100
    sb_len_ov = max(sb_len_ov, 100)
    _add_scalebar(
        ax_ov,
        x0=ov_xmin + overview_size_um * 0.04,
        y0=ov_ymin + overview_size_um * 0.03,
        length_um=sb_len_ov,
        dy=overview_size_um * 0.015,
        fontsize=scalebar_fontsize,
    )

    #if title:
        #ax_ov.set_title(title, fontsize=10, pad=6, loc="left")

    # ── Per-focal-cell inset axes (right column) ───────────────────────────
    ax_niches: List[plt.Axes] = []
    for i, (fid, (fcx, fcy)) in enumerate(zip(focal_cell_ids, focal_centroids)):
        ax_n = fig.add_subplot(gs[i, 1])
        ax_n.set_aspect("equal")
        ax_n.set_facecolor("white")
        ax_n.axis("off")

        r = niche_radius_um

        # Draw only cells whose centroid is within the circle (not just the bbox)
        shown_polys_n: Dict[str, List[np.ndarray]] = {ct: [] for ct in cell_types}
        shown_cells_n: List[Tuple[str, np.ndarray]] = []  # for transcript stats
        other_polys_n: List[np.ndarray] = []
        for cid, pts in cell_polygons.items():
            pcx, pcy = pts[:, 0].mean(), pts[:, 1].mean()
            if (pcx - fcx) ** 2 + (pcy - fcy) ** 2 > r ** 2:
                continue
            ctype = ann_map.get(cid, "Unknown")
            if ctype in shown_types_set:
                shown_polys_n[ctype].append(pts)
                shown_cells_n.append((cid, pts))
            else:
                other_polys_n.append(pts)

        if not hide_other and other_polys_n:
            ax_n.add_collection(PolyCollection(
                other_polys_n,
                facecolors=_OTHER_CELL_COLOUR, edgecolors=_OTHER_CELL_COLOUR,
                linewidths=_OTHER_CELL_LINEWIDTH, alpha=other_alpha, zorder=1,
            ))
        for ct in cell_types:
            polys = shown_polys_n.get(ct, [])
            if not polys:
                continue
            fc = colour_map[ct]
            if transcripts_active:
                ax_n.add_collection(PolyCollection(
                    polys, facecolors="none", edgecolors=fc,
                    linewidths=cell_boundary_linewidth, zorder=2,
                ))
            else:
                ec = mcolors.to_rgba(fc)
                ec = (*ec[:3], min(1.0, ec[3] * 1.5))
                ax_n.add_collection(PolyCollection(
                    polys, facecolors=fc, edgecolors=ec,
                    linewidths=0.4, alpha=cell_alpha, zorder=2,
                ))

        # Transcript overlay — mask to shown cells in this niche
        if genes and transcripts_per_gene:
            all_niche_pts = [pts for ct_list in shown_polys_n.values() for pts in ct_list]
            union_n = build_union_polygon(all_niche_pts) if HAS_SHAPELY else None
            _draw_transcripts_on_ax(
                ax_n, genes, transcripts_per_gene, gene_colours or {},
                union_n, fcx, fcy, r, transcript_dot_size, clip_circle=True,
            )

            # ── Per-cell transcript statistics over the shown window ──────
            stats = summarise_transcripts_per_cell(
                cells=shown_cells_n,
                transcripts_per_gene=transcripts_per_gene,
                genes=list(genes),
                cx=fcx, cy=fcy, radius_um=r,
            )
            if stats:
                # Console report only — stats are no longer drawn on the figure
                print(f"  Niche {fid} — per-cell transcript stats "
                      f"({len(shown_cells_n)} shown cells):")
                for g, s in stats.items():
                    print(f"    {g}: mean={s['mean']:.2f}  "
                          f"median={s['median']:.2f}  max={s['max']:.0f}  "
                          f"({s['n_transcripts_in_cells']} in cells / "
                          f"{s['n_transcripts_in_circle']} in circle)")

        # ROI circle
        circ_n = Circle((fcx, fcy), r,
                        fill=False, edgecolor="black", linewidth=1.2, zorder=6)
        ax_n.add_patch(circ_n)

        pad = r * 0.05
        ax_n.set_xlim(fcx - r - pad, fcx + r + pad)
        ax_n.set_ylim(fcy - r - pad, fcy + r + pad)

        # Scale bar
        sb_len_n = round(r * 0.4 / 10) * 10
        sb_len_n = max(sb_len_n, 10)
        _add_scalebar(
            ax_n,
            x0=fcx - r + r * 0.07,
            y0=fcy - r + r * 0.07,
            length_um=sb_len_n,
            dy=r * 0.03,
            fontsize=scalebar_fontsize,
        )

        ax_niches.append(ax_n)

    # ── Gene legend on last inset ──────────────────────────────────────────
    if genes and gene_colours:
        handles = [
            mpatches.Patch(facecolor=gene_colours.get(g, "#FFD700"),
                           edgecolor="none", label=g)
            for g in genes
        ]
        ax_niches[-1].legend(
            handles=handles,
            loc="lower left", frameon=False,
            fontsize=6, handlelength=0.9, borderpad=0.2,
        )

    return fig


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Plot spatial niches around focal cells.\n\n"
            "Always produces a composite figure:\n"
            "  Left  — large rectangular overview (overview_size_um wide) with\n"
            "          circle outlines marking each focal cell's niche.\n"
            "  Right — one circular inset per focal cell, stacked vertically,\n"
            "          with optional transcript dot overlays.\n\n"
            "Focal cell(s) specified with ONE of:\n"
            "  --focal_cell_id  (one or more cell IDs)\n"
            "  --x_um + --y_um  (single coordinate pair; nearest centroid used)"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--cell_boundaries", type=Path, required=True,
                   help="cell_boundaries.csv or cell_boundaries.csv.gz")
    p.add_argument("--annotations", type=Path, required=True,
                   help="CSV/TSV with columns: cell_id, cell_type")
    p.add_argument("--transcripts", type=Path, default=None,
                   help="Transcript file: transcripts.csv.gz, transcripts.parquet, "
                        "or a zarr transcript-coordinates store.")

    # ── focal point(s) ────────────────────────────────────────────────────
    p.add_argument("--focal_cell_id", type=str, nargs="+", default=None,
                   help="One or more cell IDs to centre niches around. "
                        "Multiple IDs trigger the composite overview+insets layout. "
                        "Mutually exclusive with --x_um / --y_um.")
    p.add_argument("--x_um", type=float, default=None,
                   help="X coordinate (µm) of niche centre (single point). "
                        "Nearest cell centroid is used. Pair with --y_um.")
    p.add_argument("--y_um", type=float, default=None,
                   help="Y coordinate (µm) of niche centre (pair with --x_um).")

    # ── region sizes ──────────────────────────────────────────────────────
    p.add_argument("--radius_um", type=float, default=300.0,
                   help="Radius of each circular niche inset in µm (default: 300). "
                        "Recommended 200–500 µm.")
    p.add_argument("--overview_size_um", type=float, default=1800.0,
                   help="Width/height of the rectangular overview panel in µm "
                        "(default: 1800). Recommended 1500–2000 µm. "
                        "Only used when multiple focal cells are given.")

    p.add_argument(
        "--cell_types", type=str, nargs="+", default=None,
        help="Cell types to colour and show. Others are grey outlines. "
             "Example: --cell_types Malignant myCAF TAM",
    )
    p.add_argument("--hide_other", action="store_true",
                   help="Hide cells not in --cell_types entirely")
    p.add_argument("--genes", type=str, nargs="+", default=None,
                   help="Genes to overlay as transcript dots (requires --transcripts). "
                        "Example: --genes ERBB2 MUC5AC")
    p.add_argument("--gene_colors", type=str, default=None,
                   help="Comma-separated HEX colours for --genes in order. "
                        "Auto-assigned if omitted.")
    p.add_argument("--min_qv", type=float, default=20.0,
                   help="Minimum transcript Phred QV (default: 20)")
    p.add_argument("--palette", type=str, default=None,
                   help="Comma-separated HEX colours for --cell_types. "
                        "Overrides --color_theme.")
    p.add_argument("--color_theme", type=int, default=1, choices=[1, 2],
                   help="Cell colour theme: 1=default high-contrast, "
                        "2=biology-informed (default: 1)")
    p.add_argument("--figsize", type=str, default=None,
                   help="Figure size as W,H inches (auto-sized if omitted)")
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--outfile", type=Path, default=Path("niche_plot.png"))
    p.add_argument("--px_scale", type=float, default=None)
    p.add_argument("--transcript_dot_size", type=float, default=2.0)
    p.add_argument("--cell_boundary_linewidth", type=float, default=0.8,
                   help="Thickness of the cell-boundary line used when "
                        "transcripts are overlaid (cells become hollow with this "
                        "outline in the cell-type colour). Default: 0.8. "
                        "Ignored when no --genes are plotted.")
    p.add_argument("--scalebar_fontsize", type=float, default=9.0,
                   help="Font size for the scale-bar label (e.g. '100 µm') "
                        "on both the overview and the per-niche insets. "
                        "Default: 9. Bump up to 11–14 for posters or talks.")
    p.add_argument("--cell_alpha", type=float, default=0.75)
    p.add_argument("--other_alpha", type=float, default=0.25)
    p.add_argument("--show_circle", type=lambda x: x.lower() != "false",
                   default=True,
                   help="Draw ROI circle on single-cell plot (default: true)")
    p.add_argument("--no_connectors", action="store_true",
                   help="Suppress dashed connector lines in composite layout")
    p.add_argument("--title", type=str, default=None)
    return p.parse_args()


def _resolve_focal_cells(
    args,
    cell_polygons: Dict[str, np.ndarray],
) -> List[str]:
    """
    Return a list of focal cell IDs from CLI args.
    Accepts --focal_cell_id (one or more) OR --x_um + --y_um (single point).
    """
    has_ids = args.focal_cell_id is not None and len(args.focal_cell_id) > 0
    has_xy  = args.x_um is not None and args.y_um is not None
    one_xy  = (args.x_um is None) != (args.y_um is None)

    if one_xy:
        sys.exit("ERROR: --x_um and --y_um must both be supplied together.")
    if has_ids and has_xy:
        sys.exit("ERROR: supply either --focal_cell_id or --x_um/--y_um, not both.")
    if not has_ids and not has_xy:
        sys.exit("ERROR: supply --focal_cell_id (one or more) or --x_um + --y_um.")

    if has_ids:
        missing = [fid for fid in args.focal_cell_id if fid not in cell_polygons]
        if missing:
            sys.exit(f"ERROR: focal cell IDs not found in boundaries: {missing}")
        return list(args.focal_cell_id)

    # Single x/y coordinate → nearest cell
    fid, dist = find_nearest_cell(cell_polygons, args.x_um, args.y_um)
    print(f"  Nearest cell to ({args.x_um}, {args.y_um}) µm: "
          f"cell_id={fid}  (distance={dist:.2f} µm)")
    return [fid]


def main():
    args = parse_args()

    # ── Load data ─────────────────────────────────────────────────────────
    print(f"Loading cell boundaries: {args.cell_boundaries}")
    boundaries_df = load_cell_boundaries(args.cell_boundaries)
    print(f"  {len(boundaries_df):,} boundary rows, "
          f"{boundaries_df['cell_id'].nunique():,} unique cells")

    print(f"Loading annotations: {args.annotations}")
    annotations = load_annotations(args.annotations)
    print(f"  {len(annotations):,} annotated cells  "
          f"({annotations['cell_type'].nunique()} types)")

    px_scale = args.px_scale or detect_px_scale(boundaries_df)
    print(f"  Coordinate scale: {px_scale} µm/unit")

    print("Building cell polygons …")
    cell_polygons = build_cell_polygons(boundaries_df, px_scale)

    focal_cell_ids = _resolve_focal_cells(args, cell_polygons)
    print(f"  Focal cells ({len(focal_cell_ids)}): {focal_cell_ids}")

    # ── Cell types & colours ──────────────────────────────────────────────
    if args.cell_types is None:
        cell_types = sorted(annotations["cell_type"].unique().tolist())
        print(f"  Showing all {len(cell_types)} annotated types")
    else:
        cell_types = args.cell_types

    palette = [c.strip() for c in args.palette.split(",")] if args.palette else None
    colour_map = build_colour_map(cell_types, palette, theme=args.color_theme)
    print(f"  Colour theme: {args.color_theme}")

    # ── Figure size ───────────────────────────────────────────────────────
    if args.figsize:
        try:
            fw, fh = [float(v) for v in args.figsize.split(",")]
            figsize: Optional[Tuple[float, float]] = (fw, fh)
        except ValueError:
            figsize = None
    else:
        figsize = None

    # ── Load transcripts ──────────────────────────────────────────────────
    transcripts_per_gene: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    gene_colours: Dict[str, str] = {}

    if args.genes:
        if args.transcripts is None:
            warnings.warn("--genes given but --transcripts not supplied.")
        else:
            raw_gc = [c.strip() for c in args.gene_colors.split(",")] \
                if args.gene_colors else []
            for i, g in enumerate(args.genes):
                gene_colours[g] = raw_gc[i] if i < len(raw_gc) \
                    else _DEFAULT_GENE_COLOURS[i % len(_DEFAULT_GENE_COLOURS)]

            # Load over the full region encompassing all focal cells + radius
            all_cx = [cell_centroid(cell_polygons[fid])[0] for fid in focal_cell_ids]
            all_cy = [cell_centroid(cell_polygons[fid])[1] for fid in focal_cell_ids]
            # Use the larger of overview or niche radius as the load window
            load_r = max(args.overview_size_um / 2.0, args.radius_um)
            roi_xmin = min(all_cx) - load_r
            roi_xmax = max(all_cx) + load_r
            roi_ymin = min(all_cy) - load_r
            roi_ymax = max(all_cy) + load_r

            for g in args.genes:
                print(f"Loading transcripts for '{g}' …")
                tx, ty = load_transcripts_for_gene(
                    transcript_path=args.transcripts,
                    gene=g,
                    roi_xmin=roi_xmin, roi_xmax=roi_xmax,
                    roi_ymin=roi_ymin, roi_ymax=roi_ymax,
                    px_scale=px_scale,
                    min_qv=args.min_qv,
                )
                print(f"  {len(tx):,} transcripts  [{gene_colours[g]}]")
                transcripts_per_gene[g] = (tx, ty)

    # ── Plot ──────────────────────────────────────────────────────────────
    n = len(focal_cell_ids)
    print(f"\nComposite layout: overview {args.overview_size_um} µm  +  "
          f"{n} inset(s) @ r={args.radius_um} µm")

    # Use outfile stem as title if no explicit title given
    plot_title = args.title if args.title else args.outfile.stem

    fig = plot_overview_with_niches(
        cell_polygons=cell_polygons,
        annotations=annotations,
        focal_cell_ids=focal_cell_ids,
        niche_radius_um=args.radius_um,
        overview_size_um=args.overview_size_um,
        cell_types=cell_types,
        colour_map=colour_map,
        hide_other=args.hide_other,
        genes=args.genes,
        transcripts_per_gene=transcripts_per_gene or None,
        gene_colours=gene_colours or None,
        figsize=figsize,
        dpi=args.dpi,
        cell_alpha=args.cell_alpha,
        other_alpha=args.other_alpha,
        transcript_dot_size=args.transcript_dot_size,
        cell_boundary_linewidth=args.cell_boundary_linewidth,
        scalebar_fontsize=args.scalebar_fontsize,
        title=plot_title,
        show_connectors=False,
    )

    args.outfile.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.outfile, bbox_inches="tight")
    print(f"✓  Saved: {args.outfile}")
    plt.close(fig)


if __name__ == "__main__":
    main()

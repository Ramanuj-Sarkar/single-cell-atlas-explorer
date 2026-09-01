#!/usr/bin/env python
"""Single-cell RNA-seq analysis of the 10x Genomics "3k PBMCs from a Healthy Donor" dataset.

This module implements the canonical Scanpy workflow for a raw 10x Genomics gene
expression matrix, packaged so that the *same* functions drive both

  * the standalone command-line script  (``pbmc3k_scanpy.py --step ...``), and
  * the narrated notebook                (``pbmc3k_scanpy.ipynb``),

and so that each step can be executed independently (which is how the Nextflow
pipeline orchestrates them) or as one shot with ``--step all``.

Pipeline steps
--------------
1.  ``qc``        Quality control: gene/cell filtering + mitochondrial read fraction
2.  ``normalize`` Normalization (library-size scaling + log1p) + HVG selection
3.  ``reduce``    Dimensionality reduction: PCA -> kNN graph -> UMAP
4.  ``cluster``   Leiden clustering + per-cluster marker genes (Wilcoxon)
5.  ``annotate``  Cell-type annotation from curated marker panels + proportions

Every step also writes MultiQC-compatible custom content (``*_mqc.json``) so the
QC metrics and cell-type composition can be aggregated into a single HTML report.

Example
-------
    python pbmc3k_scanpy.py --input ../../data/pbmc3k_raw.h5ad --outdir ../../results --step all
    python pbmc3k_scanpy.py --input filtered.h5ad --outdir out --step cluster
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scanpy import logging as logg

# --------------------------------------------------------------------------- #
#  Tunable parameters (shared by script and notebook)
# --------------------------------------------------------------------------- #

QC_MIN_GENES = 200        # keep cells with >= this many detected genes
QC_MAX_GENES = 2500       # keep cells with <= this many detected genes (doublets)
QC_MIN_CELLS = 3          # keep genes detected in >= this many cells
QC_MAX_MT_PCT = 5.0       # drop cells with > this % of reads mapping to MT genes
NORM_TARGET_SUM = 1e4     # library-size scaling target (CPM-like)
HVG_KWARGS = dict(min_mean=0.0125, max_mean=3, min_disp=0.5)  # Seurat-style HVG
PCA_N_COMPONENTS = 50     # number of PCs to compute
NEIGHBORS_KWARGS = dict(n_neighbors=10, n_pcs=40)
LEIDEN_RESOLUTION = 0.8  # recovers all 8 canonical PBMC populations incl. DCs (see notebook)
N_MARKER_GENES = 25       # top markers reported per cluster (Wilcoxon)

# Curated marker panels (classic PBMC biology; see the notebook for references).
# Genes are scored against the log-normalized (``.raw``) matrix.
MARKER_PANELS = {
    "CD4 T cells": ["CD3D", "CD3E", "CD4", "IL7R", "LEF1", "TCF7"],
    "CD8 T cells": ["CD3D", "CD3E", "CD8A", "CD8B", "GZMK", "CCL5"],
    "NK cells": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCAM1", "GZMB"],
    "B cells": ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1"],
    "CD14+ monocytes": ["CD14", "LYZ", "S100A8", "S100A9", "CST3", "FCN1"],
    "FCGR3A+ monocytes": ["FCGR3A", "MS4A7", "LST1", "IFITM3", "CTSS"],
    "Dendritic cells": ["FCER1A", "CST3", "CLEC10A", "HLA-DPA1", "HLA-DPB1", "CD1C"],
    "Platelets": ["PPBP", "PF4", "GNG11", "SDPR", "NRGN"],
}

DEFAULT_MARKER_GENES_FOR_PLOTS = [
    "CD3D", "CD8A", "NKG7", "MS4A1", "CD14", "FCGR3A", "FCER1A", "PPBP",
]

SAMPLE_NAME = "pbmc3k"  # dataset id shown in MultiQC


# --------------------------------------------------------------------------- #
#  Small helpers
# --------------------------------------------------------------------------- #

def ensure_dirs(outdir: Path) -> None:
    """Create the output directory layout (objects / tables / figures)."""
    for sub in ("objects", "tables", "figures"):
        (Path(outdir) / sub).mkdir(parents=True, exist_ok=True)


def _savefig(fig, path: Path, dpi: int = 150) -> None:
    """Save the current matplotlib figure and close it (Agg-safe)."""
    import matplotlib.pyplot as plt

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    logg.info(f"  wrote {path}")


def _write_mqc_json(path: Path, payload: dict) -> None:
    """Write a MultiQC custom-content JSON file (auto-detected via the *_mqc.json suffix)."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2)
    logg.info(f"  wrote {path}")


def _write_tsv(df: pd.DataFrame, path: Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    logg.info(f"  wrote {path}")


# --------------------------------------------------------------------------- #
#  Step 1 - Quality control
# --------------------------------------------------------------------------- #

def run_qc(
    adata: sc.AnnData,
    outdir: Path,
    *,
    sample: str = SAMPLE_NAME,
    min_genes: int = QC_MIN_GENES,
    max_genes: int = QC_MAX_GENES,
    min_cells: int = QC_MIN_CELLS,
    max_mt_pct: float = QC_MAX_MT_PCT,
) -> sc.AnnData:
    """Filter genes/cells, flag mitochondrial reads, and report QC metrics."""
    outdir = Path(outdir)
    ensure_dirs(outdir)

    n_cells_raw = adata.n_obs

    # --- mitochondrial fraction ------------------------------------------- #
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    # --- QC figures (before filtering) ------------------------------------ #
    import matplotlib.pyplot as plt

    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False,
    )
    _savefig(plt.gcf(), outdir / "figures" / "qc_violins_raw.png")

    # --- filtering -------------------------------------------------------- #
    logg.info("Filtering cells and genes")
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, max_genes=max_genes)
    adata = adata[adata.obs["pct_counts_mt"] < max_mt_pct, :].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells)

    n_cells_out = adata.n_obs
    n_genes_out = adata.n_vars

    # --- QC figures (after filtering) ------------------------------------- #
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False,
    )
    _savefig(plt.gcf(), outdir / "figures" / "qc_violins_filtered.png")

    # --- tables ----------------------------------------------------------- #
    metrics = {
        "n_cells_raw": n_cells_raw,
        "n_cells_after_qc": n_cells_out,
        "n_genes_after_qc": n_genes_out,
        "n_cells_removed": n_cells_raw - n_cells_out,
        "median_genes_per_cell": float(adata.obs["n_genes_by_counts"].median()),
        "median_counts_per_cell": float(adata.obs["total_counts"].median()),
        "median_pct_counts_mt": float(adata.obs["pct_counts_mt"].median()),
        "qc_min_genes": min_genes,
        "qc_max_genes": max_genes,
        "qc_max_mt_pct": max_mt_pct,
    }
    _write_tsv(pd.DataFrame([metrics]), outdir / "tables" / "qc_metrics.tsv")

    _write_mqc_json(
        outdir / "tables" / "multiqc_scanpy_qc_mqc.json",
        {
            "id": "scanpy_qc_summary",
            "section_name": "Scanpy QC summary",
            "description": (
                "Quality-control metrics for the 10x Genomics PBMC3k dataset "
                "(before/after cell & gene filtering)."
            ),
            "plot_type": "table",
            "pconfig": {
                "id": "scanpy_qc_summary_plot",
                "title": "Scanpy QC summary",
                "sortRows": False,
            },
            "data": {sample: metrics},
        },
    )

    adata.write(outdir / "objects" / "qc.h5ad")
    logg.info(
        f"QC done: {n_cells_raw} -> {n_cells_out} cells, {n_genes_out} genes "
        f"({100 * (n_cells_raw - n_cells_out) / n_cells_raw:.1f}% of cells removed)"
    )
    return adata


# --------------------------------------------------------------------------- #
#  Step 2 - Normalization + highly variable genes
# --------------------------------------------------------------------------- #

def run_normalize(
    adata: sc.AnnData,
    outdir: Path,
    *,
    target_sum: float = NORM_TARGET_SUM,
    hvg_kwargs: dict | None = None,
) -> sc.AnnData:
    """Normalize to a common library size, log-transform, select HVGs."""
    outdir = Path(outdir)
    ensure_dirs(outdir)

    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.raw = adata  # keep full log-normalized matrix for marker ranking

    sc.pp.highly_variable_genes(adata, **(hvg_kwargs or HVG_KWARGS))
    n_hvg = int(adata.var["highly_variable"].sum())

    adata_hvg = adata[:, adata.var["highly_variable"]].copy()

    hvg_table = (
        adata.var[["highly_variable", "means", "dispersions_norm"]]
        .sort_values("dispersions_norm", ascending=False)
        .reset_index()
        .rename(columns={"index": "gene"})
    )
    _write_tsv(hvg_table, outdir / "tables" / "hvg_genes.tsv")

    # HVG dispersion plot
    import matplotlib.pyplot as plt

    sc.pl.highly_variable_genes(adata, show=False)
    _savefig(plt.gcf(), outdir / "figures" / "hvg_dispersion.png")

    adata_hvg.write(outdir / "objects" / "normalize.h5ad")
    logg.info(f"Normalization done: {n_hvg} highly variable genes selected")
    return adata_hvg


# --------------------------------------------------------------------------- #
#  Step 3 - Dimensionality reduction (PCA -> kNN graph -> UMAP)
# --------------------------------------------------------------------------- #

def run_reduce(
    adata: sc.AnnData,
    outdir: Path,
    *,
    n_comps: int = PCA_N_COMPONENTS,
    neighbors_kwargs: dict | None = None,
) -> sc.AnnData:
    """PCA, neighborhood graph and UMAP embedding."""
    outdir = Path(outdir)
    ensure_dirs(outdir)

    sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack")
    sc.pp.neighbors(adata, **(neighbors_kwargs or NEIGHBORS_KWARGS))
    sc.tl.umap(adata)

    import matplotlib.pyplot as plt

    # variance explained
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(
        np.arange(1, len(adata.uns["pca"]["variance_ratio"]) + 1),
        adata.uns["pca"]["variance_ratio"],
        marker="o",
        markersize=3,
        lw=1,
    )
    ax.set_xlabel("PC")
    ax.set_ylabel("Variance explained")
    ax.set_title("PCA variance explained")
    _savefig(fig, outdir / "figures" / "pca_variance.png")

    # UMAP colored by total counts (a sanity check for technical artifacts)
    fig, ax = plt.subplots(figsize=(6.5, 5))
    sc.pl.umap(adata, color="total_counts", show=False, ax=ax, color_map="viridis")
    _savefig(fig, outdir / "figures" / "umap_total_counts.png")

    pca_var = pd.DataFrame(
        {
            "pc": np.arange(1, len(adata.uns["pca"]["variance_ratio"]) + 1),
            "variance_ratio": adata.uns["pca"]["variance_ratio"],
        }
    )
    _write_tsv(pca_var, outdir / "tables" / "pca_variance_ratio.tsv")

    adata.write(outdir / "objects" / "reduce.h5ad")
    logg.info(
        f"Dimensionality reduction done: {n_comps} PCs, "
        f"{adata.uns['neighbors']['params']['n_neighbors']} neighbors"
    )
    return adata


# --------------------------------------------------------------------------- #
#  Step 4 - Clustering + marker genes
# --------------------------------------------------------------------------- #

def run_cluster(
    adata: sc.AnnData,
    outdir: Path,
    *,
    resolution: float = LEIDEN_RESOLUTION,
    n_marker_genes: int = N_MARKER_GENES,
) -> sc.AnnData:
    """Leiden clustering on the kNN graph + Wilcoxon marker-gene ranking."""
    outdir = Path(outdir)
    ensure_dirs(outdir)

    sc.tl.leiden(adata, resolution=resolution, flavor="igraph", directed=False)
    n_clusters = adata.obs["leiden"].nunique()

    sc.tl.rank_genes_groups(
        adata, groupby="leiden", method="wilcoxon", use_raw=True, pts=True
    )
    marker_df = (
        sc.get.rank_genes_groups_df(adata, group=None)
        .query("pvals_adj < 0.05")
        .groupby("group", sort=True)
        .head(n_marker_genes)
    )
    _write_tsv(marker_df, outdir / "tables" / "marker_genes.tsv")

    import matplotlib.pyplot as plt

    # UMAP colored by Leiden cluster
    fig, ax = plt.subplots(figsize=(6.5, 5))
    sc.pl.umap(
        adata, color="leiden", show=False, ax=ax, palette="tab20", legend_loc="on data"
    )
    _savefig(fig, outdir / "figures" / "umap_leiden.png")

    # heatmap of the top markers per cluster
    top = marker_df.groupby("group").head(5)["names"].tolist()
    sc.pl.heatmap(
        adata, var_names=top, groupby="leiden", use_raw=True, swap_axes=True, show=False
    )
    _savefig(plt.gcf(), outdir / "figures" / "marker_heatmap_top5.png")

    # MultiQC table: top 5 markers per cluster
    top5 = marker_df.groupby("group").head(5)
    mqc_data = {}
    for cluster, sub in top5.groupby("group", sort=True):
        mqc_data[str(cluster)] = {
            f"marker_{i + 1}": g for i, g in enumerate(sub["names"].tolist())
        }
    _write_mqc_json(
        outdir / "tables" / "multiqc_scanpy_markers_mqc.json",
        {
            "id": "scanpy_marker_genes",
            "section_name": "Top marker genes per cluster",
            "description": (
                "Top 5 differentially expressed genes per Leiden cluster "
                "(Wilcoxon rank-sum, adjusted p < 0.05)."
            ),
            "plot_type": "table",
            "pconfig": {
                "id": "scanpy_marker_genes_plot",
                "title": "Top marker genes per cluster",
            },
            "data": mqc_data,
        },
    )

    adata.write(outdir / "objects" / "cluster.h5ad")
    logg.info(f"Clustering done: {n_clusters} Leiden clusters at resolution {resolution}")
    return adata


# --------------------------------------------------------------------------- #
#  Step 5 - Cell-type annotation
# --------------------------------------------------------------------------- #

def score_marker_panels(adata: sc.AnnData, panels: dict) -> pd.DataFrame:
    """Per-cluster mean log-normalized expression of each marker panel.

    Returns a DataFrame indexed by cluster with one column per panel. The panel
    with the highest score is the proposed cell type.
    """
    expr = adata.raw.to_adata() if adata.raw is not None else adata
    scores = {}
    for celltype, genes in panels.items():
        present = [g for g in genes if g in expr.var_names]
        if not present:
            scores[celltype] = np.nan
            continue
        sub = expr[:, present]
        X = sub.X
        if hasattr(X, "toarray"):  # sparse
            X = X.toarray()
        mean_expr = pd.Series(np.asarray(X).mean(axis=1), index=expr.obs_names)
        scores[celltype] = mean_expr.groupby(expr.obs["leiden"]).mean()
    return pd.DataFrame(scores)


def run_annotate(
    adata: sc.AnnData,
    outdir: Path,
    *,
    sample: str = SAMPLE_NAME,
    panels: dict | None = None,
) -> sc.AnnData:
    """Assign cell-type labels from curated marker panels and summarise composition."""
    outdir = Path(outdir)
    ensure_dirs(outdir)

    panels = panels or MARKER_PANELS
    panel_scores = score_marker_panels(adata, panels)

    # For each cluster, pick the best-scoring panel (guard against NaN panels)
    best = panel_scores.idxmax(axis=1)
    second = panel_scores.apply(lambda r: r.nlargest(2).index[-1], axis=1)

    annotation = pd.DataFrame(
        {
            "cluster": panel_scores.index.astype(str),
            "cell_type": best.values,
            "top_score": panel_scores.max(axis=1).round(3).values,
            "second_best": second.values,
            "score_gap": (
                panel_scores.max(axis=1) - panel_scores.apply(lambda r: r.nlargest(2).iloc[-1], axis=1)
            ).round(3).values,
        }
    ).sort_values("cluster")

    celltype_map = dict(zip(annotation["cluster"], annotation["cell_type"]))
    adata.obs["cell_type"] = adata.obs["leiden"].astype(str).map(celltype_map).astype("category")

    # composition table
    counts = adata.obs["cell_type"].value_counts()
    composition = pd.DataFrame(
        {
            "cell_type": counts.index,
            "n_cells": counts.values,
            "fraction": (counts.values / adata.n_obs).round(4),
        }
    ).sort_values("n_cells", ascending=False)
    _write_tsv(composition, outdir / "tables" / "celltype_proportions.tsv")
    _write_tsv(annotation, outdir / "tables" / "celltype_annotation.tsv")

    import matplotlib.pyplot as plt

    # UMAP colored by annotated cell type
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    sc.pl.umap(
        adata,
        color="cell_type",
        show=False,
        ax=ax,
        palette="tab20",
        legend_loc="on data",
        frameon=False,
    )
    _savefig(fig, outdir / "figures" / "umap_celltype.png")

    # dotplot of canonical markers per cell type
    genes = DEFAULT_MARKER_GENES_FOR_PLOTS
    sc.pl.dotplot(
        adata,
        var_names=genes,
        groupby="cell_type",
        use_raw=True,
        show=False,
        standard_scale="var",
    )
    _savefig(plt.gcf(), outdir / "figures" / "dotplot_markers.png")

    # MultiQC barplot: cell-type composition
    _write_mqc_json(
        outdir / "tables" / "multiqc_scanpy_celltypes_mqc.json",
        {
            "id": "scanpy_celltype_proportions",
            "section_name": "Cell-type proportions",
            "description": (
                "Fraction of cells assigned to each annotated cell type "
                "(curated marker-panel scoring, Leiden clusters)."
            ),
            "plot_type": "barplot",
            "pconfig": {
                "id": "scanpy_celltype_proportions_plot",
                "title": "Cell-type composition",
                "cpswitch": False,
                "ylab": "Fraction of cells",
            },
            "data": {
                sample: {
                    row.cell_type: float(row.fraction)
                    for row in composition.itertuples()
                }
            },
        },
    )

    adata.write(outdir / "objects" / "annotate.h5ad")
    logg.info(
        "Annotation done: "
        + ", ".join(f"{ct} n={n}" for ct, n in counts.items())
    )
    return adata


# --------------------------------------------------------------------------- #
#  Full pipeline + CLI
# --------------------------------------------------------------------------- #

STEPS = ("qc", "normalize", "reduce", "cluster", "annotate")


def run_pipeline(
    input_h5ad: str | Path,
    outdir: str | Path,
    steps: tuple[str, ...] = STEPS,
) -> sc.AnnData:
    """Run the requested steps sequentially over the input h5ad file."""
    adata = sc.read_h5ad(input_h5ad)
    for step in steps:
        logg.info(f"===== Step: {step} =====")
        if step == "qc":
            adata = run_qc(adata, outdir)
        elif step == "normalize":
            adata = run_normalize(adata, outdir)
        elif step == "reduce":
            adata = run_reduce(adata, outdir)
        elif step == "cluster":
            adata = run_cluster(adata, outdir)
        elif step == "annotate":
            adata = run_annotate(adata, outdir)
        else:
            raise ValueError(f"Unknown step: {step!r}")
    return adata


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", required=True, help="Input raw h5ad file")
    parser.add_argument("--outdir", default="results", help="Output directory")
    parser.add_argument(
        "--step",
        choices=("all", *STEPS),
        default="all",
        help="Run a single step (requires the previous step's h5ad as --input), or 'all'.",
    )
    parser.add_argument(
        "--resolution",
        type=float,
        default=LEIDEN_RESOLUTION,
        help="Leiden clustering resolution (used by the 'cluster' step).",
    )
    parser.add_argument(
        "--verbosity",
        type=int,
        default=3,
        help="Scanpy verbosity level (0=error ... 4=debug)",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)
    import matplotlib

    matplotlib.use("Agg")  # headless figure rendering
    sc.settings.verbosity = args.verbosity
    sc.set_figure_params(dpi=120, frameon=False)

    steps = STEPS if args.step == "all" else (args.step,)
    adata = run_pipeline(args.input, args.outdir, steps)

    if args.step in ("all", "cluster"):
        logg.info(
            f"Leiden clustering used resolution={args.resolution}; "
            f"{adata.obs['leiden'].nunique()} clusters found"
        )
    logg.info("Pipeline finished successfully.")


if __name__ == "__main__":
    main()

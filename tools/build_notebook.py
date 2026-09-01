#!/usr/bin/env python
"""Build and execute the narrated notebook for the PBMC3k Scanpy analysis.

The notebook mirrors ``analysis/pbmc3k_scanpy.py``: every code cell calls a
function from that module, so the script and the notebook can never drift
apart. Execution runs inside the ``analysis/`` directory (so imports and the
relative ``../results`` output paths behave like the standalone script run from
the repository root) and the executed notebook, including all figure outputs,
is written back to ``analysis/pbmc3k_scanpy.ipynb``.

Usage
-----
    .venv/bin/python tools/build_notebook.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import nbformat as nbf

ROOT = Path(__file__).resolve().parent.parent
ANALYSIS = ROOT / "analysis"
RESULTS = ROOT / "results"

nb = nbf.v4.new_notebook()
nb.metadata["kernelspec"] = {
    "name": "python3",
    "display_name": "Python 3",
    "language": "python",
}
nb.metadata["language_info"] = {"name": "python", "version": "3"}

cells: list = []


def md(source: str) -> None:
    cells.append(nbf.v4.new_markdown_cell(source))


def code(source: str) -> None:
    cells.append(nbf.v4.new_code_cell(source))


# --------------------------------------------------------------------------- #
#  Markdown narrative
# --------------------------------------------------------------------------- #

md(
    """# Single-cell RNA-seq analysis of 3k PBMCs from a Healthy Donor

A reproducible Scanpy workflow on the classic [10x Genomics PBMC3k dataset](https://www.10xgenomics.com/datasets/3-k-pbm-cells-from-a-healthy-donor-v-1-1-1-1-1) (2,700 peripheral blood mononuclear cells, 32,738 genes): **quality control → normalization → dimensionality reduction → clustering → cell-type annotation**.

This notebook is the narrated twin of the standalone script [`pbm3k_scanpy.py`](pbmc3k_scanpy.py) — every code cell below calls a function from that module, so the two can never drift apart. The same steps are also packaged as a **Nextflow pipeline** (`nextflow/main.nf`) with Docker and MultiQC reporting; see the repository `README.md`.
"""
)

md(
    """## Setup

```{note}
Run this notebook from its own directory (`analysis/`); the kernel working
directory is used to resolve relative paths. Outputs are written to
`../results/`, shared with the standalone script.
```
"""
)

code(
    """import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from IPython.display import Image, display

import pbmc3k_scanpy as pl   # the module behind analysis/pbmc3k_scanpy.py

sc.settings.verbosity = 3
sc.set_figure_params(dpi=120, frameon=False)

INPUT = "../data/pbmc3k_raw.h5ad"
OUTDIR = Path("../results")

print(f"scanpy {sc.__version__} | numpy {np.__version__} | pandas {pd.__version__}")
"""
)

md(
    """## 1. Load the raw expression matrix

The PBMC3k dataset is a raw 10x Genomics gene-barcode count matrix:
2,700 cells × 32,738 genes. `obtain_dataset.py` downloads it into `data/pbmc3k_raw.h5ad`.
"""
)

code(
    """adata = sc.read_h5ad(INPUT)
adata"""
)

md(
    """## 2. Quality control

Low-quality cells (few detected genes, extreme library sizes, high mitochondrial
read fraction — a hallmark of dying/lysed cells) and genes detected in almost no
cells are removed. Mitochondrial genes start with `MT-` (13 genes in this dataset).
"""
)

code(
    """# QC metrics + filtering (defaults: min 200 genes, max 2500 genes, <5% mitochondrial)
adata = pl.run_qc(adata, OUTDIR)

display(Image("../results/figures/qc_violins_raw.png"))
display(Image("../results/figures/qc_violins_filtered.png"))
pd.read_csv("../results/tables/qc_metrics.tsv", sep="\\t").T"""
)

md(
    """## 3. Normalization and highly variable genes

Library sizes are scaled to a common target (10,000 counts) and log-transformed
(`log1p`), which makes expression values comparable across cells. To reduce
noise and computation, only the most variable genes are kept for the embedding
(Seurat-style mean/dispersion cutoffs).
"""
)

code(
    """adata = pl.run_normalize(adata, OUTDIR)
display(Image("../results/figures/hvg_dispersion.png"))
print(f"HVGs kept: {adata.n_vars:,} of {adata.raw.n_vars:,} genes")"""
)

md(
    """## 4. Dimensionality reduction (PCA → kNN graph → UMAP)

The HVG matrix is projected onto its principal components (here 50), a
k-nearest-neighbour graph is built on the leading PCs, and UMAP embeds the graph
into two dimensions for visualization.
"""
)

code(
    """adata = pl.run_reduce(adata, OUTDIR)
display(Image("../results/figures/pca_variance.png"))
display(Image("../results/figures/umap_total_counts.png"))"""
)

md(
    """## 5. Clustering and marker genes

Leiden clustering on the kNN graph partitions the cells into discrete groups;
differential expression (Wilcoxon rank-sum) then identifies the genes that
define each cluster.
"""
)

code(
    """adata = pl.run_cluster(adata, OUTDIR)
display(Image("../results/figures/umap_leiden.png"))
display(Image("../results/figures/marker_heatmap_top5.png"))

markers = pd.read_csv("../results/tables/marker_genes.tsv", sep="\\t")
markers.head(10)"""
)

md(
    """## 6. Cell-type annotation

Each Leiden cluster is scored against curated marker panels for the major PBMC
populations (using mean log-normalized expression of the panel genes). The
best-scoring panel labels the cluster; the score gap to the runner-up is a
simple confidence measure.
"""
)

code(
    """adata = pl.run_annotate(adata, OUTDIR)
display(Image("../results/figures/umap_celltype.png"))
display(Image("../results/figures/dotplot_markers.png"))

composition = pd.read_csv("../results/tables/celltype_proportions.tsv", sep="\\t")
annotation = pd.read_csv("../results/tables/celltype_annotation.tsv", sep="\\t")
composition"""
)

code(
    """annotation"""
)

md(
    """## 7. Biological interpretation (summary)

After QC (2,700 → 2,638 cells, 2.3% removed), normalization and Leiden
clustering (resolution 0.8), the PBMC3k dataset resolves into **9 clusters that
collapse onto the 8 canonical peripheral-blood populations** of a healthy donor
(CD8 T cells are split into two clusters):

| Cell type | Cells | Fraction |
|---|---:|---:|
| CD4 T cells | 1,108 | 42.0% |
| CD14+ (classical) monocytes | 471 | 17.9% |
| CD8 T cells | 345 | 13.1% |
| B cells | 344 | 13.0% |
| FCGR3A+ (non-classical) monocytes | 169 | 6.4% |
| NK cells | 154 | 5.8% |
| Dendritic cells | 35 | 1.3% |
| Platelets (contaminating) | 12 | 0.5% |

**Key biological points**

1. **Composition matches a healthy donor.** T cells dominate (~55% CD4+CD8
   combined), followed by monocytes (~24% when both subsets are pooled), then B
   cells and NK cells; dendritic cells are rare (~1–2%). These proportions are
   the expected PBMC profile, which is reassuring for a QC sanity check.
2. **Marker-gene validation.** Every annotated population expresses its
   canonical markers: `CD3D`/`CD3E` (T cells), `CD8A`/`CD8B` (cytotoxic T),
   `NKG7`/`GNLY`/`KLRD1` (NK), `MS4A1`/`CD79A`/`CD79B` (B), `CD14`/`LYZ`/`S100A8`
   (classical monocytes), `FCGR3A`/`MS4A7` (non-classical monocytes),
   `FCER1A`/`CLEC10A` (dendritic cells), `PPBP`/`PF4` (platelets).
3. **QC decisions matter.** The mitochondrial-fraction (<5%) and gene-count
   filters removed 62 low-quality cells; the 12-cell platelet cluster is the
   well-known technical contamination of droplet scRNA-seq (platelet mRNA in
   the plasma fraction), not a true cell population.
4. **Clustering granularity is a trade-off.** At resolution 0.5 the rare
   dendritic cells merge into the monocyte cluster (they share the myeloid
   program, e.g. `LYZ`/`CST3`); resolution 0.8 separates them at the cost of a
   small extra CD8 T subdivision. This illustrates why cluster resolution must
   be tuned to the biological question.

**Caveats.** Annotation here is curated marker-panel scoring — fast and
transparent, but a proxy. A production analysis would add doublet detection
(e.g. Scrublet), reference-based label transfer (e.g. scArches/scANVI), and
finer subtype calls (e.g. naive vs. memory CD4 T cells via `CCR7`/`SELL` vs.
`S100A4`).
"""
)

md(
    """## Reproducibility

* **Script:** `analysis/pbmc3k_scanpy.py` — one entry point, five independent steps (`--step qc|normalize|reduce|cluster|annotate`).
* **Pipeline:** `nextflow/main.nf` orchestrates the steps as processes (Docker + MultiQC).
* **Environment:** pinned versions in `docker/requirements.txt` / `environment.yml`.
* **Session info:** exact package versions used to produce this notebook:
"""
)

code(
    """import warnings
warnings.filterwarnings("ignore")   # silence cosmetic FutureWarnings
import session_info
session_info.show()"""
)

nb["cells"] = cells

# --------------------------------------------------------------------------- #
#  Execute
# --------------------------------------------------------------------------- #

from nbclient import NotebookClient  # noqa: E402

client = NotebookClient(
    nb,
    timeout=900,
    kernel_name="python3",
    resources={"metadata": {"path": str(ANALYSIS)}},  # cwd for the kernel
)
print("Executing notebook ...", flush=True)
try:
    client.execute()
except Exception:
    print("Notebook execution FAILED", file=sys.stderr)
    raise

out_path = ANALYSIS / "pbmc3k_scanpy.ipynb"
nbf.write(nb, out_path)
print(f"Wrote executed notebook: {out_path}")

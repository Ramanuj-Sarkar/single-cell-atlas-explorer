"""Download the 10x Genomics "3k PBMCs from a Healthy Donor" dataset as an h5ad file.

The dataset is fetched through Scanpy (cached by the 10x Genomics public
bucket) and written to disk for downstream analysis. This script is the data
ingestion step of the Nextflow pipeline, and can also be run standalone.

Usage
-----
    python obtain_dataset.py [--output data/pbmc3k_raw.h5ad]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import scanpy as sc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data") / "pbmc3k_raw.h5ad",
        help="Where to write the raw AnnData (.h5ad) file.",
    )
    args = parser.parse_args()

    print("Downloading PBMC3k (3k PBMCs from a Healthy Donor) ...")
    adata = sc.datasets.pbmc3k()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    adata.write(args.output)
    print(f"Saved raw matrix {adata.n_obs:,} cells x {adata.n_vars:,} genes -> {args.output}")


if __name__ == "__main__":
    main()

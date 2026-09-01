#!/usr/bin/env nextflow
/*
 * PBMC3k Scanpy pipeline
 * ----------------------
 * Single-cell RNA-seq analysis of the 10x Genomics "3k PBMCs from a Healthy
 * Donor" dataset: quality control -> normalization -> dimensionality reduction
 * -> clustering -> cell-type annotation, with a MultiQC report.
 *
 * Every process wraps a step of `analysis/pbmc3k_scanpy.py`, which also exists
 * as an executed notebook (`analysis/pbmc3k_scanpy.ipynb`).
 *
 * Usage
 * -----
 *   nextflow run nextflow/main.nf -profile standard          # local execution
 *   nextflow run nextflow/main.nf -profile docker            # containerized
 *   nextflow run nextflow/main.nf -profile docker -resume    # cached rerun
 *
 * Typical parameters
 * ------------------
 *   --input    raw h5ad file (default: data/pbmc3k_raw.h5ad)
 *   --outdir   results directory (default: results_nextflow)
 *   --resolution  Leiden clustering resolution (default: 0.8)
 */

nextflow.enable.dsl = 2

include { OBTAIN }    from './modules/obtain'
include { QC }        from './modules/qc'
include { NORMALIZE } from './modules/normalize'
include { REDUCE }    from './modules/reduce'
include { CLUSTER }   from './modules/cluster'
include { ANNOTATE }  from './modules/annotate'
include { MULTIQC }   from './modules/multiqc'

workflow {

    // 1. raw matrix (download if missing)
    OBTAIN( Channel.fromPath(params.input) )

    // 2-6. analysis steps, chained through h5ad files
    QC( OBTAIN.out.raw )
    NORMALIZE( QC.out.h5ad )
    REDUCE( NORMALIZE.out.h5ad )
    CLUSTER( REDUCE.out.h5ad, params.resolution )
    ANNOTATE( CLUSTER.out.h5ad )

    // 7. MultiQC report aggregating the step-wise custom content
    MULTIQC(
        Channel.empty()
            .mix( QC.out.mqc, CLUSTER.out.mqc, ANNOTATE.out.mqc )
            .collect(),
        Channel.fromPath(params.multiqc_config)
    )

    // report the pipeline results
    ANNOTATE.out.h5ad    | view { "Annotated object:  $it" }
    MULTIQC.out.report   | view { "MultiQC report:    $it" }
}

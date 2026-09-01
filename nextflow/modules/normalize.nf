#!/usr/bin/env nextflow
/*
 * Step 2 - Normalization (library-size scaling + log1p) and HVG selection.
 */

process NORMALIZE {

    tag "NORMALIZE ${h5ad_in.name}"
    label 'cpu_1'
    publishDir "${params.outdir}/normalize", mode: 'copy'

    input:
    path h5ad_in

    output:
    path 'objects/normalize.h5ad',     emit: h5ad
    path 'tables/hvg_genes.tsv',       emit: hvg
    path 'figures/hvg_dispersion.png', emit: hvg_plot

    script:
    """
    ${params.python} ${params.script} \
        --input ${h5ad_in} \
        --outdir . \
        --step normalize
    """
}

#!/usr/bin/env nextflow
/*
 * Step 3 - Dimensionality reduction: PCA -> kNN graph -> UMAP.
 */

process REDUCE {

    tag "REDUCE ${h5ad_in.name}"
    label 'cpu_1'
    publishDir "${params.outdir}/reduce", mode: 'copy'

    input:
    path h5ad_in

    output:
    path 'objects/reduce.h5ad',                emit: h5ad
    path 'figures/pca_variance.png',           emit: pca_plot
    path 'figures/umap_total_counts.png',      emit: umap_plot

    script:
    """
    ${params.python} ${params.script} \
        --input ${h5ad_in} \
        --outdir . \
        --step reduce
    """
}

#!/usr/bin/env nextflow
/*
 * Step 4 - Leiden clustering + per-cluster marker genes (Wilcoxon).
 */

process CLUSTER {

    tag "CLUSTER ${h5ad_in.name}"
    label 'cpu_1'
    publishDir "${params.outdir}/cluster", mode: 'copy'

    input:
    path h5ad_in
    val resolution

    output:
    path 'objects/cluster.h5ad',                     emit: h5ad
    path 'tables/marker_genes.tsv',                  emit: markers
    path 'tables/multiqc_scanpy_markers_mqc.json',   emit: mqc
    path 'figures/umap_leiden.png',                  emit: umap_plot

    script:
    """
    ${params.python} ${params.script} \
        --input ${h5ad_in} \
        --outdir . \
        --step cluster \
        --resolution ${resolution}
    """
}

#!/usr/bin/env nextflow
/*
 * Step 5 - Cell-type annotation from curated marker panels + composition.
 */

process ANNOTATE {

    tag "ANNOTATE ${h5ad_in.name}"
    label 'cpu_1'
    publishDir "${params.outdir}/annotate", mode: 'copy'

    input:
    path h5ad_in

    output:
    path 'objects/annotate.h5ad',                        emit: h5ad
    path 'tables/celltype_proportions.tsv',              emit: proportions
    path 'tables/celltype_annotation.tsv',               emit: annotation
    path 'tables/multiqc_scanpy_celltypes_mqc.json',     emit: mqc
    path 'figures/umap_celltype.png',                    emit: umap_plot

    script:
    """
    ${params.python} ${params.script} \
        --input ${h5ad_in} \
        --outdir . \
        --step annotate
    """
}

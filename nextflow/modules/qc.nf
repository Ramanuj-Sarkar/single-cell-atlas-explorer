#!/usr/bin/env nextflow
/*
 * Step 1 - Quality control: cell/gene filtering + mitochondrial read fraction.
 * Writes the filtered h5ad plus QC tables, figures and MultiQC custom content.
 */

process QC {

    tag "QC ${h5ad_in.name}"
    label 'cpu_1'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path h5ad_in

    output:
    path 'objects/qc.h5ad',                          emit: h5ad
    path 'tables/multiqc_scanpy_qc_mqc.json',        emit: mqc
    path 'figures/*',                                emit: figures

    script:
    """
    ${params.python} ${params.script} \
        --input ${h5ad_in} \
        --outdir . \
        --step qc
    """
}

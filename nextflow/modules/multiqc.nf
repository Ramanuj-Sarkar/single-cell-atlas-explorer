#!/usr/bin/env nextflow
/*
 * Aggregate all step-wise MultiQC custom content into a single HTML report.
 */

process MULTIQC {

    tag "MultiQC"
    label 'cpu_1'
    publishDir "${params.outdir}/multiqc", mode: 'copy', pattern: 'multiqc_report.html'

    input:
    path mqc_files
    path multiqc_config

    output:
    path 'multiqc_report.html', emit: report

    script:
    """
    ${params.multiqc} ${mqc_files} \\
        --config ${multiqc_config} \\
        --filename multiqc_report \\
        --outdir .
    """
}

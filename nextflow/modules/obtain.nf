#!/usr/bin/env nextflow
/*
 * Obtain the raw expression matrix.
 * If the input file already exists it is passed through; otherwise the
 * dataset is downloaded via obtain_dataset.py.
 */

process OBTAIN {

    tag "OBTAIN ${input_path.name}"
    label 'cpu_1'
    publishDir "${params.outdir}/obtain", mode: 'copy', pattern: 'raw.h5ad'

    input:
    path input_path

    output:
    path 'raw.h5ad', emit: raw

    script:
    """
    if [ -f "${input_path}" ]; then
        cp "${input_path}" raw.h5ad
    else
        ${params.python} ${params.obtain_script} --output raw.h5ad
    fi
    """
}

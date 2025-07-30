nextflow.enable.dsl=2

process FASTQC {
    tag { "${family_id}" }
    label 'process_low'
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    //tuple val(sample_id), path(reads)
    tuple val(family_id), val(members), path(reads)

    output:
    path "*.{html,zip}", emit: fastqc_results
    path "*.{html,zip}", optional: true, emit: fastqc_results_raw

    script:
    """
    fastqc --quiet --threads ${task.cpus} ${reads}
    """
}

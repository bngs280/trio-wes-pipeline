process TRIOCNV_ANNOTATION {
    tag "${family}"
    publishDir "${params.outdir}/CNV_annotation/${family}", mode: 'copy'

    input:
    tuple val(family), val(sample_id), path(input_vcf)  // Takes the VCF from TRIOCNV
    path(cnv_prioritise_py)

    output:
    tuple val(family), val(sample_id), path("${sample_id}_cnv_annotated.tsv"), path("${sample_id}_CNV_annotated_ranked.xlsx"), emit: annotated_cnv

    script:
    """
    /usr/src/app/AnnotSV/bin/./AnnotSV -SVinputFile ${input_vcf} -outputFile ${sample_id}_cnv_annotated.tsv
    cp *AnnotSV/${sample_id}_cnv_annotated.tsv ${sample_id}_cnv_annotated.tsv
    /usr/bin/python3.10 ${cnv_prioritise_py} --annotSR ${sample_id}_cnv_annotated.tsv --output ${sample_id}_CNV_annotated_ranked.xlsx
    """
}

//path("${sample_id}_CNV_annotated_ranked.xlsx"),

process SV_ANNOTATION {
    tag "${family}"
    publishDir "${params.outdir}/SV_annotation/${family}", mode: 'copy'
    
    input:
    tuple val(family), val(sample_id), path(input_vcf)
    path(sv_prioritise_py)

    output:
    tuple val(family), val(sample_id), path("${sample_id}_sv_annotated.tsv"), path("${sample_id}_SV_annotated_ranked.xlsx"), emit: annotated_sv
    
    script:
    """
    # Run AnnotSV on SV VCF
    /usr/src/app/AnnotSV/bin/AnnotSV -SVinputFile ${input_vcf} -outputFile ${sample_id}_sv_annotated.tsv
        
    # Alternative if AnnotSV creates dated folders:
    cp *AnnotSV/${sample_id}_sv_annotated.tsv ${sample_id}_sv_annotated.tsv

    # Prioritization
    /usr/bin/python3 ${sv_prioritise_py} --annotSR ${sample_id}_sv_annotated.tsv --output ${sample_id}_SV_annotated_ranked.xlsx
    """
}

//path("${sample_id}_SV_annotated_ranked.xlsx"),

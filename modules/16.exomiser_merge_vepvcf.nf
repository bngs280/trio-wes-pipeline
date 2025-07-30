process mergeExomiserVEPResults {
    tag "Final P2G merge for ${family_id}"
    publishDir "${params.outdir}/final_results", mode: 'copy'

    conda 'pandas numpy'

    input:
    tuple val(family_id), path(variants_tsv), path(variants_json), path(genes_tsv)  // From run_prioritizer (exact match)
    tuple val(family_id), path(vep_annotated_vcf)  // From trioG2Pprocessing
    path(p2g_exomisermerging_script)   

    output:
    tuple val(family_id), path("${family_id}_P2G.csv"), emit: final_p2g_results
    
    script:
    """
    # Create required directories
    mkdir -p annotation_dir
    mkdir -p patient_dir

    # Debug: Show input files
    echo "Using VEP file: ${vep_annotated_vcf}"
    echo "Using variants TSV: ${variants_tsv}"
    echo "Using variants JSON: ${variants_json}"
    echo "Using genes TSV: ${genes_tsv}"

    # Run the postprocessing script
    /usr/bin/python3 ${p2g_exomisermerging_script} \
        --vep        ${vep_annotated_vcf} \
        --json       ${variants_json} \
        --variants   ${variants_tsv} \
        --output     ${family_id}

    # Verify outputs were created
    if [ ! -f "${family_id}_P2G.csv" ]; then
        echo "ERROR: Failed to generate CSV output" >&2
        exit 1
    fi
    """
}

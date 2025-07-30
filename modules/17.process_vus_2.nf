// Process VEP Output
process processVEPOutput {
    tag "Process VEP output for ${family_id}"
    publishDir "${params.outdir}/10.VUSprioritised", mode: 'copy'
    cpus 14

    input:
    tuple val(family_id), path(vus_vep)

    output:
    tuple val(family_id), path("${family_id}_processed.txt"), emit: processed_vcf

    script:
    """
    # Process VEP output in one step without temporary file
    grep -v '^##' ${vus_vep} | \
    sed '1s/Eigen-phred_coding/Eigen-pred_coding/' > ${family_id}_processed.txt
    """
}

// Run VUS Prioritization
process runVusPrioritization {
    tag "Run VUSprize for ${family_id}"
    publishDir "${params.outdir}/10.VUSprioritised", mode: 'copy'
    cpus 14

    input:
    tuple val(family_id), path(processed_vus)
    val(params.vusprize_script)
    path(params.dependent_script)

    output:
    tuple val(family_id), path("${family_id}_vusScore.txt"), emit: vusprize_output

    script:
    """
    # Set up environment
    export WORKDIR=\$(pwd)
    export INPUT_FILE=${processed_vus}
    
    cp ${params.vusprize_script} .
    cp ${params.dependent_script} .
 
    /opt/miniconda3/envs/vusprize/bin/python \
      VusPrize.py \
      "\$INPUT_FILE" \
      "${family_id}_vusScore.txt"
  
    """
}

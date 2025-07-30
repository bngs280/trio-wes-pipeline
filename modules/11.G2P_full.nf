process fullG2PPrioritisation {
    label 'python'
    label 'high_memory'
    publishDir "${params.outdir}/Final_Prioritisation", mode: 'copy'
    
    input:
    tuple val(family_id), path(annotated_tsv)
    val(params.omim_clingene_file)
    val(params.mondo_owl)
    val(params.child_sex)
    val(params.family_history)
    tuple val("g2p_script"), path(g2p_script)
    
    output:
    tuple val(family_id), path("${family_id}_firstlist_pvs1.csv"), emit: firstlist
    tuple val(family_id), path("${family_id}_top_PVS1.csv"), emit: toplist
    tuple val(family_id), path("${family_id}_missense_final.csv"), emit: missense
    
    script:
    """
    /usr/bin/python3.10 ${g2p_script} \\
        --vcfvep                 ${annotated_tsv} \\
        --omim_clingene          ${params.omim_clingene_file} \\
        --mondo_owl              ${params.mondo_owl} \\
        --child_sex              ${params.child_sex} \\
        --family_history         ${params.family_history} \\
        --firstlist_output_csv   ${family_id}_firstlist_pvs1.csv \\
        --output_csv             ${family_id}_top_PVS1.csv \\
        --missense_output_csv    ${family_id}_missense_final.csv
    """
}

//--child_sex ${child_sex} \\


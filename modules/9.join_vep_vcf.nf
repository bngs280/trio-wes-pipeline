process trioG2Pprocessing { 
    label 'python'
    label 'high_memory'  // Adjust based on your script's requirements
    publishDir "${params.outdir}/fastqc/VEPVCF_Results", mode:'copy'

    input:
    tuple val(family_id), path(vep_txt)
    tuple val(family_id), path(raw_vcf)
    tuple val(family_id), path(acmg)
    path(params.triosEXOped)
    tuple val("clingen"), path(clingen_file)
    tuple val("script"), path(vep_vcf_script)
    //val(params.ped_file)
    //val(params.clingen_file)
    //val(params.vep_vcf_script)
    
    output:
    tuple val(family_id), path("${family_id}_annotated_trio.tsv"), emit: vcfvep_output
    path("python_logs/*.log"), optional: true, emit: logs
    
    script:
    """
    mkdir -p python_logs
    
    # Run your specific trio processing script
    /usr/bin/python3 ${vep_vcf_script} \\
        --vep ${vep_txt} \\
        --clingen ${params.clingen_file} \\
        --ped ${params.triosEXOped} \\
        --vcf ${raw_vcf} \\
        --acmg ${acmg} \\
        --output ${family_id}_annotated_trio.tsv
    """
}

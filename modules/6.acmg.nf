process acmgpre {
    tag { "${family_id}" }
    publishDir "${params.outdir}/variant_calling/acmg", mode: 'copy'

    input:
    tuple val(family_id), path(vcf_file)
    path(params.cache_bias)
    val(params.sd_bias)
    val(params.ref_bias)

    output:
    tuple val(family_id), path("${family_id}.json.gz"), emit: acmg_pre

    script:
    """
    dotnet /usr/src/app/Nirvana-v3.18.1/Nirvana.dll \
        --cache ${params.cache_bias} \
        --sd ${params.sd_bias} \
        --ref ${params.ref_bias} \
        --in ${vcf_file} \
        --out ${family_id}
    """
}

process acmgclass {
    tag { "${family_id}" }
    publishDir "${params.outdir}/variant_calling/acmg", mode: 'copy'

    input:
    tuple val(family_id), path(json_file)
    path(params.parameter_bias)

    output:
    tuple val(family_id), path("${family_id}_ACMG.tsv"), emit: acmg_class

    script:
    """
    python3 /usr/src/app/BIAS-2015/bias_2015.py \
        ${json_file} \
        ${params.parameter_bias} \
        ${family_id}_ACMG.tsv
    """
}

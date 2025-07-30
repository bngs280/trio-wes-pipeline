process Denom {
    //container 'google/deepvariant:1.9.0'
    tag { "${family_id}" }
    publishDir "${params.outdir}/variant_calling", mode: 'copy'
    cpus 16

    input:
    tuple val(family_id), val(members), path(vcf_file)
    //path ped_file
    path(triosEXOped)


    output:

    tuple val(family_id), path("${family_id}_${members.proband}_cohort.dnm.vcf"), emit: denovo_vcf
    script:
    """
    bcftools +trio-dnm2 \
        -P ${triosEXOped} \
        --use-NAIVE \
        ${vcf_file} > ${family_id}_${members.proband}_cohort.dnm.vcf
    """
}

process phasingDnm {
    //container 'google/deepvariant:1.9.0'
    tag { "${family_id}" }
    publishDir "${params.outdir}/variant_calling", mode: 'copy'
    cpus 16

    input:
    tuple val(family_id), path(denovo_vcf)
    tuple val(family_id),
          val(members),
          path(proband_bam),
          path(proband_bai),
          path(father_bam),
          path(father_bai),
          path(mother_bam),
          path(mother_bai)
    path(params.ref)
    path(params.triosEXOped)

    output:

    tuple val(family_id), path("${family_id}_${members.proband}_cohort.dnm.phas.vcf"), emit: phased_vcf
    tuple val(family_id), path("${family_id}_${members.proband}_cohort.dnm.phas.vcf.gz"), emit: final_vcf_phased
    
    script:
    """
    /root/.local/bin/./whatshap phase \
        --reference ${params.ref} \
        --ped ${params.triosEXOped} \
        --output ${family_id}_${members.proband}_cohort.dnm.phas.vcf \
        ${denovo_vcf} \
        ${proband_bam} \
        ${father_bam} \
        ${mother_bam}
    
    # gzip
    bgzip -c ${family_id}_${members.proband}_cohort.dnm.phas.vcf > ${family_id}_${members.proband}_cohort.dnm.phas.vcf.gz

    # index
    bcftools index ${family_id}_${members.proband}_cohort.dnm.phas.vcf.gz
    """
}

//, path("${child_sample}_cohort_filtered_PASS.vcf")
// custom filter python script replacemtn of bcf tool PASS
// F M > yes Proband No
// GQ > 20
//DP > 10
// if GT ./. in all 


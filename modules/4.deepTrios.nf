process DEEPTRIO {
    //container 'google/deepvariant:1.9.0'
    tag { "${family_id}" }
    publishDir "${params.outdir}/variant_calling", mode: 'copy'
    cpus 16

    input:
    tuple val(family_id),
          val(members),
          path(proband_bam),
          path(proband_bai),
          path(father_bam),
          path(father_bai),
          path(mother_bam),
          path(mother_bai)
    path(params.ref)
    path(params.bed)

    output:
    tuple val(family_id),
          path("${family_id}_${members.proband}.vcf.gz"),
          path("${family_id}_${members.proband}.vcf.gz.tbi"),
          path("${family_id}_${members.father}.vcf.gz"),
          path("${family_id}_${members.father}.vcf.gz.tbi"),
          path("${family_id}_${members.mother}.vcf.gz"),
          path("${family_id}_${members.mother}.vcf.gz.tbi"),
          emit: vcfs
    tuple val(family_id),
          path("${family_id}_${members.proband}.g.vcf.gz"),
          path("${family_id}_${members.proband}.g.vcf.gz.tbi"),
          path("${family_id}_${members.father}.g.vcf.gz"),
          path("${family_id}_${members.father}.g.vcf.gz.tbi"),
          path("${family_id}_${members.mother}.g.vcf.gz"),
          path("${family_id}_${members.mother}.g.vcf.gz.tbi"),
          emit: gvcfs
    tuple val(family_id), path("${family_id}_${members.proband}_cohort.vcf.gz"), emit: cohortvcfs
    tuple val(family_id), path("${family_id}_${members.proband}_cohort_filtered_PASS.vcf.gz"), emit: final_cohort_vcfs

    script:
    """
    rm -rf make_examples* call_variants_output* gvcf*
    /opt/deepvariant/bin/deeptrio/run_deeptrio \
    --model_type WES \
    --ref ${params.ref} \
    --reads_child ${proband_bam} \
    --reads_parent1 ${father_bam} \
    --reads_parent2 ${mother_bam} \
    --output_vcf_child ${family_id}_${members.proband}.vcf.gz \
    --output_vcf_parent1 ${family_id}_${members.father}.vcf.gz \
    --output_vcf_parent2 ${family_id}_${members.mother}.vcf.gz \
    --output_gvcf_child ${family_id}_${members.proband}.g.vcf.gz \
    --output_gvcf_parent1 ${family_id}_${members.father}.g.vcf.gz \
    --output_gvcf_parent2 ${family_id}_${members.mother}.g.vcf.gz \
    --num_shards ${task.cpus} \
    --sample_name_child ${members.proband} \
    --sample_name_parent1 ${members.father} \
    --sample_name_parent2 ${members.mother} \
    ${params.bed} \
    --intermediate_results_dir deeptrio_${family_id}

      ## glnexus_cli
      /usr/local/bin/glnexus_cli --config DeepVariantWES \
      --threads ${task.cpus} \
      ${family_id}_${members.proband}.g.vcf.gz \
      ${family_id}_${members.father}.g.vcf.gz \
      ${family_id}_${members.mother}.g.vcf.gz \
      | bcftools view -Ou - \
      | bcftools view -Oz -o ${family_id}_${members.proband}_cohort.vcf.gz

      bcftools view -f PASS ${family_id}_${members.proband}_cohort.vcf.gz > ${family_id}_${members.proband}_cohort_filtered_PASS.vcf
      bgzip -c ${family_id}_${members.proband}_cohort_filtered_PASS.vcf > ${family_id}_${members.proband}_cohort_filtered_PASS.vcf.gz
      tabix -p vcf ${family_id}_${members.proband}_cohort_filtered_PASS.vcf.gz
    """
}
//, path("${child_sample}_cohort_filtered_PASS.vcf")
// custom filter python script replacemtn of bcf tool PASS
// F M > yes Proband No
// GQ > 20
//DP > 10
// if GT ./. in all 
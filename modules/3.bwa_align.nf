// process BWA_ALIGN {
//     tag { "${family_id}" }
//     label 'bwa_align'
//     publishDir "${params.outdir}/bwa", mode: 'copy'

//     input:
//     tuple val(family_id),
//           path(proband_reads),
//           path(father_reads),
//           path(mother_reads)

//     output:
//     tuple val(family_id), path("${family_id}_${sample}.bam"), emit: aligned_bams

//     script:
//     def bam_outputs = []
//     def samples = ["proband", "father", "mother"]
//     def all_reads = [proband_reads, father_reads, mother_reads]

//     """
//     set -e

//     bwa_index="${params.bwa_index}"

//     # Proband
//     bwa mem \$bwa_index ${proband_reads[0]} ${proband_reads[1]} | \
//         samtools sort -o ${family_id}_proband.bam

//     # Father
//     bwa mem \$bwa_index ${father_reads[0]} ${father_reads[1]} | \
//         samtools sort -o ${family_id}_father.bam

//     # Mother
//     bwa mem \$bwa_index ${mother_reads[0]} ${mother_reads[1]} | \
//         samtools sort -o ${family_id}_mother.bam
//     """
// }
// worked

process BWA {
    tag { "${family_id}" }
    label 'process_high'
    publishDir "${params.outdir}/alignment", mode: 'copy'
    cpus 16
    
    input:
    tuple val(family_id), 
          val(members),
          path(proband_R1),
          path(proband_R2),
          path(father_R1),
          path(father_R2),
          path(mother_R1),
          path(mother_R2)
    path(params.ref)

    output:
    tuple val(family_id),
          val(members),
          path("${family_id}_${members.proband}.sorted.md.bam"),
          path("${family_id}_${members.proband}.sorted.md.bam.bai"),
          path("${family_id}_${members.father}.sorted.md.bam"),
          path("${family_id}_${members.father}.sorted.md.bam.bai"),
          path("${family_id}_${members.mother}.sorted.md.bam"),
          path("${family_id}_${members.mother}.sorted.md.bam.bai"),
          emit: bams

    script:
    """
    # Align proband
    bwa mem -t ${task.cpus} -M \\
        -R "@RG\\tID:${family_id}_${members.proband}\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:${family_id}_${members.proband}" \\
        ${params.ref} ${proband_R1} ${proband_R2} > ${family_id}_${members.proband}.sam
    samtools view -@ ${task.cpus} -b ${family_id}_${members.proband}.sam | \\
        samtools sort -@ ${task.cpus} -o ${family_id}_${members.proband}.sorted.bam -
    sambamba markdup -r -t ${task.cpus} ${family_id}_${members.proband}.sorted.bam ${family_id}_${members.proband}.sorted.md.bam 
    samtools index ${family_id}_${members.proband}.sorted.md.bam 
    rm ${family_id}_${members.proband}.sam
    
    # Align father
    bwa mem -t ${task.cpus} -M \\
        -R "@RG\\tID:${family_id}_${members.father}\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:${family_id}_${members.father}" \\
        ${params.ref} ${father_R1} ${father_R2} > ${family_id}_${members.father}.sam
    samtools view -@ ${task.cpus} -b ${family_id}_${members.father}.sam | \\
        samtools sort -@ ${task.cpus} -o ${family_id}_${members.father}.sorted.bam -
    sambamba markdup -r -t ${task.cpus} ${family_id}_${members.father}.sorted.bam ${family_id}_${members.father}.sorted.md.bam
    samtools index ${family_id}_${members.father}.sorted.md.bam
    rm ${family_id}_${members.father}.sam
    
    # Align mother
    bwa mem -t ${task.cpus} -M \\
        -R "@RG\\tID:${family_id}_${members.mother}\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:${family_id}_${members.mother}" \\
        ${params.ref} ${mother_R1} ${mother_R2} > ${family_id}_${members.mother}.sam
    samtools view -@ ${task.cpus} -b ${family_id}_${members.mother}.sam | \\
        samtools sort -@ ${task.cpus} -o ${family_id}_${members.mother}.sorted.bam -
    sambamba markdup -r -t ${task.cpus} ${family_id}_${members.mother}.sorted.bam ${family_id}_${members.mother}.sorted.md.bam
    samtools index ${family_id}_${members.mother}.sorted.md.bam
    rm ${family_id}_${members.mother}.sam
    """
}
// process PROCESS_FASTQ {
//     tag { "${sample_id}" }
//     label 'process_medium'
//     publishDir "${params.outdir}/fastp", mode: 'copy'
    
//     input:
//     tuple val(sample_id), path(reads)
    
//     output:
//     tuple val(sample_id), path("*.fastq.gz"), emit: processed_reads
//     path "*.log", emit: logs
    
//     script:
//     def read1 = reads[0]
//     def read2 = reads.size() > 1 ? reads[1] : ''
//     def read3 = reads.size() > 2 ? reads[2] : ''
//     def read4 = reads.size() > 3 ? reads[3] : ''
//     def read5 = reads.size() > 4 ? reads[4] : ''
//     def read6 = reads.size() > 5 ? reads[5] : ''
    
//     def output_prefix = "trimmed_${sample_id}"
    
//     """
//     # For trio WES (3 pairs)
//     cutadapt \
//         -p ${output_prefix}_R2.fastq.gz \
//         -o ${output_prefix}_R1.fastq.gz \
//         ${read1} ${read2} \
//         > ${output_prefix}.cutadapt.log 2>&1
        
//     cutadapt \
//         -p ${output_prefix}_R4.fastq.gz \
//         -o ${output_prefix}_R3.fastq.gz \
//         ${read3} ${read4} \
//         >> ${output_prefix}.cutadapt.log 2>&1
        
//     cutadapt \
//         -p ${output_prefix}_R6.fastq.gz \
//         -o ${output_prefix}_R5.fastq.gz \
//         ${read5} ${read6} \
//         >> ${output_prefix}.cutadapt.log 2>&1
//     """
// }

// In fastq_process.nf
// worekd
// process PROCESS_FASTQ {
//     tag { "${family_id}" }  // Changed from sample_id to family_id
//     label 'process_medium'
//     publishDir "${params.outdir}/fastp", mode: 'copy'
    
//     input:
//     //tuple val(family_id), val(members), path(reads)  // Now receives members
//     tuple val(family_id), val(members), path(reads)
    
//     output:
//     tuple val(family_id), 
//           path("trimmed_${family_id}_${members.proband}_R{1,2}.fastq.gz"),
//           path("trimmed_${family_id}_${members.father}_R{1,2}.fastq.gz"),
//           path("trimmed_${family_id}_${members.mother}_R{1,2}.fastq.gz"),
//           emit: processed_reads
//     path "*.log", emit: logs
    
//     script:
//     """
//     # Process proband (original R1/R2)
//     cutadapt \
//         -o trimmed_${family_id}_${members.proband}_R1.fastq.gz \
//         -p trimmed_${family_id}_${members.proband}_R2.fastq.gz \
//         ${reads[0]} ${reads[1]} \
//         > ${family_id}.cutadapt.log 2>&1
        
//     # Process father
//     cutadapt \
//         -o trimmed_${family_id}_${members.father}_R1.fastq.gz \
//         -p trimmed_${family_id}_${members.father}_R2.fastq.gz \
//         ${reads[2]} ${reads[3]} \
//         >> ${family_id}.cutadapt.log 2>&1
        
//     # Process mother
//     cutadapt \
//         -o trimmed_${family_id}_${members.mother}_R1.fastq.gz \
//         -p trimmed_${family_id}_${members.mother}_R2.fastq.gz \
//         ${reads[4]} ${reads[5]} \
//         >> ${family_id}.cutadapt.log 2>&1
//     """
// }

process PROCESS_FASTQ {
    tag { "${family_id}" }  // Changed from sample_id to family_id
    label 'process_medium'
    publishDir "${params.outdir}/fastp", mode: 'copy'
    
    input:
    //tuple val(family_id), val(members), path(reads)  // Now receives members
    tuple val(family_id), val(members), path(reads)
    
    output:
    tuple val(family_id),
        path("trimmed_${family_id}_${members.proband}_R1.fastq.gz"),
        path("trimmed_${family_id}_${members.proband}_R2.fastq.gz"),
        path("trimmed_${family_id}_${members.father}_R1.fastq.gz"),
        path("trimmed_${family_id}_${members.father}_R2.fastq.gz"),
        path("trimmed_${family_id}_${members.mother}_R1.fastq.gz"),
        path("trimmed_${family_id}_${members.mother}_R2.fastq.gz"),
        emit: processed_reads
    path "*.log", emit: logs
    // output:
    // tuple val(family_id), 
    //       path("trimmed_${family_id}_${members.proband}_R{1,2}.fastq.gz"),
    //       path("trimmed_${family_id}_${members.father}_R{1,2}.fastq.gz"),
    //       path("trimmed_${family_id}_${members.mother}_R{1,2}.fastq.gz"),
    //       emit: processed_reads
    // path "*.log", emit: logs
    
    script:
    """
    # Process proband (original R1/R2)
    cutadapt \
        -o trimmed_${family_id}_${members.proband}_R1.fastq.gz \
        -p trimmed_${family_id}_${members.proband}_R2.fastq.gz \
        ${reads[0]} ${reads[1]} \
        > ${family_id}.cutadapt.log 2>&1
        
    # Process father
    cutadapt \
        -o trimmed_${family_id}_${members.father}_R1.fastq.gz \
        -p trimmed_${family_id}_${members.father}_R2.fastq.gz \
        ${reads[2]} ${reads[3]} \
        >> ${family_id}.cutadapt.log 2>&1
        
    # Process mother
    cutadapt \
        -o trimmed_${family_id}_${members.mother}_R1.fastq.gz \
        -p trimmed_${family_id}_${members.mother}_R2.fastq.gz \
        ${reads[4]} ${reads[5]} \
        >> ${family_id}.cutadapt.log 2>&1
    """
}
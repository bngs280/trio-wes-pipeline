process SV_CALLING {
    tag "${family}"
    publishDir "${params.outdir}/SV_calling/${family}", mode: 'copy'
    
    input:
    tuple val(family), val(members), val(proband_bam), val(father_bam), val(mother_bam)
    path(params.ref)
    
    output:
    tuple val(family), val(sample_id), path("${family}_manta/results/variants/diploidSV.vcf.gz"), emit: manta_vcf
    
    script:
    sample_id = proband_bam.getBaseName().replace('.sorted.md.bam', '')
    """
    # Increase file handle limit
    ulimit -n 4096
    
    # Activate conda environment
    source /opt/miniconda3/etc/profile.d/conda.sh
    conda activate p2Manta
    
    # 1. Configure Manta
    configManta.py \\
        --bam ${proband_bam} \\
        --bam ${father_bam} \\
        --bam ${mother_bam} \\
        --referenceFasta ${params.ref} \\
        --runDir ${family}_manta
    
    # 2. Run Manta
    ${family}_manta/runWorkflow.py -j 8
    
    # 3. Verify output
    if [ ! -f "${family}_manta/results/variants/diploidSV.vcf.gz" ]; then
        echo "ERROR: Manta failed to generate output VCF"
        exit 1
    fi
    """
}

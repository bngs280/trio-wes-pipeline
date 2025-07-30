process TRIOCNV {
    tag "${family}"
    publishDir "${params.outdir}/fastqc/CNV_calling/${family}", mode: 'copy'

    input:
    tuple val(family), val(members), path(proband_bam), path(father_bam), path(mother_bam)
    path(ref)
    path(mappability)
    path(triosEXOped)
    path(cnv2vcf_py)

    output:
    // path * "${family}_output"
    tuple val(family), val(sample_id), path("${sample_id}_cnv.vcf"), emit: cnv_vcf

    script:
    sample_id = proband_bam.getBaseName().replace('.sorted.md.bam', '')
    """
    # Remove the header line from the PED file (creates a temp file)
    grep -v '^#' ${triosEXOped} > ped_noheader.ped
  
    mkdir -p "${family}_output"
    # Get absolute paths using readlink -f (works in bash)
    # This resolves the symlinks to their actual absolute paths
    PROBAND_ABS=\$(readlink -f ${proband_bam})
    FATHER_ABS=\$(readlink -f ${father_bam})
    MOTHER_ABS=\$(readlink -f ${mother_bam})
    
    echo "=== Resolved absolute paths ==="
    echo "Proband: \$PROBAND_ABS"
    echo "Father: \$FATHER_ABS"
    echo "Mother: \$MOTHER_ABS"

    # Create bamList.txt with absolute paths
    # TrioCNV2 expects absolute paths, not staged filenames
    printf "\$PROBAND_ABS\\n" > bamList.txt
    printf "\$FATHER_ABS\\n" >> bamList.txt
    printf "\$MOTHER_ABS\\n" >> bamList.txt

    # Debug output
    echo "=== bamList.txt contents (with absolute paths) ==="
    cat bamList.txt
    
    # Extract sample IDs from BAM files (remove .bam extension)
    PROBAND_ID=\$(basename "${proband_bam}" .sorted.md.bam)
    FATHER_ID=\$(basename "${father_bam}" .sorted.md.bam)  
    MOTHER_ID=\$(basename "${mother_bam}" .sorted.md.bam)

    # Update the PED file in one step
    awk -v pid="\$PROBAND_ID" -v fid="\$FATHER_ID" -v mid="\$MOTHER_ID" '
    BEGIN {OFS="\t"} 
    {
        if (NR == 1) {\$2 = pid; \$3 = fid; \$4 = mid}
        if (NR == 2) {\$2 = fid}
        if (NR == 3) {\$2 = mid}
        print
    }' ped_noheader.ped > updated_ped.ped

    # Verify result
    echo "Updated PED file:"
    cat updated_ped.ped

    # Verify the absolute path files exist
    echo "=== Verifying absolute path BAM files ==="
    while IFS=\$'\\t' read -r bam_file sample_id; do
        echo "Checking: '\$bam_file' for sample '\$sample_id'"
        if [[ -f "\$bam_file" ]]; then
            echo "  ✓ Found: \$bam_file"
            ls -la "\$bam_file"
        else
            echo "  ✗ NOT FOUND: \$bam_file"
            exit 1
        fi
    done < bamList.txt


    echo "=== Running TrioCNV2 with absolute paths ==="
    java -jar /usr/local/bin/TrioCNV2-0.1.1.jar preprocess \\
        -R ${ref} \\
        -B bamList.txt \\
        -P updated_ped.ped \\
        -M ${mappability} \\
        -O ${family}_output/

    java -jar /usr/local/bin/TrioCNV2-0.1.1.jar call -I ${family}_output -P updated_ped.ped -M ${mappability} -O ${family}_output --nt 100
    /usr/bin/python3 ${cnv2vcf_py} ${family}_output/DiscordantReadPairs.txt ${sample_id}_cnv.vcf
    """
}

//python /home/vgenomics/Germline/tools/TrioCNV2/TrioCNV2/resuCALL300625/Prioritization/TrioCNV2_vcf.py ${family}_output/call/DiscordantReadPairs.txt ${family}_cnv.vcf

process MODIFY_CONFIG {
    tag { "${family_id}" }
    publishDir "${params.outdir}/exomiser", mode: 'copy'

    input:
    path template_yml
    path ped_file
    val phenotypes
    val genome_version

    output:
    path 'modified-config.yml', emit: config

    script:
    def hpoString = phenotypes.collect { "'${it}'" }.join(", ")
    """
    #!/bin/bash
    set -euo pipefail  # Strict error handling
    
    # Create output directory
    mkdir -p "${params.outdir}/exomiser"
    
    # Extract proband ID from PED file (first individual with affected status = 2)
    PROBAND_ID=\$(awk '\$6 == "2" {print \$2; exit}' "${ped_file}")
    if [[ -z "\$PROBAND_ID" ]]; then
        echo "ERROR: No affected individual found in PED file"
        exit 1
    fi
    
    # Process the template while maintaining exact structure
    {
        # Keep header until analysis section
        grep -B 1000 '^analysis:' "${template_yml}" | head -n -1
        
        # Print analysis section with modifications
        echo "analysis:"
        grep -A 1000 '^analysis:' "${template_yml}" | \
        head -n 1 | \
        sed "s/genomeAssembly:.*/genomeAssembly: ${genome_version}/"
        
        # Print vcf and ped lines unchanged
        grep -A 1000 '^analysis:' "${template_yml}" | \
        grep -e '^    vcf:' -e '^    ped:'
        
        # Print proband line with extracted ID
        grep -A 1000 '^analysis:' "${template_yml}" | \
        grep '^    proband:' | \
        sed "s/proband:.*/proband: \$PROBAND_ID/"
        
        # Print hpoIds with new values
        grep -A 1000 '^analysis:' "${template_yml}" | \
        grep '^    hpoIds:' | \
        sed "s/hpoIds:.*/hpoIds: [${hpoString}]/"
        
        # Print everything else unchanged
        grep -A 1000 '^analysis:' "${template_yml}" | \
        tail -n +6 | \
        grep -v '^    proband:' | \
        grep -v '^    hpoIds:'
    } > modified-config.yml
    
    # Verification
    echo "=== Modified Config ==="
    echo "Proband ID set to: \$PROBAND_ID"
    grep -A 3 "proband:" modified-config.yml
    echo "=== First 15 lines ==="
    head -n 15 modified-config.yml
    """
}
// Process definition
process run_prioritizer {
    tag { "${family_id}" }
    publishDir "${params.outdir}/exomiser", mode: 'copy'

    input:
    tuple val(family_id), path(vcf_file)
    path config_yml
    val(genome)
    val(exomiser_properties)
    val(triosEXOped)
    path(exomiser_jar)

    output:
    path "${family_id}/*"
    tuple val(family_id), path("${family_id}/*variants.tsv"), path("${family_id}/*.json"), path("${family_id}/*.genes.tsv"), emit: variants_tsv_json
    //tuple val(sample_id), path("${sample_id}/*.json"), emit: variants_json

    script:
    """
    mkdir -p "${family_id}"
    java -jar ${exomiser_jar} \
        --analysis ${config_yml} \
        --output-directory ${family_id} \
        --assembly ${genome} \
        --vcf ${vcf_file} \
        --output-filename ${family_id} \
        --ped ${triosEXOped} \
        --spring.config.location=${exomiser_properties}
    """
}

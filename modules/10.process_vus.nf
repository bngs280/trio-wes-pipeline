process vusVEP {
    tag "VEP on ${family_id}"
    publishDir "${params.outdir}/9.vus_vep", mode: 'copy'
    cpus 15

    input:
    tuple val(family_id), path(vcf_file)
    val(params.ref)
    val(params.cache_dir)
    val(params.dirPlugin)
    val(params.dbnsfp4)
    val(params.LoFtool)
    val(params.CADDsnv)
    val(params.CADDindel)
    val(params.dbscSNV)
    val(params.assembly)
    path(params.maxentscan)

    output:
    tuple val(family_id), path("${family_id}_filtered_VUS_VEP.txt"), emit: vus_vep

    script:
    """
    /usr/src/app/ensembl-vep/./vep \
    --af --appris --buffer_size 50000 --cache --check_existing --distance 5000 \
    --assembly ${params.assembly} \
    --dir ${params.cache_dir} \
    --dir_plugins ${params.dirPlugin} \
    --mane --pick \
    --fasta ${params.ref} \
    --force --fork 15 \
    --input_file ${vcf_file} --format vcf \
    --output_file ${family_id}_filtered_VUS_VEP.txt --tab \
    --plugin MaxEntScan,${params.maxentscan} \
    --plugin Blosum62 \
    --plugin dbNSFP,${params.dbnsfp4},codon_degeneracy,Eigen-phred_coding,integrated_fitCons_score,GERP++_RS,phyloP100way_vertebrate \
    --plugin dbscSNV,${params.dbscSNV} \
    --plugin LoFtool,${params.LoFtool} \
    --plugin CADD,snv=${params.CADDsnv},indel=${params.CADDindel} \
    --refseq --quiet --regulatory --sift s --species homo_sapiens --symbol --transcript_version --tsl --safe --show_ref_allele --stats_text --uploaded_allele
    """
}

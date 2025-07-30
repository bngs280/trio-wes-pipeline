#!/usr/bin/env nextflow
// nextflow.enable.dsl=2

// Version 1
// // Genome reference selection
// if (params.genome == 'hg38') {
//     params.ref = params.refhg38
//     params.bed = params.bedhg38
//     params.assembly = 'GRCh38'
//     //params.cachee = params.cache_dir
//     params.dirPlugin = params.dirplugin
//     params.dbNSFP = params.dbNSFP38
//     params.CADDsnv = params.CADDsnv38
//     params.CADDindel = params.CADDindels38
//     params.dbscSNV = params.dbscSNV38
//     params.exomiser_properties = params.exomiser_propertieshg38
//     //params.vusprize_script = params.vusprize_script
//     //params.dependent_script = params.dependent_script
//     params.spliceaiindel = params.spliceaiindel38
//     params.spliceaisnv = params.spliceaisnv38
//     params.phenotypesfile = params.phenotypesfile38
//     params.loeuf = params.loeuf38
//     params.clinpred = params.clinpred38
//     params.alphamissense = params.alphamissense38
//     params.primateai = params.primateai38
//     params.revel = params.revel38
//     params.clinvar = params.clinvar38
//     params.clingen_file = params.clingen_filehg38
//     params.pli = params.plihg38
//     params.mappability = params.mappabilityhg38
//     params.dbnsfp4 = params.dbnsfp4_38
//     // ACMG
//     params.cache_bias = params.cache_biashg38
//     params.sd_bias = params.sd_biashg38
//     params.ref_bias = params.ref_biashg38
//     params.parameter_bias = params.parameter_biashg38

// } else if (params.genome == 'hg19') {
//     params.ref = params.refhg37
//     params.ref_dir = params.ref_dirhg37
//     params.reffai = params.refhg37fai
//     params.bed = params.bedhg37
//     params.gnomad = params.gnomad37
//     params.db1000g = params.db1000g37
//     params.msimodel = params.msimodelhg19
//     params.assembly = 'GRCh37'
//     params.cachee = params.cache_dir
//     params.dirplugin = params.dirPlugin
//     params.dbNSFP = params.dbnsfp37
//     params.CADDsnv = params.CADDsnv37
//     params.CADDindel = params.CADDindels37
//     params.dbscSNV = params.dbscSNV37
//     params.cnv_reference = params.cnv_reference37
//     params.mito_reference = params.mito_reference
//     params.mito_annotation = params.mito_annotation
//     params.exomiser_properties = params.exomiser_hg19_properties
//     params.template_yml = params.template_yml
//     params.clingen_file = params.clingen_filehg37
//     params.pli = params.pli37
//     params.vusprize_script = params.vusprize_script
//     params.dependent_script = params.dependent_script
//     params.spliceaiindel = params.spliceaiindel37
//     params.spliceaisnv = params.spliceaisnv37
//     params.phenotypesfile = params.phenotypesfile37
//     params.loeuf = params.loeuf37
//     params.clinpred = params.clinpred37
//     params.alphamissense = params.alphamissense37
//     params.primateai = params.primateai37
//     params.revel = params.revel37
//     params.clinvar = params.clinvar37
//     //params.dbnsfp4 = params.dbnsfp4_37
//     params.mappability = params.mappabilityhg37
//     // ACMG
//     params.cache_bias = params.cache_biashg19
//     params.sd_bias = params.sd_biashg19
//     params.ref_bias = params.ref_biashg19
//     params.parameter_bias = params.parameter_biashg19

// } else {
//     error "Unsupported genome reference: ${params.genome}. Please use 'hg38' or 'hg19'."
// }


// include { FASTQC } from './modules/fastqc_process.nf'
// include { PROCESS_FASTQ } from './modules/fastq_process.nf'
// include { BWA } from './modules/bwa_align.nf'
// include { DEEPTRIO } from './modules/deepTrios.nf'
// include { Denom; phasingDnm } from './modules/phasingdenovo.nf'
// include { acmgpre; acmgclass } from './modules/acmg.nf'
// include { MODIFY_CONFIG; run_prioritizer } from './modules/exomiserrios.nf'
// include { annotVEP } from './modules/9.VEPannotation.nf'
// include { trioG2Pprocessing } from './modules/join_vep_vcf.nf'
// include { fullG2PPrioritisation } from './modules/11.G2P_full.nf'
// include { TRIOCNV } from './modules/tricnv2.nf'
// include { TRIOCNV_ANNOTATION } from './modules/cnvannotation.nf'
// include { SV_CALLING } from './modules/mantasv.nf'
// include { SV_ANNOTATION } from './modules/sv_annotation.nf'
// include { mergeExomiserVEPResults } from './modules/exomiser_merge_vepvcf.nf'
// include { vusVEP } from './modules/10.process_vus.nf'
// include { processVEPOutput; runVusPrioritization } from './modules/11.process_vus_2.nf'

// // Helper function
// def hasPhenotypes() {
//     return !params.phenotypes.isEmpty()
// }

// // Simple PED file parser
// def parsePedFile(pedFile) {
//     def trios = [:]
//     new File(pedFile).eachLine { line ->
//         if (!line.startsWith("#") && line.trim()) {
//             def parts = line.trim().split(/\s+/)
//             if (parts.size() >= 4 && parts[2] != "0" && parts[3] != "0") {
//                 trios[parts[0]] = [
//                     proband: parts[1],
//                     father: parts[2],
//                     mother: parts[3]
//                 ]
//             }
//         }
//     }
//     return trios
// }

// // Explicit file grouping without complex operations
// workflow {
//     def clingen_file_ch = Channel.fromPath(params.clingen_file)
//         .map { file -> tuple("clingen", file) }
    
//     def vep_vcf_script_ch = Channel.fromPath(params.vep_vcf_script)
//         .map { file -> tuple("script", file) }
    
//     def omim_clingene_file_ch = Channel.fromPath(params.omim_clingene_file)
//         .map { file -> tuple("omim_clingene", file) }

//     def mondo_owl_ch = Channel.fromPath(params.mondo_owl)
//         .map { file -> tuple("mondo", file) }

//     def g2p_script_ch = Channel.fromPath(params.g2p_script)
//         .map { file -> tuple("g2p_script"), file }

//     // Parse PED file
//     trioMap = parsePedFile(params.triosEXOped)
    
//     // Create explicit file groups
//     Channel.from(trioMap.collect { family, members ->
//         def files = [
//             file("${params.input_dir}/${family}_${members.proband}_R1.fastq.gz"),
//             file("${params.input_dir}/${family}_${members.proband}_R2.fastq.gz"),
//             file("${params.input_dir}/${family}_${members.father}_R1.fastq.gz"),
//             file("${params.input_dir}/${family}_${members.father}_R2.fastq.gz"),
//             file("${params.input_dir}/${family}_${members.mother}_R1.fastq.gz"),
//             file("${params.input_dir}/${family}_${members.mother}_R2.fastq.gz")
//         ]
//         // Verify all files exist
//         if (files.every { it.exists() }) {
//             // Now includes members info in the tuple
//             [family, members, files]  // <-- This is the key change
//         } else {
//             System.err.println "Missing files for family ${family}"
//             null
//         }
//     })
//     .filter { it != null }
//     .set { trio_channel }

//     // Debug view
//     trio_channel.view { "PROCESSING: ${it[0]}\nFILES:\n${it[2].join('\n')}" }

//     // Starts pre-processes
//     fastqc_results_raw = FASTQC(trio_channel)
//     processed_reads = PROCESS_FASTQ(trio_channel)
//     processed_reads = PROCESS_FASTQ.out.processed_reads
    
//     //aligned_bams = BWA_ALIGN(PROCESS_FASTQ.out.processed_reads)
//     // Prepare BWA input channel
//     bwa_input = PROCESS_FASTQ.out.processed_reads.map { family, p1, p2, f1, f2, m1, m2 ->
//         [
//             family,
//             trioMap[family],  // members info
//             p1, p2,
//             f1, f2,
//             m1, m2
//         ]
//     }
    
//     // Run BWA alignment
//     BWA(bwa_input, params.ref)
    
//     // View results
//     BWA.out.bams.view { 
//         """
//         FAMILY: ${it[0]}
//         PROBAND: ${it[1]} 
//         FATHER: ${it[3]}
//         MOTHER: ${it[5]}
//         """ 
//     }
//     // Connect BWA to DeepTrio
//     BWA.out.bams
//         .combine( Channel.fromPath(params.ref) )
//         .combine( Channel.fromPath(params.bed) )
//         .set { deeptrio_input }
    
//     DEEPTRIO(deeptrio_input, params.ref, params.bed)
    
//     // Access outputs
//     DEEPTRIO.out.gvcfs.view { 
//         """
//         Variant calls for ${it[0]}:
//         - Proband: ${it[1]}
//         - Father: ${it[3]}
//         - Mother: ${it[5]}
//         """ 
//     }
//     cohortvcfs = DEEPTRIO.out.cohortvcfs

//     // Denovo mutation
//     cohortvcfs_patched = cohortvcfs.map { family_id, vcf_file -> 
//         def members = trioMap[family_id]
//         tuple(family_id, members, vcf_file)
//     }
    

//     Denom(cohortvcfs_patched, params.triosEXOped)
//     denovo_vcf = Denom.out.denovo_vcf
//     denovo_vcf.view { "PHASING VCF INPUT: ${it}" }


//     // phasing
//     BWA.out.bams
//         .combine( Channel.fromPath(params.ref) )
//         .combine( Channel.fromPath(params.bed) )
//         .set { deeptrio_input1 }

//     phasingDnm(denovo_vcf, deeptrio_input1, params.ref, params.triosEXOped)
//     final_vcf_phased = phasingDnm.out.final_vcf_phased
//     final_vcf_phased.view { "FINAL VCF PHASED: ${it}" }

//     // ACMG
//     acmgpre(final_vcf_phased, params.cache_bias, params.sd_bias, params.ref_bias)
//     acmg_pre = acmgpre.out.acmg_pre
//     acmgclass(acmg_pre, params.parameter_bias) 
//     acmg_class = acmgclass.out.acmg_class

//     // VEP annotation
//     annotVEP(
//         final_vcf_phased,
//         params.ref,
//         params.cache_dir,
//         params.dirPlugin,
//         params.dbNSFP,
//         params.LoFtool,
//         params.CADDsnv,
//         params.CADDindel,
//         params.dbscSNV,
//         params.assembly,
//         params.dosage_collins,
//         params.spliceaiindel,
//         params.spliceaisnv,
//         params.pli,
//         params.phenotypesfile,
//         params.loeuf,
//         params.clinpred,
//         params.alphamissense,
//         params.primateai,
//         params.revel,
//         params.clinvar
//     )
    
//         //vus reclassification
//     vusVEP(final_vcf_phased, params.cache_dir, params.dirPlugin, params.ref,
//               params.assembly, params.dbnsfp4, params.LoFtool,
//               params.CADDsnv, params.CADDindel, params.dbscSNV, params.maxentscan)
//     vus_vep = vusVEP.out.vus_vep

//     processVEPOutput(vus_vep)
//     processed_vus = processVEPOutput.out.processed_vcf

//     vusprize_script_ch = Channel.value(file(params.vusprize_script))
//     dependent_script_ch = Channel.value(file(params.dependent_script))
//     runVusPrioritization(processed_vus,
//                            vusprize_script_ch,
//                            dependent_script_ch)
//     vus_ch = runVusPrioritization.out.vusprize_output


//     // merging vep with vcf to get GT,DP,GQ. will be common in G2P and P2G
//     trioG2Pprocessing(
//         annotVEP.out.annotVEP_vcfs,  // VEP annotated output
//         final_vcf_phased,
//         acmg_class,
//         params.triosEXOped,
//         clingen_file_ch,
//         vep_vcf_script_ch
//     )
//    // CNV calling
//     BWA.out.bams.view { "BAMs for ${it[0]}: ${it[2]}, ${it[4]}, ${it[6]}" }
//     BWA.out.bams
//         .map { family, members, proband_bam, proband_bai, father_bam, father_bai, mother_bam, mother_bai ->
//             [family, members, proband_bam, father_bam, mother_bam]
//         }
//         .set { cnv_input }
//     cnv_input.view { "CNV INPUT: Family=${it[0]}, Proband=${it[2]}, Father=${it[3]}, Mother=${it[4]}" }

//     TRIOCNV(
//         cnv_input,
//         params.ref,
//         params.mappability,
//         params.triosEXOped,
//         params.cnv2vcf_py
//     )

//    // CNV annotation with annotsv
//     TRIOCNV_ANNOTATION(
//         TRIOCNV.out.cnv_vcf,
//         params.cnv_prioritise_py
//     )
 
//       // SV calling
//     SV_CALLING(
//         cnv_input, params.ref
//     )

//     SV_ANNOTATION(
//         SV_CALLING.out.manta_vcf,
//         params.sv_prioritise_py
//     )

//     // Ends pre-process

//     if (params.mode == 'G2P') {
//         fullG2PPrioritisation(
//             trioG2Pprocessing.out.vcfvep_output,
//             omim_clingene_file_ch,
//             mondo_owl_ch,
//             params.child_sex,
//             params.family_history,
//             g2p_script_ch
//         )
//     }

//     // --------------------------------
//     //       P2G: PHENOTYPE TO GENES
//     // --------------------------------
//     if (params.mode == 'P2G') {
//         MODIFY_CONFIG(params.template_yml, file(params.triosEXOped), params.phenotypes, params.genome)
//         modified_config = MODIFY_CONFIG.out.config
//         run_prioritizer(final_vcf_phased, modified_config, params.genome, params.exomiser_properties, params.triosEXOped, params.exomiser_jar)
//         variants_tsv_json = run_prioritizer.out.variants_tsv_json

//         mergeExomiserVEPResults(
//             variants_tsv_json,
//             trioG2Pprocessing.out.vcfvep_output,
//             params.p2g_exomisermerging_script
//         )
//     }
// }

// /////////////////// Version 2 worked with mode and skip_snv, skip_cnv, skip_sv parameters ////////////////////////////////


// nextflow.enable.dsl=2

// // Parameters
// params.genome = params.genome ?: 'hg38'  // Default to hg38 if not specified
// params.platform = params.platform ?: 'illumina_pe'  // Default to Illumina paired-end
// params.input_type = null
// params.vcf_type = null
// params.refhg38 = '/usr/src/app/ref38/hg38/hg381_22XYM/Homo_sapiens_assembly38cleaned.fasta' //chr1_22 X_Y_M only
// params.refhg37 = '/usr/src/app/ref38/hg19/hg19122XYM/hg19122XYM.fa' //chr1_22 X_Y_M only
// params.bedhg38 = params.bedhg38 ?:'/usr/src/app/ref38/BED/hg38_exome.bed'
// params.bedhg37 = params.bedhg37 ?:'/usr/src/app/ref38/BED/hg37_exome.bed'
// params.cnv_reference38 = '/usr/src/app/reference/reference/cnv_cnvkit/flat_referencehg38.cnn'
// params.cnv_reference37 = '/usr/src/app/reference/reference/cnv_cnvkit/flat_referencehg19_cleaned.cnn'
// params.mito_reference = '/usr/src/app/reference/reference/mitochondria/rCRS.fasta'
// params.mito_annotation = '/usr/src/app/reference/reference/mitochondria/rCRS_annotation_2020-08-20.txt'
// params.vusprize_script = '/usr/src/app/VusPrize/vusprize/VusPrize.py'
// params.dependent_script = '/usr/src/app/VusPrize/vusprize/RF_niu.joblib'

// //new
// params.secondaryfindings_gene_file = '/usr/src/app/modules/ACMG_SF3_gene_list.xlsx'
// params.clingen_file38 = '/usr/src/app/modules/ClinGen_gene_curation_list_GRCh38_cleaned.tsv'
// params.clingen_file37 = '/usr/src/app/modules/ClinGen_gene_curation_list_GRCh37_cleaned.tsv'
// params.gene_data_file = '/usr/src/app/modules/gene_data_with_lof_class.csv' //taken from clingen gene-disease association
// //params.gender = params.gender ?: 'male'
// params.moi_dk = '/usr/src/app/modules/gene.tsv'
// params.omim = '/usr/src/app/modules/merge_9314_final_formatted.csv'

// //// VEP input Databases Hg38
// params.cache_dir = '/usr/src/app/vepC/cache'
// params.dirPlugin = '/usr/src/app/ensembl-vep/Plugins'
// params.dbnsfp38 = '/usr/src/app/vepDB/DBs/hg38/dbNSFP5.1a_grch38-004.gz'
// params.dbnsfp4_38 = '/usr/src/app/vepDB/DBs/hg38/dbNSFP4.7a_grch38.gz'
// params.LoFtool = '/usr/src/app/vepPlugins/LoFtool_scores.txt'
// params.CADDsnv38 = '/usr/src/app/vepDB/DBs/hg38/whole_genome_SNVs.tsv.gz'
// params.CADDindels38 = '/usr/src/app/vepDB/DBs/hg38/gnomad.genomes.r4.0.indel_inclAnno.tsv.gz'
// params.dbscSNV38 = '/usr/src/app/vepDB/DBs/hg38/dbscSNV1.1_GRCh38.txt.gz'
// params.dosage_collins = '/usr/src/app/vepPlugins/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz'
// params.spliceaiindel38 = '/usr/src/app/vepDB/DBs/hg38/spliceai_scores.raw.indel.hg38.vcf.gz'
// params.spliceaisnv38 = '/usr/src/app/vepDB/DBs/hg38/spliceai_scores.raw.snv.hg38.vcf.gz'
// params.pli38 = '/usr/src/app/vepDB/DBs/hg38/plI_gene.txt'
// params.phenotypesfile38 = '/usr/src/app/vepDB/DBs/hg38/Phenotypes.pm_homo_sapiens_114_GRCh38.gvf.gz'
// params.loeuf38 = '/usr/src/app/vepDB/DBs/hg38/loeuf_dataset_grch38.tsv.gz'
// params.clinpred38 = '/usr/src/app/vepDB/DBs/hg38/ClinPred_hg38_sorted_tabbed.tsv.gz'
// params.alphamissense38 = '/usr/src/app/vepDB/DBs/hg38/AlphaMissense_hg38.tsv.gz'
// params.primateai38 = '/usr/src/app/vepDB/DBs/hg38/PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz'
// params.revel38 = '/usr/src/app/vepDB/DBs/hg38/new_tabbed_revel_grch38.tsv.gz'
// params.clinvar38 = '/usr/src/app/vepDB/DBs/hg38/clinvar.vcf.gz'

// //// VEP input Databases Hg37
// params.dbnsfp37 ='/usr/src/app/vepDB/DBs/hg19/dbNSFP5.1a_grch37-003.gz'
// params.dbnsfp4_37 = '/usr/src/app/vepDB/DBs/hg19/dbNSFP4.7a_grch37.gz'
// params.CADDsnv37 = '/usr/src/app/vepDB/DBs/hg19/whole_genome_SNVs_hg37.tsv.gz'
// params.CADDindels37 = '/usr/src/app/vepDB/DBs/hg19/gnomad.genomes-exomes.r4.0.indel_inclAnno_hg37.tsv.gz'
// params.dbscSNV37 = '/usr/src/app/vepDB/DBs/hg19/dbscSNV1.1_GRCh37.txt.gz'
// params.LoFtool = '/usr/src/app/vepPlugins/LoFtool_scores.txt'
// params.CADDsnv37 = '/usr/src/app/vepDB/DBs/hg19/whole_genome_SNVs_hg37.tsv.gz'
// params.CADDindels37 = '/usr/src/app/vepDB/DBs/hg19/gnomad.genomes.r4.0.indel_inclAnno_hg37.tsv.gz'
// params.dbscSNV37 = '/usr/src/app/vepDB/DBs/hg19/dbscSNV1.1_GRCh37.txt.gz'
// params.dosage_collins = '/usr/src/app/vepPlugins/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz'
// params.spliceaiindel37 = '/usr/src/app/vepDB/DBs/hg19/spliceai_scores.raw.indel.hg19.vcf.gz'
// params.spliceaisnv37 = '/usr/src/app/vepDB/DBs/hg19/spliceai_scores.raw.snv.hg19.vcf.gz'
// params.pli37 = '/usr/src/app/vepDB/DBs/hg19/plI_gene.txt'
// params.phenotypesfile37 = '/usr/src/app/vepDB/DBs/hg19/Phenotypes.pm_homo_sapiens_114_GRCh37.gvf.gz'
// params.loeuf37 = '/usr/src/app/vepDB/DBs/hg19/loeuf_dataset_hg37.tsv.gz'
// params.clinpred37 = '/usr/src/app/vepDB/DBs/hg19/ClinPred_tabbed.tsv.gz'
// params.alphamissense37 = '/usr/src/app/vepDB/DBs/hg19/AlphaMissense_hg19.tsv.gz'
// params.primateai37 = '/usr/src/app/vepDB/DBs/hg19/PrimateAI_scores_v0.2_hg19.tsv.bgz'
// params.revel37 = '/usr/src/app/vepDB/DBs/hg19/new_tabbed_revel_grch37.tsv.gz'
// params.clinvar37 = '/usr/src/app/vepDB/DBs/hg19/clinvar.vcf.gz'

// /// Exomiser parameters
// params.template_yml = params.template_yml ?:'/usr/src/app/exomiser/db/exomiser-template.yml'  // Path inside the container
// params.phenotypes = []                        // List of phenotypes
// params.exomiser_jar = '/usr/src/app/exomiser/exomiser-cli-14.0.0/exomiser-cli-14.0.0.jar'   // Correct jar path in Docker
// params.exomiser_hg38_properties = '/usr/src/app/vepDB/DBs/exomiser/hg38application.properties' // Properties path
// params.exomiser_hg19_properties = '/usr/src/app/vepDB/DBs/exomiser/hg37application.properties'

// /// CNV baseline sample dir
// params.cnv_baseline = '/usr/src/app/ref38/CNV_baseline_samples/hg38/'


//     // ACMG
// params.cache_biashg38 = '/usr/src/app/ref38/ACMG/data/GRCh38/Cache/GRCh38/Both'
// params.sd_biashg38 = '/usr/src/app/ref38/ACMG/data/GRCh38/SupplementaryAnnotation/GRCh38/'
// params.ref_biashg38 = '/usr/src/app/ref38/ACMG/data/GRCh38/References/Homo_sapiens.GRCh38.Nirvana.dat'
// params.parameter_biashg38 = '/usr/src/app/ref38/ACMG/hg38_required_paths.json'
    
// params.cache_biashg19 = '/usr/src/app/ref38/ACMG/data/GRCh37/Cache/GRCh37/Both'
// params.sd_biashg19 = '/usr/src/app/ref38/ACMG/data/GRCh37/SupplementaryAnnotation/GRCh37/'
// params.ref_biashg19 = '/usr/src/app/ref38/ACMG/data/GRCh37/References/Homo_sapiens.GRCh37.Nirvana.dat'
// params.parameter_biashg19 = '/usr/src/app/ref38/ACMG/hg19_required_paths.json'


// // Genome reference selection
// if (params.genome == 'hg38') {
//     params.ref = params.refhg38
//     params.bed = params.bedhg38
//     params.assembly = 'GRCh38'
//     //params.cachee = params.cache_dir
//     params.dirPlugin = params.dirplugin
//     params.dbNSFP = params.dbNSFP38
//     params.CADDsnv = params.CADDsnv38
//     params.CADDindel = params.CADDindels38
//     params.dbscSNV = params.dbscSNV38
//     params.exomiser_properties = params.exomiser_propertieshg38
//     //params.vusprize_script = params.vusprize_script
//     //params.dependent_script = params.dependent_script
//     params.spliceaiindel = params.spliceaiindel38
//     params.spliceaisnv = params.spliceaisnv38
//     params.phenotypesfile = params.phenotypesfile38
//     params.loeuf = params.loeuf38
//     params.clinpred = params.clinpred38
//     params.alphamissense = params.alphamissense38
//     params.primateai = params.primateai38
//     params.revel = params.revel38
//     params.clinvar = params.clinvar38
//     params.clingen_file = params.clingen_filehg38
//     params.pli = params.plihg38
//     params.mappability = params.mappabilityhg38
//     params.dbnsfp4 = params.dbnsfp4_38
//     // ACMG
//     params.cache_bias = params.cache_biashg38
//     params.sd_bias = params.sd_biashg38
//     params.ref_bias = params.ref_biashg38
//     params.parameter_bias = params.parameter_biashg38

// } else if (params.genome == 'hg19') {
//     params.ref = params.refhg37
//     params.ref_dir = params.ref_dirhg37
//     params.reffai = params.refhg37fai
//     params.bed = params.bedhg37
//     params.gnomad = params.gnomad37
//     params.db1000g = params.db1000g37
//     params.msimodel = params.msimodelhg19
//     params.assembly = 'GRCh37'
//     params.cachee = params.cache_dir
//     params.dirplugin = params.dirPlugin
//     params.dbNSFP = params.dbnsfp37
//     params.CADDsnv = params.CADDsnv37
//     params.CADDindel = params.CADDindels37
//     params.dbscSNV = params.dbscSNV37
//     params.cnv_reference = params.cnv_reference37
//     params.mito_reference = params.mito_reference
//     params.mito_annotation = params.mito_annotation
//     params.exomiser_properties = params.exomiser_hg19_properties
//     params.template_yml = params.template_yml
//     params.clingen_file = params.clingen_filehg37
//     params.pli = params.pli37
//     params.vusprize_script = params.vusprize_script
//     params.dependent_script = params.dependent_script
//     params.spliceaiindel = params.spliceaiindel37
//     params.spliceaisnv = params.spliceaisnv37
//     params.phenotypesfile = params.phenotypesfile37
//     params.loeuf = params.loeuf37
//     params.clinpred = params.clinpred37
//     params.alphamissense = params.alphamissense37
//     params.primateai = params.primateai37
//     params.revel = params.revel37
//     params.clinvar = params.clinvar37
//     //params.dbnsfp4 = params.dbnsfp4_37
//     params.mappability = params.mappabilityhg37
//     // ACMG
//     params.cache_bias = params.cache_biashg19
//     params.sd_bias = params.sd_biashg19
//     params.ref_bias = params.ref_biashg19
//     params.parameter_bias = params.parameter_biashg19

// } else {
//     error "Unsupported genome reference: ${params.genome}. Please use 'hg38' or 'hg19'."
// }


// include { FASTQC } from './modules/fastqc_process.nf'
// include { PROCESS_FASTQ } from './modules/fastq_process.nf'
// include { BWA } from './modules/bwa_align.nf'
// include { DEEPTRIO } from './modules/deepTrios.nf'
// include { Denom; phasingDnm } from './modules/phasingdenovo.nf'
// include { acmgpre; acmgclass } from './modules/acmg.nf'
// include { MODIFY_CONFIG; run_prioritizer } from './modules/exomiserrios.nf'
// include { annotVEP } from './modules/9.VEPannotation.nf'
// include { trioG2Pprocessing } from './modules/join_vep_vcf.nf'
// include { fullG2PPrioritisation } from './modules/11.G2P_full.nf'
// include { TRIOCNV } from './modules/tricnv2.nf'
// include { TRIOCNV_ANNOTATION } from './modules/cnvannotation.nf'
// include { SV_CALLING } from './modules/mantasv.nf'
// include { SV_ANNOTATION } from './modules/sv_annotation.nf'
// include { mergeExomiserVEPResults } from './modules/exomiser_merge_vepvcf.nf'
// include { vusVEP } from './modules/10.process_vus.nf'
// include { processVEPOutput; runVusPrioritization } from './modules/11.process_vus_2.nf'

// // Helper function
// def hasPhenotypes() {
//     return !params.phenotypes.isEmpty()
// }

// // Simple PED file parser
// def parsePedFile(pedFile) {
//     def trios = [:]
//     new File(pedFile).eachLine { line ->
//         if (!line.startsWith("#") && line.trim()) {
//             def parts = line.trim().split(/\s+/)
//             if (parts.size() >= 4 && parts[2] != "0" && parts[3] != "0") {
//                 trios[parts[0]] = [
//                     proband: parts[1],
//                     father: parts[2],
//                     mother: parts[3]
//                 ]
//             }
//         }
//     }
//     return trios
// }

// // Explicit file grouping without complex operations
// workflow {
//     def clingen_file_ch = Channel.fromPath(params.clingen_file)
//         .map { file -> tuple("clingen", file) }
    
//     def vep_vcf_script_ch = Channel.fromPath(params.vep_vcf_script)
//         .map { file -> tuple("script", file) }
    
//     def omim_clingene_file_ch = Channel.fromPath(params.omim_clingene_file)
//         .map { file -> tuple("omim_clingene", file) }

//     def mondo_owl_ch = Channel.fromPath(params.mondo_owl)
//         .map { file -> tuple("mondo", file) }

//     def g2p_script_ch = Channel.fromPath(params.g2p_script)
//         .map { file -> tuple("g2p_script"), file }

//     // Parse PED file
//     trioMap = parsePedFile(params.triosEXOped)
    
//     // Create explicit file groups
//     Channel.from(trioMap.collect { family, members ->
//         def files = [
//             file("${params.input_dir}/${family}_${members.proband}_R1.fastq.gz"),
//             file("${params.input_dir}/${family}_${members.proband}_R2.fastq.gz"),
//             file("${params.input_dir}/${family}_${members.father}_R1.fastq.gz"),
//             file("${params.input_dir}/${family}_${members.father}_R2.fastq.gz"),
//             file("${params.input_dir}/${family}_${members.mother}_R1.fastq.gz"),
//             file("${params.input_dir}/${family}_${members.mother}_R2.fastq.gz")
//         ]
//         // Verify all files exist
//         if (files.every { it.exists() }) {
//             // Now includes members info in the tuple
//             [family, members, files]  // <-- This is the key change
//         } else {
//             System.err.println "Missing files for family ${family}"
//             null
//         }
//     })
//     .filter { it != null }
//     .set { trio_channel }

//     // Debug view
//     trio_channel.view { "PROCESSING: ${it[0]}\nFILES:\n${it[2].join('\n')}" }

//     // Starts pre-processes
//     fastqc_results_raw = FASTQC(trio_channel)
//     processed_reads = PROCESS_FASTQ(trio_channel)
//     processed_reads = PROCESS_FASTQ.out.processed_reads
    
//     //aligned_bams = BWA_ALIGN(PROCESS_FASTQ.out.processed_reads)
//     // Prepare BWA input channel
//     bwa_input = PROCESS_FASTQ.out.processed_reads.map { family, p1, p2, f1, f2, m1, m2 ->
//         [
//             family,
//             trioMap[family],  // members info
//             p1, p2,
//             f1, f2,
//             m1, m2
//         ]
//     }
    
//     // Run BWA alignment
//     BWA(bwa_input, params.ref)

//     // ----------------------------
//     //       OPTIONAL: SNV
//     // ----------------------------
//     if (!params.skip_snv) {

//         BWA.out.bams.combine(Channel.fromPath(params.ref)).combine(Channel.fromPath(params.bed)).set { deeptrio_input }
//         DEEPTRIO(deeptrio_input, params.ref, params.bed)

//         cohortvcfs = DEEPTRIO.out.cohortvcfs
//         cohortvcfs_patched = cohortvcfs.map { family_id, vcf_file -> tuple(family_id, trioMap[family_id], vcf_file) }
//         Denom(cohortvcfs_patched, params.triosEXOped)
//         denovo_vcf = Denom.out.denovo_vcf

//         BWA.out.bams.combine(Channel.fromPath(params.ref)).combine(Channel.fromPath(params.bed)).set { deeptrio_input1 }
//         phasingDnm(denovo_vcf, deeptrio_input1, params.ref, params.triosEXOped)
//         final_vcf_phased = phasingDnm.out.final_vcf_phased

//         acmgpre(final_vcf_phased, params.cache_bias, params.sd_bias, params.ref_bias)
//         acmg_pre = acmgpre.out.acmg_pre

//         acmgclass(acmg_pre, params.parameter_bias)
//         acmg_class = acmgclass.out.acmg_class

//         annotVEP(final_vcf_phased,
//                  params.ref,
//                  params.cache_dir,
//                  params.dirPlugin,
//                  params.dbNSFP,
//                  params.LoFtool,
//                  params.CADDsnv,
//                  params.CADDindel,
//                  params.dbscSNV,
//                  params.assembly,
//                  params.dosage_collins,
//                  params.spliceaiindel,
//                  params.spliceaisnv,
//                  params.pli,
//                  params.phenotypesfile,
//                  params.loeuf,
//                  params.clinpred,
//                  params.alphamissense,
//                  params.primateai,
//                  params.revel,
//                  params.clinvar)

//         vusVEP(final_vcf_phased, params.cache_dir, params.dirPlugin, params.ref,
//                params.assembly, params.dbnsfp4, params.LoFtool, params.CADDsnv,
//                params.CADDindel, params.dbscSNV, params.maxentscan)
//         vus_vep = vusVEP.out.vus_vep
//         processVEPOutput(vus_vep)
//         processed_vus = processVEPOutput.out.processed_vcf

//         vusprize_script_ch = Channel.value(file(params.vusprize_script))
//         dependent_script_ch = Channel.value(file(params.dependent_script))
        
//         runVusPrioritization(processed_vus,
//                             vusprize_script_ch,
//                             dependent_script_ch)

//         vus_ch = runVusPrioritization.out.vusprize_output


//         trioG2Pprocessing(
//             annotVEP.out.annotVEP_vcfs,
//             final_vcf_phased,
//             acmg_class,
//             params.triosEXOped,
//             clingen_file_ch,
//             vep_vcf_script_ch)
//     }

//     // ----------------------------
//     //       OPTIONAL: CNV
//     // ----------------------------
//     if (!params.skip_cnv) {
//         cnv_input = BWA.out.bams.map { family, members, proband_bam, proband_bai, father_bam, father_bai, mother_bam, mother_bai -> [family, members, proband_bam, father_bam, mother_bam] }
//         TRIOCNV(cnv_input, params.ref, params.mappability, params.triosEXOped, params.cnv2vcf_py)
//         TRIOCNV_ANNOTATION(TRIOCNV.out.cnv_vcf, params.cnv_prioritise_py)

//     }

//     // ----------------------------
//     //       OPTIONAL: SV
//     // ----------------------------
//     if (!params.skip_sv) {
//         sv_input = BWA.out.bams.map { family, members, proband_bam, proband_bai, father_bam, father_bai, mother_bam, mother_bai -> [family, members, proband_bam, father_bam, mother_bam] }
//         SV_CALLING(sv_input, params.ref)
//         SV_ANNOTATION(SV_CALLING.out.manta_vcf, params.sv_prioritise_py)
//     }

//     // ----------------------------
//     //       G2P Mode
//     // ----------------------------
//     if (params.mode == 'G2P' && !params.skip_snv) {
//         fullG2PPrioritisation(
//             trioG2Pprocessing.out.vcfvep_output,
//             omim_clingene_file_ch,
//             mondo_owl_ch,
//             params.child_sex,
//             params.family_history,
//             g2p_script_ch
//         )
//     }

//     // ----------------------------
//     //       P2G Mode
//     // ----------------------------
//     if (params.mode == 'P2G' && !params.skip_snv) {
//         MODIFY_CONFIG(params.template_yml, file(params.triosEXOped), params.phenotypes, params.genome)
//         modified_config = MODIFY_CONFIG.out.config

//         run_prioritizer(final_vcf_phased, modified_config, params.genome, params.exomiser_properties, params.triosEXOped, params.exomiser_jar)
//         mergeExomiserVEPResults(
//             run_prioritizer.out.variants_tsv_json,
//             trioG2Pprocessing.out.vcfvep_output,
//             params.p2g_exomisermerging_script
//         )
//     }
// }


/////////////////// Version 3  ////////////////////////////////


nextflow.enable.dsl=2


// Genome reference selection
if (params.genome == 'hg38') {
    params.ref = params.refhg38
    params.bed = params.bedhg38
    params.assembly = 'GRCh38'
    //params.cachee = params.cache_dir
    params.dirPlugin = params.dirplugin
    params.dbNSFP = params.dbNSFP38
    params.CADDsnv = params.CADDsnv38
    params.CADDindel = params.CADDindels38
    params.dbscSNV = params.dbscSNV38
    params.exomiser_properties = params.exomiser_propertieshg38
    //params.vusprize_script = params.vusprize_script
    //params.dependent_script = params.dependent_script
    params.spliceaiindel = params.spliceaiindel38
    params.spliceaisnv = params.spliceaisnv38
    params.phenotypesfile = params.phenotypesfile38
    params.loeuf = params.loeuf38
    params.clinpred = params.clinpred38
    params.alphamissense = params.alphamissense38
    params.primateai = params.primateai38
    params.revel = params.revel38
    params.clinvar = params.clinvar38
    params.clingen_file = params.clingen_filehg38
    params.pli = params.plihg38
    params.mappability = params.mappabilityhg38
    params.dbnsfp4 = params.dbnsfp4_38
    // ACMG
    params.cache_bias = params.cache_biashg38
    params.sd_bias = params.sd_biashg38
    params.ref_bias = params.ref_biashg38
    params.parameter_bias = params.parameter_biashg38

} else if (params.genome == 'hg19') {
    params.ref = params.refhg37
    params.ref_dir = params.ref_dirhg37
    params.reffai = params.refhg37fai
    params.bed = params.bedhg37
    params.gnomad = params.gnomad37
    params.db1000g = params.db1000g37
    params.msimodel = params.msimodelhg19
    params.assembly = 'GRCh37'
    params.cachee = params.cache_dir
    params.dirPlugin = params.dirPlugin
    params.dbNSFP = params.dbNSFP37
    params.CADDsnv = params.CADDsnv37
    params.CADDindel = params.CADDindels37
    params.dbscSNV = params.dbscSNV37
    params.cnv_reference = params.cnv_reference37
    params.mito_reference = params.mito_reference
    params.mito_annotation = params.mito_annotation
    params.exomiser_properties = params.exomiser_propertieshg37
    params.template_yml = params.template_yml
    params.clingen_file = params.clingen_filehg37
    params.pli = params.plihg37
    params.vusprize_script = params.vusprize_script
    params.dependent_script = params.dependent_script
    params.spliceaiindel = params.spliceaiindel37
    params.spliceaisnv = params.spliceaisnv37
    params.phenotypesfile = params.phenotypesfile37
    params.loeuf = params.loeuf37
    params.clinpred = params.clinpred37
    params.alphamissense = params.alphamissense37
    params.primateai = params.primateai37
    params.revel = params.revel37
    params.clinvar = params.clinvar37
    params.dbnsfp4 = params.dbnsfp4_37
    params.mappability = params.mappabilityhg37
    // ACMG
    params.cache_bias = params.cache_biashg19
    params.sd_bias = params.sd_biashg19
    params.ref_bias = params.ref_biashg19
    params.parameter_bias = params.parameter_biashg19

} else {
    error "Unsupported genome reference: ${params.genome}. Please use 'hg38' or 'hg19'."
}


include { FASTQC } from './modules/fastqc_process.nf'
include { PROCESS_FASTQ } from './modules/fastq_process.nf'
include { BWA } from './modules/bwa_align.nf'
include { DEEPTRIO } from './modules/deepTrios.nf'
include { Denom; phasingDnm } from './modules/phasingdenovo.nf'
include { acmgpre; acmgclass } from './modules/acmg.nf'
include { MODIFY_CONFIG; run_prioritizer } from './modules/exomiserrios.nf'
include { annotVEP } from './modules/9.VEPannotation.nf'
include { trioG2Pprocessing } from './modules/join_vep_vcf.nf'
include { fullG2PPrioritisation } from './modules/11.G2P_full.nf'
include { TRIOCNV } from './modules/tricnv2.nf'
include { TRIOCNV_ANNOTATION } from './modules/cnvannotation.nf'
include { SV_CALLING } from './modules/mantasv.nf'
include { SV_ANNOTATION } from './modules/sv_annotation.nf'
include { mergeExomiserVEPResults } from './modules/exomiser_merge_vepvcf.nf'
include { vusVEP } from './modules/10.process_vus.nf'
include { processVEPOutput; runVusPrioritization } from './modules/11.process_vus_2.nf'

// Helper function
def hasPhenotypes() {
    return !params.phenotypes.isEmpty()
}

// Simple PED file parser
def parsePedFile(pedFile) {
    def trios = [:]
    new File(pedFile).eachLine { line ->
        if (!line.startsWith("#") && line.trim()) {
            def parts = line.trim().split(/\s+/)
            if (parts.size() >= 4 && parts[2] != "0" && parts[3] != "0") {
                trios[parts[0]] = [
                    proband: parts[1],
                    father: parts[2],
                    mother: parts[3]
                ]
            }
        }
    }
    return trios
}

// Explicit file grouping without complex operations
workflow {
    def clingen_file_ch = Channel.fromPath(params.clingen_file)
        .map { file -> tuple("clingen", file) }
    
    def vep_vcf_script_ch = Channel.fromPath(params.vep_vcf_script)
        .map { file -> tuple("script", file) }
    
    def omim_clingene_file_ch = Channel.fromPath(params.omim_clingene_file)
        .map { file -> tuple("omim_clingene", file) }

    def mondo_owl_ch = Channel.fromPath(params.mondo_owl)
        .map { file -> tuple("mondo", file) }

    def g2p_script_ch = Channel.fromPath(params.g2p_script)
        .map { file -> tuple("g2p_script"), file }

    // Parse PED file
    trioMap = parsePedFile(params.triosEXOped)
    
    // Create explicit file groups
    Channel.from(trioMap.collect { family, members ->
        def files = [
            file("${params.input_dir}/${family}_${members.proband}_R1.fastq.gz"),
            file("${params.input_dir}/${family}_${members.proband}_R2.fastq.gz"),
            file("${params.input_dir}/${family}_${members.father}_R1.fastq.gz"),
            file("${params.input_dir}/${family}_${members.father}_R2.fastq.gz"),
            file("${params.input_dir}/${family}_${members.mother}_R1.fastq.gz"),
            file("${params.input_dir}/${family}_${members.mother}_R2.fastq.gz")
        ]
        // Verify all files exist
        if (files.every { it.exists() }) {
            // Now includes members info in the tuple
            [family, members, files]  // <-- This is the key change
        } else {
            System.err.println "Missing files for family ${family}"
            null
        }
    })
    .filter { it != null }
    .set { trio_channel }

    // -------------
    // FASTQ input
    // -------------
    if (params.input_type == 'fastq') {
        fastqc_results_raw = FASTQC(trio_channel)
        processed_reads = PROCESS_FASTQ(trio_channel)
        processed_reads = PROCESS_FASTQ.out.processed_reads
        bwa_input = processed_reads.map { family, p1, p2, f1, f2, m1, m2 -> [family, trioMap[family], p1, p2, f1, f2, m1, m2] }
        BWA(bwa_input, params.ref)
    }

    // -------------
    // BAM input
    // -------------
    if (params.input_type == 'bam') {
        // Skip FASTQ → BAM; just read BAMs
        BWA.out.bams = Channel.from(trioMap.collect { family, members ->
            [
                family,
                members,
                file("${params.input_dir}/${family}_${members.proband}.bam"),
                file("${params.input_dir}/${family}_${members.proband}.bai"),
                file("${params.input_dir}/${family}_${members.father}.bam"),
                file("${params.input_dir}/${family}_${members.father}.bai"),
                file("${params.input_dir}/${family}_${members.mother}.bam"),
                file("${params.input_dir}/${family}_${members.mother}.bai")
            ]
        })
    }

    // -------------
    // VCF input
    // -------------
    if (params.input_type == 'vcf') {
    final_vcf_phased = Channel
        .fromPath("${params.input_dir}/*_cohort.dnm.phas.vcf.gz")
        .map { file_path ->
            def fname = file_path.getName()  // e.g., fam1_SRR26402599_cohort.dnm.phas.vcf.gz
            def family_id = fname.tokenize('_')[0]  // Extracts 'fam1'
            tuple(family_id, file_path)
        }
    }
    
    
    // ------------------------------------------
    // SNV: Only if not skipped and not vcf input
    // ------------------------------------------
    if (!params.skip_snv && params.input_type in ['fastq', 'bam']) {
        BWA.out.bams.combine(Channel.fromPath(params.ref)).combine(Channel.fromPath(params.bed)).set { deeptrio_input }
        DEEPTRIO(deeptrio_input, params.ref, params.bed)

        cohortvcfs = DEEPTRIO.out.cohortvcfs
        cohortvcfs_patched = cohortvcfs.map { family_id, vcf_file -> tuple(family_id, trioMap[family_id], vcf_file) }

        Denom(cohortvcfs_patched, params.triosEXOped)
        denovo_vcf = Denom.out.denovo_vcf

        BWA.out.bams.combine(Channel.fromPath(params.ref)).combine(Channel.fromPath(params.bed)).set { deeptrio_input1 }
        phasingDnm(denovo_vcf, deeptrio_input1, params.ref, params.triosEXOped)
        final_vcf_phased = phasingDnm.out.final_vcf_phased
    }

    // ------------------------------------------
    // SNV Downstream: If VCF is present
    // ------------------------------------------
    if (!params.skip_snv && final_vcf_phased) {
        acmgpre(final_vcf_phased, params.cache_bias, params.sd_bias, params.ref_bias)
        acmg_pre = acmgpre.out.acmg_pre
        acmgclass(acmg_pre, params.parameter_bias) 
        acmg_class = acmgclass.out.acmg_class

        annotVEP(final_vcf_phased,
            params.ref,
            params.cache_dir,
            params.dirPlugin,
            params.dbNSFP,
            params.LoFtool,
            params.CADDsnv,
            params.CADDindel,
            params.dbscSNV,
            params.assembly,
            params.dosage_collins,
            params.spliceaiindel,
            params.spliceaisnv,
            params.pli,
            params.phenotypesfile,
            params.loeuf,
            params.clinpred,
            params.alphamissense,
            params.primateai,
            params.revel,
            params.clinvar)

        //         //vus reclassification
        vusVEP(final_vcf_phased, params.cache_dir, params.dirPlugin, params.ref,
                params.assembly, params.dbnsfp4, params.LoFtool,
                params.CADDsnv, params.CADDindel, params.dbscSNV, params.maxentscan)
        vus_vep = vusVEP.out.vus_vep

        processVEPOutput(vus_vep)
        processed_vus = processVEPOutput.out.processed_vcf

        vusprize_script_ch = Channel.value(file(params.vusprize_script))
        dependent_script_ch = Channel.value(file(params.dependent_script))
        runVusPrioritization(processed_vus,
                               vusprize_script_ch,
                               dependent_script_ch)
        vus_ch = runVusPrioritization.out.vusprize_output

        trioG2Pprocessing(
            annotVEP.out.annotVEP_vcfs,
            final_vcf_phased,
            acmg_class,
            params.triosEXOped,
            clingen_file_ch,
            vep_vcf_script_ch)
    }

    // ------------------
    // CNV & SV conditionals
    // ------------------
    if (!params.skip_cnv && params.input_type in ['fastq', 'bam']) {
        cnv_input = BWA.out.bams.map { family, members, proband_bam, proband_bai, father_bam, father_bai, mother_bam, mother_bai -> [family, members, proband_bam, father_bam, mother_bam] }
        TRIOCNV(cnv_input, params.ref, params.mappability, params.triosEXOped, params.cnv2vcf_py)
        TRIOCNV_ANNOTATION(TRIOCNV.out.cnv_vcf, params.cnv_prioritise_py)
    }

    if (!params.skip_sv && params.input_type in ['fastq', 'bam']) {
        sv_input = BWA.out.bams.map { family, members, proband_bam, proband_bai, father_bam, father_bai, mother_bam, mother_bai -> [family, members, proband_bam, father_bam, mother_bam] }
        SV_CALLING(sv_input, params.ref)
        SV_ANNOTATION(SV_CALLING.out.manta_vcf, params.sv_prioritise_py)
    }

    // --------------
    // G2P Mode
    // --------------
    if (params.mode == 'G2P' && !params.skip_snv) {
        fullG2PPrioritisation(
            trioG2Pprocessing.out.vcfvep_output,
            omim_clingene_file_ch,
            mondo_owl_ch,
            params.child_sex,
            params.family_history,
            g2p_script_ch
        )
    }

    // --------------
    // P2G Mode
    // --------------
    if (params.mode == 'P2G' && !params.skip_snv) {
        MODIFY_CONFIG(params.template_yml, file(params.triosEXOped), params.phenotypes, params.genome)
        modified_config = MODIFY_CONFIG.out.config

        run_prioritizer(final_vcf_phased, modified_config, params.genome, params.exomiser_properties, params.triosEXOped, params.exomiser_jar)
        mergeExomiserVEPResults(
            run_prioritizer.out.variants_tsv_json,
            trioG2Pprocessing.out.vcfvep_output,
            params.p2g_exomisermerging_script
        )
    }
}    
import pandas as pd
import numpy as np
import re
import subprocess
import argparse
import os

def run_script(args, msg):
    print(f"Running: {' '.join(args)}")
    result = subprocess.run(args, capture_output=True, text=True)
    print(msg)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr)
        raise RuntimeError(f"Script {args[0]} failed!")

def is_in_domain(row):
    if pd.isna(row['Interpro_domain']) or pd.isna(row['Protein_position']):
        return False
    try:
        pos = float(row['Protein_position'])
        domains = []
        for entry in str(row['Interpro_domain']).split(','):
            if ':' in entry and '-' in entry.split(':')[2]:
                start, end = map(float, entry.split(':')[2].split('-'))
                domains.append((start, end))
        return any(start <= pos <= end for (start, end) in domains)
    except Exception:
        return False

def annotate_missense_acmg(df):
    df['PM2'] = df['gnomADe_AF'].apply(lambda af: True if pd.isna(af) or float(af) < 0.001 else False)
    df['PP3'] = (
        (df['CADD_PHRED'] >= 20) | 
        (df['REVEL_score'] >= 0.7) |
        df['Polyphen2_HDIV_pred'].str.contains('damaging', na=False) |
        df['Polyphen2_HVAR_pred'].str.contains('damaging', na=False)
    )
    df['PM1'] = df.apply(is_in_domain, axis=1)
    df['PM5'] = df['ClinVar_CLNSIG'].str.contains('pathogenic', case=False, na=False)
    df['BS1'] = df['gnomADe_AF'] > 0.05
    df['BP4'] = (
        (df['CADD_PHRED'] < 10) | 
        (df['REVEL_score'] < 0.3) |
        df['Polyphen2_HDIV_pred'].str.contains('benign', na=False)
    )
    return df

def main(args):
    # 1. Initial filtering
    df = pd.read_csv(args.vcfvep, sep="\t", engine="c")
    is_syn = df['Consequence'].str.contains(r'\bsynonymous_variant\b', na=False)
    in_splice = df['Consequence'].str.contains(r'\bsplice_region_variant\b', na=False)
    df_patho_syn = df[ is_syn & in_splice ]
    df_patho_syn_hig = df_patho_syn[df_patho_syn['IMPACT'] == 'LOW']
    df_patho_syn_hig.to_csv("synonymous_splice_region_LOW.csv", index=False)

    df = df[(df['IMPACT']=='HIGH') | (df['IMPACT']=='MODERATE') | (df['IMPACT']=='MODIFIER')]
    filter_out = [
        'intron_variant','downstream_gene_variant','upstream_gene_variant','3_prime_UTR_variant','5_prime_UTR_variant',
        'non_coding_transcript_exon_variant','intergenic_variant'
    ]
    pattern = r'\b(?:' + '|'.join(map(re.escape, filter_out)) + r')\b'
    mask_bad = df['Consequence'].str.contains(pattern, regex=True, na=False)
    df = df[~mask_bad]
    df_filtered = df[
        (df['gnomAD4.1_joint_SAS_AF'] <= 0.001)
        | df['gnomAD4.1_joint_SAS_AF'].isna()
    ]
    df_filtered.to_csv("initial_filtered.csv", index=False)

    # 2. Run primary_checking_ofvariants.py (inheritance checking)
    run_script(
        ["/usr/bin/python3", "/home/vgenomics/Germline/Pipeline/Trios_Pipeline/samples/testTrios_g2p/modules/Final_G2P_trios/primary_checking_ofvariants.py", "initial_filtered.csv", "--child_sex", args.child_sex, "--prefix", "initialfiltered_results"],
        "Ran inheritance checking"
    )

    # 3. Run add_MOI_trios.py
    run_script(
        ["/usr/bin/python3", "/home/vgenomics/Germline/Pipeline/Trios_Pipeline/samples/testTrios_g2p/modules/Final_G2P_trios/add_MOI_trios.py", "initialfiltered_results_all_zyg.csv", args.omim_clingene, args.mondo_owl, "trioPVS1_with_MOI_g2p.csv"],
        "Ran MOI annotation"
    )

    # 4. Run get_firstimportantlist_PVS1.py
    run_script(
        ["/usr/bin/python3", "/home/vgenomics/Germline/Pipeline/Trios_Pipeline/samples/testTrios_g2p/modules/Final_G2P_trios/get_firstimportantlist_PVS1.py", "trioPVS1_with_MOI_g2p.csv", "--firstlist_output_csv", args.firstlist_output_csv],
        "Ran PVS1 script"
    )
    firstlist_output_csv = args.firstlist_output_csv

    # 5. Run get_secondimportantlist.py
    #run_script(
    #    ["/usr/bin/python3", "/home/vgenomics/Germline/Pipeline/Trios_Pipeline/samples/testTrios_g2p/modules/Final_G2P_trios/get_secondimportantlist.py", "imp1.csv", "pvs1_df.csv", "--output_csv", "missense_and_splice_candidates_updated.csv"],
    #    "Ran second important list script"
    #)

    #args.firstlist_output_csv
    #  # Value from previous step
    subprocess.run([
        "/usr/bin/python3",
        "/home/vgenomics/Germline/Pipeline/Trios_Pipeline/samples/testTrios_g2p/modules/Final_G2P_trios/get_secondimportantlist.py",
        "imp1.csv",
        firstlist_output_csv,
        "--output_csv", "missense_and_splice_candidates_updated.csv"
    ], check=True)
    print("Ran second important list script")

    # 6. Run prioritise_firstimportantlist.py
    #run_script(
    #    ["/usr/bin/python3", "/home/vgenomics/Germline/Pipeline/Trios_Pipeline/samples/testTrios_g2p/modules/Final_G2P_trios/prioritise_firstimportantlist.py", "pvs1_df.csv", "--output_csv", args.output_csv],
    #    "Ran prioritisation of first important list"
    #)
    subprocess.run([
        "/usr/bin/python3",
        "/home/vgenomics/Germline/Pipeline/Trios_Pipeline/samples/testTrios_g2p/modules/Final_G2P_trios/prioritise_firstimportantlist.py",
        firstlist_output_csv,
        "--output_csv", args.output_csv
    ], check=True)
    print("Ran prioritisation of first important list")

    # 7. ACMG guideline follow for missense list
    print("Annotating ACMG codes for missense/splice candidates...")
    final = pd.read_csv("missense_and_splice_candidates_updated.csv")
    final = annotate_missense_acmg(final)
    final.to_csv(args.missense_output_csv, index=False)
    print("Annotated missense variants saved as args.missense_output_csv)")

    # 8. Run g2p prioritization script (optional, AD/AR/XL must be specified)
    if args.family_history:
        run_script(
            ["/usr/bin/python3", "/home/vgenomics/Germline/Pipeline/Trios_Pipeline/samples/testTrios_g2p/modules/Final_G2P_trios/script_g2p.py", "--family_history", args.family_history, "--input", args.missense_output_csv],
            "Ran final gene2phenotype prioritization"
        )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Full variant analysis pipeline for trio data.")
    parser.add_argument("--vcfvep", required=True, help="VEP-annotated input variant file")
    parser.add_argument("--omim_clingene", required=True, help="OMIM clingene file")
    parser.add_argument("--mondo_owl", required=True, help="MONDO .owl file")
    parser.add_argument("--child_sex", default="FEMALE", help="Child sex (MALE or FEMALE)")
    parser.add_argument("--family_history", choices=['AD','AR','XL'], default=None, help="Family history: AD, AR, XL (for final prioritization)")
    parser.add_argument("--firstlist_output_csv", default="pvs1_df.csv", help="Output CSV for all variants list (PVS1)")
    parser.add_argument("--output_csv", default="prioritized_variants.csv", help="Output CSV file (default: G2P_prioritized_variants1.csv)")
    parser.add_argument("--missense_output_csv", default="input_to_g2p_priority.csv", help="Output CSV for annotated missense variants")

    args = parser.parse_args()
    main(args)

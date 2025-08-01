import pandas as pd
import numpy as np
import argparse

def main(imp1_csv, firstlist_output_csv, output_csv):
    # Load data
    imp1 = pd.read_csv(imp1_csv)
    pvs1_df = pd.read_csv(firstlist_output_csv)

    # ——— 1) Select missense/splice candidates not in PVS1 ———
    other = imp1[~imp1['Location'].isin(pvs1_df['Location'])]
    mask_mis    = other['Consequence'].str.contains('missense_variant',       na=False)
    mask_splice = other['Consequence'].str.contains('splice_region_variant', na=False)
    cand = other[mask_mis | mask_splice].copy()

    # ——— 2) Parse & vectorize SpliceAI, exon‐distance, exon‐fraction ———
    def get_spliceai_max(s):
        if pd.isna(s):
            return 0.0
        parts = s.split('|')
        scores = []
        for idx in range(1, 5):  # DS_AG, DS_AL, DS_DG, DS_DL
            try:
                scores.append(float(parts[idx]))
            except (IndexError, ValueError):
                continue
        return max(scores) if scores else 0.0

    def get_ne_dist(s):
        if pd.isna(s):
            return np.nan
        parts = s.split('|')
        try:
            return float(parts[1])
        except (IndexError, ValueError):
            return np.nan

    # apply parsing functions
    cand['SpliceAI_max'] = cand['SpliceAI_pred'].apply(get_spliceai_max)
    cand['NE_dist']      = cand['NearestExonJB'].apply(get_ne_dist)

    # parse EXON into numeric exon number/total and fraction
    ex_split = cand['EXON'].str.split('/', expand=True)
    cand['EX_num']  = pd.to_numeric(ex_split[0], errors='coerce')
    cand['EX_tot']  = pd.to_numeric(ex_split[1], errors='coerce')
    cand['EX_frac'] = cand['EX_num'] / cand['EX_tot']

    # ——— 3) PM2: rare in population ———
    for col in ['CADD_PHRED','REVEL_score','gnomAD4.1_joint_SAS_AF']:
        cand[col] = pd.to_numeric(cand[col], errors='coerce')
        
    MAF_THR = 0.001
    cand['PM2'] = cand['gnomAD4.1_joint_SAS_AF'].fillna(0) < MAF_THR

    # ——— 4) PP3: in silico support ———
    cand['PP3'] = (
           (cand['CADD_PHRED'] >= 20)
        | (cand['REVEL_score'] >= 0.7)
        |  cand['Polyphen2_HDIV_pred'] .str.contains('damaging', na=False)
        |  cand['Polyphen2_HVAR_pred'] .str.contains('damaging', na=False)
    )

    # ——— 5) splice_support: rescue non‐canonical splice ———
    cand['splice_support'] = (
           (cand['SpliceAI_max'] >= 0.5)
        | (cand['NE_dist'] <= 3)
    )

    # ——— 6) PM1: in‐domain via InterPro_domain + Protein_position ———
    def parse_interpro(dom_str):
        ranges = []
        if pd.notna(dom_str):
            for entry in dom_str.split(','):
                parts = entry.split(':')
                # parts[2] is "start-end"
                if len(parts) >= 3 and '-' in parts[2]:
                    start, end = parts[2].split('-', 1)
                    try:
                        ranges.append((int(start), int(end)))
                    except ValueError:
                        pass
        return ranges

    cand['domain_ranges'] = cand['Interpro_domain'].apply(parse_interpro)
    cand['Protein_position'] = pd.to_numeric(cand['Protein_position'], errors='coerce')

    def in_domain(row):
        pos = row['Protein_position']
        return any(start <= pos <= end for start, end in row['domain_ranges'])

    cand['PM1'] = cand.apply(in_domain, axis=1)

    # ——— 7) PM5: exact ClinVar pathogenic match ———
    cand['PM5'] = cand['ClinVar_CLNSIG'].str.contains('pathogenic', case=False, na=False)

    # ——— 8) Build final filter ———
    # require PM2 AND at least one of {PM1, PP3, PM5, splice_support}
    keep = (
        cand['PM2'] &
        cand[['PM1','PP3','PM5','splice_support']].any(axis=1)
    )
    final = cand[keep].copy()

    # ——— 9) Rank & output ———
    final = final.sort_values(
        by=['PM2','PM1','PP3','PM5','SpliceAI_max','CADD_PHRED'],
        ascending=[False,False,False,False,False,False]
    )
    final.to_csv(output_csv, sep=',', index=False)
    print(f"Output written to: {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter missense/splice candidates for ACMG evidence.")
    parser.add_argument("imp1_csv", help="CSV file with initial impact-filtered variants (imp1)")
    parser.add_argument("firstlist_output_csv", help="CSV file with previously PVS1-annotated variants")
    parser.add_argument("--output_csv", default="missense_and_splice_candidates_updated.tsv", help="Output CSV file (default: missense_and_splice_candidates_updated.csv)")
    args = parser.parse_args()
    main(args.imp1_csv, args.firstlist_output_csv, args.output_csv)

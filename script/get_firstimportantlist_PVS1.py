import pandas as pd
import numpy as np
import argparse

def main(input_csv, firstlist_output_csv):
    # Load data
    var_moi = pd.read_csv(input_csv)

    # Filter for impactful variants
    imp1 = var_moi[(var_moi['IMPACT'] == 'HIGH') | (var_moi['IMPACT'] == 'MODERATE')]

    # Include inframe indels
    low_SO_terms = {'inframe_deletion', 'inframe_insertion'}
    contains_low = var_moi['Consequence'].str.split(',').apply(
        lambda terms: any(t in low_SO_terms for t in terms)
    )
    low_df = var_moi[contains_low]

    # Combine all variants of interest
    imp2 = pd.concat([low_df, imp1], ignore_index=True).copy()  # Explicit copy to avoid fragmentation

    # Parse NearestExonJB in one operation
    ne_cols = ['NE_trans', 'NE_dist', 'NE_region', 'NE_exonlen']
    ne_data = imp2['NearestExonJB'].str.split('|', expand=True)
    ne_data.columns = ne_cols
    ne_data['NE_dist'] = pd.to_numeric(ne_data['NE_dist'], errors='coerce')

    # Parse EXON in one operation
    exon_data = imp2['EXON'].str.split('/', expand=True)
    exon_data.columns = ['EX_num', 'EX_tot']
    exon_data = exon_data.apply(pd.to_numeric, errors='coerce')

    # Combine all parsed data at once
    imp2 = pd.concat([
        imp2.drop(columns=['NearestExonJB', 'EXON']),
        ne_data,
        exon_data
    ], axis=1)

    # Predict NMD
    def predict_nmd(row):
        if pd.isna(row['EX_num']) or pd.isna(row['EX_tot']):
            return False
        # 1) last exon → no NMD
        if row['EX_num'] == row['EX_tot']:
            return False
        # 2) penultimate exon within last 50 bp → no NMD
        if row['EX_num'] == row['EX_tot'] - 1 and row['NE_dist'] <= 50:
            return False
        # else → triggers NMD
        return True

    # HI Eligibility (3 OR 35)
    def check_hi_eligibility(score):
        if pd.isna(score):
            return False
        # Accept either 3 (standard HI) or 35 (developmental HI)
        return score in {3, 35}

    # Add all new columns at once
    new_cols = pd.DataFrame({
        'predicted_NMD': imp2.apply(predict_nmd, axis=1),
        'HI_eligible': imp2['Haploinsufficiency Score'].apply(check_hi_eligibility),
        'PVS1': np.nan  # Initialize
    })

    # Assign PVS1 after HI eligibility is set
    new_cols['PVS1'] = new_cols.apply(
        lambda row: 'PVS1_VeryStrong' if row['HI_eligible'] and row['predicted_NMD'] else 
                   'PVS1_Moderate' if row['HI_eligible'] else np.nan,
        axis=1
    )

    # Final combine
    imp2 = pd.concat([imp2, new_cols], axis=1)

    # Filter and sort PVS1 variants
    pvs1_df = imp2[imp2['PVS1'].notna()].copy()
    pvs1_df = pvs1_df.sort_values(
        by=['predicted_NMD', 'CADD_PHRED'],
        ascending=[False, False]
    ).drop_duplicates()

    # Output to CSV
    pvs1_df.to_csv(firstlist_output_csv, index=False)
    print(f"Output written to: {firstlist_output_csv}")
    imp1.to_csv("imp1.csv", index=False)
    print(f"imp1 is saved.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter and annotate PVS1 variants from a variant table.")
    parser.add_argument("input_csv", help="Path to input trioPVS1_with_MOI_g2p.csv")
    parser.add_argument("--firstlist_output_csv", default="pvs1_variants.csv", help="Output CSV file for PVS1 variants (default: pvs1_variants.csv)")
    args = parser.parse_args()
    main(args.input_csv, args.firstlist_output_csv)

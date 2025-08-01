import pandas as pd
import argparse
from typing import Optional

# --------------------------
# 1. PRIORITIZATION FUNCTIONS
# --------------------------
def prioritize_moi(moi: str, family_history: Optional[str]) -> int:
    """Convert MOI string to priority score (1=highest)."""
    if pd.isna(moi):
        return 4
    
    moi = str(moi).upper()
    
    # Family history overrides default priority
    if family_history == "AD" and ("AD" in moi or "DOMINANT" in moi):
        return 1
    elif family_history == "AR" and ("AR" in moi or "RECESSIVE" in moi):
        return 1
    elif family_history == "XL" and ("X-LINKED" in moi or "XL" in moi):
        return 1
    
    # Default priority if no family history
    if "AD" in moi:
        return 1
    elif "X-LINKED" in moi or "XL" in moi:
        return 2
    elif "AR" in moi:
        return 3
    else:
        return 4

def assign_pathogenicity_tier(row):
    """Tier variants by ClinVar status and supporting evidence."""
    clnsig = str(row['ClinVar_CLNSIG']).lower()
    
    if "pathogenic" in clnsig or "likely_pathogenic" in clnsig:
        return 1
    elif "conflicting" in clnsig and "pathogenic" in clnsig:
        return 2
    elif row['PP3'] or row['PM1'] or row['PM5'] or row.get('splice_support', False):
        return 3
    else:
        return 4

def check_moi_compatibility(row, family_history: Optional[str]) -> int:
    """Check if variant's MOI matches family history."""
    moi = str(row['MOI']).upper()
    if family_history == "AD" and ("AD" in moi or "DOMINANT" in moi):
        return 1
    elif family_history == "AR" and ("AR" in moi or "RECESSIVE" in moi):
        return 1
    elif family_history == "XL" and ("X-LINKED" in moi or "XL" in moi):
        return 1
    else:
        return 0

def assign_functional_impact(row):
    """Score functional impact (1=highest)."""
    if 'splice' in str(row['Consequence']) and row.get('SpliceAI_max', 0) >= 0.8:
        return 1
    elif 'missense' in str(row['Consequence']) and (
        row.get('CADD_PHRED', 0) >= 30 or row.get('REVEL_score', 0) >= 0.75
    ):
        return 2
    elif row.get('PM1', False):
        return 3
    else:
        return 4

def assign_frequency_score(row):
    """Score rarity (1=highest)."""
    af = row.get('gnomAD4.1_joint_SAS_AF', None)
    if pd.isna(af) or af == 0:
        return 1
    elif row.get('PM2', False):
        return 2
    else:
        return 3

# --------------------------
# 2. MAIN PIPELINE
# --------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--family_history", choices=["AD", "AR", "XL"], default=None,
                      help="Inheritance pattern: AD (dominant), AR (recessive), XL (X-linked)")
    parser.add_argument("--input", type=str, default="your_variants.csv",
                      help="Input CSV file path")
    parser.add_argument("--missense_output_csv", default="input_to_g2p_priority.csv", 
                        help="Output CSV for annotated missense variants")
    args = parser.parse_args()
    
    # Load data
    df = pd.read_csv(args.input)
    df['gnomAD4.1_joint_SAS_AF'] = pd.to_numeric(df['gnomAD4.1_joint_SAS_AF'], errors='coerce')
    # Apply classifications (if not already present)
    if 'PM1' not in df:
        df['PM1'] = df.apply(lambda row: any(start <= row['Protein_position'] <= end 
                                           for start, end in row.get('domain_ranges', [])), axis=1)
    if 'PM2' not in df:
        df['PM2'] = df['gnomAD4.1_joint_SAS_AF'].apply(lambda af: pd.isna(af) or af < 0.001)
    
    # Add prioritization columns
    df['Pathogenicity_Tier'] = df.apply(assign_pathogenicity_tier, axis=1)
    df['MOI_Compatibility'] = df.apply(
        lambda row: check_moi_compatibility(row, args.family_history), axis=1
    )
    df['Functional_Impact'] = df.apply(assign_functional_impact, axis=1)
    df['Frequency_Score'] = df.apply(assign_frequency_score, axis=1)
    
    # Sort
    df_sorted = df.sort_values(
        by=[
            'MOI_Compatibility',  # Family history match first
            'Pathogenicity_Tier',
            'Functional_Impact',
            'Frequency_Score',
            'CADD_PHRED',
            'REVEL_score'
        ],
        ascending=[False, True, True, True, False, False]
    )
    
    # Save
    output_cols = [
        'Location', 'SYMBOL', 'Consequence', 'MOI',
        'ClinVar_CLNSIG', 'CADD_PHRED', 'REVEL_score',
        'Pathogenicity_Tier', 'MOI_Compatibility'
    ]
    #output_file = f"prioritized_variants_{args.family_history or 'sporadic'}.csv"
    df_sorted.to_csv(args.missense_output_csv, index=False)
    print(f"Missense Prioritized variants saved to {args.missense_output_csv}")
    
    print(f"Prioritized {len(df)} variants | Family history: {args.family_history or 'None'}")
    #print(f"Top variant: {df_sorted.iloc[0]['SYMBOL']} (MOI={df_sorted.iloc[0]['MOI']})")
    if not df_sorted.empty:
        print(f"Top variant: {df_sorted.iloc[0]['SYMBOL']} (MOI={df_sorted.iloc[0]['MOI']})")
    else:
        print("No prioritized variants found.")

if __name__ == "__main__":
    main()
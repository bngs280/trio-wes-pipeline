import pandas as pd
import argparse

def prioritize_clinvar(clinsig):
    if pd.isna(clinsig):
        return 4  # Unknown
    clinsig = str(clinsig).lower()
    if "pathogenic" in clinsig:
        return 1
    elif "likely_pathogenic" in clinsig or "conflicting" in clinsig:
        return 2
    elif "uncertain" in clinsig or "not_provided" in clinsig:
        return 3
    elif "benign" in clinsig or "likely_benign" in clinsig:
        return 5  # Lowest
    else:
        return 4

def prioritize_af(af):
    if pd.isna(af) or af == 0:
        return 1  # Highest priority
    else:
        return 2  # Lower priority

def prioritize_moi(moi):
    """Convert MOI string to priority score (1=highest) with better handling of complex MOI strings."""
    if pd.isna(moi):
        return 4  # Unknown (lowest priority)
    moi = str(moi).strip().upper()  # Ensure consistent formatting
    if "AD" in moi or "AUTOSOMAL DOMINANT" in moi:
        return 1
    elif "X-LINKED" in moi or "XL" in moi:
        return 2
    elif "AR" in moi or "AUTOSOMAL RECESSIVE" in moi:
        return 3
    else:
        return 4  # Unknown (lowest priority)

def main(firstlist_output_csv, output_csv):
    # Load the data
    df = pd.read_csv(firstlist_output_csv)

    # Apply priority tiers
    df['ClinVar_Priority'] = df['ClinVar_CLNSIG'].apply(prioritize_clinvar)
    df['AF_Priority'] = df['gnomADe_AF'].apply(prioritize_af)
    df["MOI_Priority"] = df["MOI"].apply(prioritize_moi)

    # Sort by: MOI > AF=0 > ClinVar > HI > CADD > REVEL
    df_sorted = df.sort_values(
        by=[
            "MOI_Priority",
            "AF_Priority",
            "ClinVar_Priority",
            "Haploinsufficiency Score",
            "CADD_PHRED",
            "REVEL_score"
        ],
        ascending=[True, True, True, False, False, False]
    )

    # Filter out benign/likely benign (optional)
    df_prioritized = df_sorted[df_sorted["ClinVar_Priority"] < 5]

    # Save to CSV (include MOI and priority columns)
    output_cols = [
        "Uploaded_variation",
        "SYMBOL",
        "MOI",
        "MOI_Priority",
        "gnomADe_AF",
        "ClinVar_CLNSIG",
        "CADD_PHRED",
        "REVEL_score",
        "PVS1"
    ]
    df_prioritized.to_csv(output_csv, index=False)

    # Print top variants
    print("Top 5 Prioritized Variants:")
    print(df_prioritized[output_cols].head())

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prioritize PVS1 variants by MOI, AF, ClinVar, CADD, REVEL.")
    parser.add_argument("firstlist_output_csv", help="Input CSV file (PVS1 variants)")
    parser.add_argument("--output_csv", default="prioritized_variants.csv", help="Output CSV file (default: G2P_prioritized_variants1.csv)")
    args = parser.parse_args()
    main(args.firstlist_output_csv, args.output_csv)

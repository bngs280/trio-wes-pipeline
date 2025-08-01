import pandas as pd
import re
import argparse

# --- minimal GT parsing helpers ---
def split_alleles(gt):
    return re.split(r"[\/|]", str(gt))

def is_het(gt):
    alleles = split_alleles(gt)
    if len(alleles) != 2:
        return False
    a, b = alleles
    if '.' in (a, b):
        other = b if a == '.' else a
        return other not in ('0', '.')
    return a != b

def is_hom(gt):
    alleles = split_alleles(gt)
    if len(alleles) != 2:
        return False
    a, b = alleles
    if '.' in (a, b):
        return False
    return a == b and a != "0"

def main(input_csv, child_sex, prefix):
    # Load the data
    df_filtered = pd.read_csv(input_csv)
    
    # 1) standardize the child_sex variable
    child_sex = child_sex.upper()  # e.g. 'MALE' or 'FEMALE'
    is_female = (child_sex == 'FEMALE')
    is_male   = (child_sex == 'MALE')

    # 2) build chromosome masks using the scalar flags
    mask_auto = ~df_filtered['#CHROM'].isin(['X','Y'])
    mask_X    = df_filtered['#CHROM'] == 'X'
    mask_Y    = df_filtered['#CHROM'] == 'Y'

    mask_auto_or_femX = mask_auto | (mask_X & is_female)
    mask_X_male       = mask_X    & is_male
    mask_Y_male       = mask_Y    & is_male

    # 3) define genotype masks
    mask_child_het = df_filtered['Child_GT'].apply(is_het)
    mask_child_hom = df_filtered['Child_GT'].apply(is_hom)
    mask_child_alt = mask_child_het | mask_child_hom
    mask_child_ref = ~mask_child_alt

    mask_mom_het   = df_filtered['Mother_GT'].apply(is_het)
    mask_mom_hom   = df_filtered['Mother_GT'].apply(is_hom)
    mask_mom_alt   = mask_mom_het | mask_mom_hom
    mask_mom_ref   = ~mask_mom_alt

    mask_dad_het   = df_filtered['Father_GT'].apply(is_het)
    mask_dad_hom   = df_filtered['Father_GT'].apply(is_hom)
    mask_dad_alt   = mask_dad_het | mask_dad_hom
    mask_dad_ref   = ~mask_dad_alt

    # 4) apply filters

    # De novo
    dn_auto = df_filtered[
        mask_auto_or_femX &
        mask_child_alt &
        mask_dad_ref &
        mask_mom_ref
    ]
    dn_X_male = df_filtered[
        mask_X_male &
        mask_child_alt &
        mask_mom_ref
    ]
    dn_Y_male = df_filtered[
        mask_Y_male &
        mask_child_alt &
        mask_dad_ref
    ]
    de_novo_df = pd.concat([dn_auto, dn_X_male, dn_Y_male], ignore_index=True)

    # Autosomal recessive hom
    auto_rec_hom_df = df_filtered[
        df_filtered['#CHROM'].str.isdigit() &
        mask_child_hom &
        mask_dad_het &
        mask_mom_het
    ].reset_index(drop=True)

    # Dominant
    dom_auto = df_filtered[
        mask_auto_or_femX &
        mask_child_alt &
        (mask_dad_alt | mask_mom_alt)
    ]
    dom_X_male = df_filtered[
        mask_X_male &
        mask_child_alt &
        mask_mom_alt
    ]
    dom_Y_male = df_filtered[
        mask_Y_male &
        mask_child_alt &
        mask_dad_alt
    ]
    dominant_df = pd.concat([dom_auto, dom_X_male, dom_Y_male], ignore_index=True)

    # Paternal‐only
    pat_auto = df_filtered[
        mask_auto_or_femX &
        mask_child_alt &
        mask_dad_alt &
        mask_mom_ref
    ]
    pat_Y_male = df_filtered[
        mask_Y_male &
        mask_child_alt &
        mask_dad_alt
    ]
    pat_inherited_df = pd.concat([pat_auto, pat_Y_male], ignore_index=True)

    # Maternal‐only
    mat_auto = df_filtered[
        mask_auto_or_femX &
        mask_child_alt &
        mask_mom_alt &
        mask_dad_ref
    ]
    mat_X_male = df_filtered[
        mask_X_male &
        mask_child_alt &
        mask_mom_alt
    ]
    mat_inherited_df = pd.concat([mat_auto, mat_X_male], ignore_index=True)

    # Wild‐type
    ref_auto = df_filtered[
        mask_auto_or_femX &
        mask_child_ref &
        mask_dad_ref &
        mask_mom_ref
    ]
    ref_X_male = df_filtered[
        mask_X_male &
        mask_child_ref &
        mask_mom_ref
    ]
    ref_Y_male = df_filtered[
        mask_Y_male &
        mask_child_ref &
        mask_dad_ref
    ]
    all_ref_df = pd.concat([ref_auto, ref_X_male, ref_Y_male], ignore_index=True)

    # Print summary counts
    print("De novo:", len(de_novo_df))
    print("Autosomal recessive hom:", len(auto_rec_hom_df))
    print("Dominant:", len(dominant_df))
    print("Paternal only:", len(pat_inherited_df))
    print("Maternal only:", len(mat_inherited_df))
    print("All ref:", len(all_ref_df))

    # Optionally, output results as CSVs
    if prefix:
        de_novo_df.to_csv(f"{prefix}_de_novo.csv", index=False)
        auto_rec_hom_df.to_csv(f"{prefix}_auto_rec_hom.csv", index=False)
        dominant_df.to_csv(f"{prefix}_dominant.csv", index=False)
        pat_inherited_df.to_csv(f"{prefix}_pat_inherited.csv", index=False)
        mat_inherited_df.to_csv(f"{prefix}_mat_inherited.csv", index=False)
        all_ref_df.to_csv(f"{prefix}_all_ref.csv", index=False)
        print(f"\nFiles written with prefix: {prefix}_*.csv")

    # Concatenate all zygosity groups except "all_ref"
    all_zyg = pd.concat([
        de_novo_df,
        auto_rec_hom_df,
        dominant_df,
        pat_inherited_df,
        mat_inherited_df
    ], ignore_index=True)
    if prefix:
        all_zyg.to_csv(f"{prefix}_all_zyg.csv", index=False)
        print(f"Zygosity-union file written: {prefix}_all_zyg.csv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter trio genotypes into inheritance categories.")
    parser.add_argument("input_csv", help="Input trio genotype CSV file")
    parser.add_argument("--child_sex", default="FEMALE", help="Child sex (MALE or FEMALE, default: FEMALE)")
    parser.add_argument("--prefix", default=None, help="Output file prefix (if set, writes category CSVs)")
    args = parser.parse_args()
    main(args.input_csv, args.child_sex, args.prefix)

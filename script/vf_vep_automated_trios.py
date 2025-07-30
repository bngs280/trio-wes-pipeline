# ### worked
# ## Usage:  python3 script2.py --vep cohort.phased.merged.annotatedvep.txt --clingen ClinGen_gene_curation_list_GRCh38_cleaned.tsv --ped trioPRJ8766.ped --vcf cohort.phased.merged.vcf.gz --output annotated_trio050725.tsv
# import argparse
# import pandas as pd
# import numpy as np
# import gzip
# import sys
# import re

# from mendelian_inheritance3 import mendelian_tag

# # ------------------ Block 1: VEP File ------------------

# def load_vep(vep_file):
#     headers = []
#     with open(vep_file, "r") as f:
#         for line in f:
#             if line.startswith("##"):
#                 headers.append(line)
#             elif line.startswith("#"):
#                 headers.append(line)
#                 break
#     meta = headers[:-1]
#     df = pd.read_csv(vep_file, sep="\t", skiprows=len(meta), low_memory=False)
#     df.columns = [c.lstrip("#") for c in df.columns]
#     exclude = {"Location", "Uploaded_variation", "Allele", "UPLOADED_ALLELE"}
#     cols = [c for c in df.columns if c not in exclude]
#     df.loc[:, cols] = df.loc[:, cols].replace("-", np.nan)

# #    def adjust_location(r):
# #        try:
# #            allele = str(r.get("Allele","")).strip()
# #            loc    = str(r["Location"]).strip()
# #            uv     = str(r["Uploaded_variation"]).strip()
# #            chrom, pos = loc.split(":")
# #            left = int(pos.split("-")[0])
# #            up   = int(uv.split("_")[1])
# #            if len(allele)==1 and allele.upper() in "ACGT":
# #                return loc if "-" not in pos else f"{chrom}:{pos.split('-')[0]}"
# #            return f"{chrom}:{left-1}" if left != up-1 else f"{chrom}:{left}"
# #        except:
# #            return r.get("Location","")
# #    df["Adjusted_Location"] = df.apply(adjust_location, axis=1)
#     def adjust_location(r):
#         allele = str(r.get("Allele", "")).strip().upper()
#         loc = str(r["Location"]).strip()
#         uv = str(r["Uploaded_variation"]).strip()  # e.g., "chr1_1005757_A_AT"
    
#         # Extract VCF's POS from Uploaded_variation (authoritative)
#         chrom_uv, pos_uv, ref_uv, alt_uv = uv.split("_")[:4]
#         pos_uv = int(pos_uv)  # VCF's 1-based POS
    
#         # Case 1: Deletions (use VCF's POS)
#         if "DEL" in allele or (len(ref_uv) > len(alt_uv)):
#             return f"{chrom_uv}:{pos_uv}"
    
#         # Case 2: Insertions (use VCF's POS)
#         elif len(ref_uv) < len(alt_uv):  # Insertion
#             return f"{chrom_uv}:{pos_uv}"
    
#         # Case 3: SNPs (trim ranges if needed)
#         elif len(allele) == 1 and allele in "ACGT":
#             return loc.split("-")[0] if "-" in loc else loc
    
#         # Default: Keep unchanged
#         else:
#             return loc

#     df["Adjusted_Location"] = df.apply(adjust_location, axis=1)
#     return df

# # ------------------ Block 2: ClinGen File ------------------

# def load_clingen(path, df):
#     df_clin = pd.read_csv(path, sep="\t", low_memory=False, dtype= str)
#     haplo_pm = [c for c in df.columns if c.startswith("Haploinsufficiency PMID")]
#     triplo_pm= [c for c in df.columns if c.startswith("Triplosensitivity PMID")]
#     def collapse_pmids(row, cols):
#         pmids = []
#         for col in cols:
#             val = row.get(col)
#             if pd.notna(val):
#                 try:
#                     pmids.append(str(int(float(val))))
#                 except (ValueError, TypeError):
#                     pass
#         return ';'.join(sorted(set(pmids))) if pmids else np.nan
#     df_clin['Haploinsufficiency_PMIDs'] = df_clin.apply(lambda r: collapse_pmids(r, haplo_pm), axis=1)
#     df_clin['Triplosensitivity_PMIDs']  = df_clin.apply(lambda r: collapse_pmids(r, triplo_pm), axis=1)
#     def combine_series_semis(series):
#         items = set()
#         for cell in series.dropna():
#             items.update(cell.split(';'))
#         return ';'.join(sorted(items)) if items else np.nan
#     df_clin_agg = df_clin.groupby('Gene Symbol', as_index=False).agg({
#         'Haploinsufficiency Score':     lambda s: pd.to_numeric(s, errors='coerce').max(),
#         'Triplosensitivity Score':      lambda s: pd.to_numeric(s, errors='coerce').max(),
#         'Haploinsufficiency Disease ID':lambda s: ';'.join(sorted(set(s.dropna()))),
#         'Triplosensitivity Disease ID': lambda s: ';'.join(sorted(set(s.dropna()))),
#         'Genomic Location':             lambda s: ';'.join(sorted(set(s.dropna()))),
#         'Haploinsufficiency_PMIDs':     lambda s: combine_series_semis(s),
#         'Triplosensitivity_PMIDs':      lambda s: combine_series_semis(s),
#     })
#     return df_clin_agg

# # ------------------ Block 3: PED File ------------------

# def load_ped_file(ped_path):
#     ped = pd.read_csv(
#         ped_path,
#         sep=r'\s+',
#         engine='python',
#         comment='#',
#         header=None,
#         usecols=[1, 2, 3, 4],  # IID, Father_ID, Mother_ID, Sex
#         names=['IID', 'Father', 'Mother', 'Sex']
#     )
#     parents = set(ped['Father']).union(set(ped['Mother']))
#     child_ids = [iid for iid in ped['IID'] if iid not in parents]
#     if len(child_ids) != 1:
#         raise ValueError(f"Expected exactly one proband, found {len(child_ids)}")
#     child_id = child_ids[0]
#     mother_id = ped.loc[ped['IID'] == child_id, 'Mother'].iloc[0]
#     father_id = ped.loc[ped['IID'] == child_id, 'Father'].iloc[0]
#     sex_map = {}
#     for iid, sex in zip(ped['IID'], ped['Sex']):
#         if iid == child_id:
#             role = 'child'
#         elif iid == mother_id:
#             role = 'mother'
#         elif iid == father_id:
#             role = 'father'
#         else:
#             continue
#         if sex in (1, '1'):
#             sex_map[(role, iid)] = 'male'
#         elif sex in (2, '2'):
#             sex_map[(role, iid)] = 'female'
#         else:
#             sex_map[(role, iid)] = None
#     if sex_map.get(('mother', mother_id)) == 'male':
#         raise ValueError(f"Mother {mother_id} is marked as male")
#     if sex_map.get(('father', father_id)) == 'female':
#         raise ValueError(f"Father {father_id} is marked as female")
#     return {
#         'child': child_id,
#         'mother': mother_id,
#         'father': father_id,
#         'sex_map': sex_map
#     }

# # ------------------ Block 4: VCF File ------------------

# def read_vcf(vcf):
#     vcf_headers = []
#     vcf_lines   = []
#     with gzip.open(vcf, "rt") as f:
#         for line in f:
#             if line.startswith("##"):
#                 vcf_headers.append(line)
#             elif line.startswith("#CHROM"):
#                 vcf_headers.append(line)
#                 vcf_columns = line.rstrip().split("\t")
#             else:
#                 vcf_lines.append(line.rstrip("\n").split("\t"))
#     fixed_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
#     sample_cols = vcf_columns[9:]
#     df_vcf = pd.DataFrame(vcf_lines, columns=fixed_columns + sample_cols)
#     required_cols = ['#CHROM', 'POS', 'REF', 'ALT', 'FORMAT']
#     if not all(col in df_vcf.columns for col in required_cols):
#         print("Error: VCF missing required columns.", file=sys.stderr)
#         sys.exit(1)
#     df_vcf["CHROM_POS"] = df_vcf["#CHROM"] + ":" + df_vcf["POS"].astype(str)
#     format_fields = df_vcf["FORMAT"].str.split(":", expand=False)
#     sample_dfs = []
#     for sample in sample_cols:
#         geno_splits = df_vcf[sample].str.split(":", expand=False)
#         keys = [k for sub in format_fields for k in sub]
#         vals = [v for sub in geno_splits for v in sub]
#         lengths = format_fields.str.len().to_numpy()
#         rows_i = np.repeat(df_vcf.index.values, lengths)
#         long = pd.DataFrame({"row": rows_i, "key": keys, "val": vals})
#         sample_df = (
#             long.pivot(index="row", columns="key", values="val")
#             .reindex(index=df_vcf.index)
#             .fillna(".")
#         )
#         sample_df.columns = [(sample, col) for col in sample_df.columns]
#         sample_dfs.append(sample_df)
#     sample_combined = pd.concat(sample_dfs, axis=1)
#     vcf = pd.concat([df_vcf.drop(columns=["FORMAT"] + sample_cols), sample_combined], axis=1)
#     samples = sample_cols
#     for s in samples:
#         if ('DP', s) not in vcf.columns:
#             if ('DPI', s) in vcf.columns:
#                 vcf[(s, 'DP')] = vcf[(s, 'DPI')]
#             elif ('AD', s) in vcf.columns:
#                 a0a1 = vcf[(s, 'AD')].str.split(',', expand=True).astype(float)
#                 vcf[(s, 'DP')] = a0a1.sum(axis=1)
#         if ('DPI', s) in vcf.columns:
#             vcf.drop(columns=[('DPI', s)], inplace=True)
#         if ('GQ', s) not in vcf.columns:
#             def calc_gq(r, sample):
#                 pl = None
#                 if ('PL', sample) in r and pd.notna(r[('PL', sample)]):
#                     raw = list(map(int, str(r[('PL', sample)]).split(",")))
#                     base = min(raw)
#                     pl = [p - base for p in raw]
#                 elif ('GL', sample) in r and pd.notna(r[('GL', sample)]):
#                     raw = [round(-10 * float(x)) for x in str(r[('GL', sample)]).split(",")]
#                     base = min(raw)
#                     pl = [p - base for p in raw]
#                 if pl is None:
#                     return np.nan
#                 gt = str(r[('GT', sample)]).replace("|", "/")
#                 idx = {'0/0': 0, '0/1': 1, '1/1': 2}.get(gt)
#                 if idx is None:
#                     return np.nan
#                 return min(p for i, p in enumerate(pl) if i != idx)
#             vcf[(s, 'GQ')] = vcf.apply(lambda r: calc_gq(r, s), axis=1)
#         if ('PL', s) not in vcf.columns and ('GL', s) in vcf.columns:
#             def gl_to_pl(gl):
#                 raw = [round(-10 * float(x)) for x in str(gl).split(",")]
#                 base = min(raw)
#                 return ",".join(str(p - base) for p in raw)
#             vcf[(s, 'PL')] = vcf[(s, 'GL')].apply(gl_to_pl)
#         if (s,'VAF') not in vcf.columns and (s,'AD') in vcf.columns and (s,'DP') in vcf.columns:
#             def calc_vaf(ad, dp):
#                 try:
#                     a0, a1 = map(int, ad.split(","))
#                     return a1 / (a0 + a1) if (a0 + a1) > 0 else np.nan
#                 except:
#                     return np.nan
#             vcf[(s, 'VAF')] = vcf.apply(lambda r: calc_vaf(r[(s, 'AD')], r[(s, 'DP')]), axis=1)
#         for field in ('DP', 'GQ', 'VAF'):
#             if (field, s) in vcf.columns:
#                 vcf[(field, s)] = pd.to_numeric(vcf[(field, s)], errors="coerce")
#     return vcf

# # ------------------ Block 5: Mendelian Inheritance ------------------

# def add_mendelian_inheritance(vcf, ped_data, gq_thresh=20, dp_thresh=10, vaf_mosaic=(0.05, 0.35)):
#     child_id = ped_data['child']
#     mother_id = ped_data['mother']
#     father_id = ped_data['father']
#     sex_map = ped_data['sex_map']
#     child_sex = sex_map[('child', child_id)]
#     mother_sex = sex_map.get(('mother', mother_id), 'female')
#     father_sex = sex_map.get(('father', father_id), 'male')
#     vcf['INHERITANCE'] = vcf.apply(
#         lambda row: mendelian_tag(
#             chrom=row['#CHROM'].replace('chr', '').upper(),
#             child_gt=row[(child_id, 'GT')],
#             mother_gt=row[(mother_id, 'GT')],
#             father_gt=row[(father_id, 'GT')],
#             child_sex=child_sex,
#             mother_sex=mother_sex,
#             father_sex=father_sex,
#             child_gq   =  pd.to_numeric(row.get((child_id,  'GQ')),   errors='coerce'),
#             mother_gq  =  pd.to_numeric(row.get((mother_id, 'GQ')),   errors='coerce'),
#             father_gq  =  pd.to_numeric(row.get((father_id, 'GQ')),   errors='coerce'),
#             child_dp   =  pd.to_numeric(row.get((child_id,  'DP')),   errors='coerce'),
#             mother_dp  =  pd.to_numeric(row.get((mother_id, 'DP')),   errors='coerce'),
#             father_dp  =  pd.to_numeric(row.get((father_id, 'DP')),   errors='coerce'),
#             child_vaf =  pd.to_numeric(row.get((child_id, 'VAF')),   errors='coerce')
#         ),
#         axis=1
#     )
#     vcf['label'] = vcf['INHERITANCE']
#     return vcf

# def simplify_columns(vcf, ped_data):
#     column_map = {}
#     role_names = {'mother': 'Mother', 'father': 'Father', 'child': 'Child'}
#     relationships = {'child': ped_data['child'], 'mother': ped_data['mother'], 'father': ped_data['father']}
#     id_to_role = {v: k for k, v in relationships.items()}
#     for col in vcf.columns:
#         if isinstance(col, str):
#             column_map[col] = col
#         else:
#             sample_id, field = col
#             role = id_to_role.get(sample_id)
#             if role:
#                 new_name = f"{role_names[role]}_{field.upper()}"
#                 column_map[col] = new_name
#             else:
#                 column_map[col] = col
#     return vcf.rename(columns=column_map)

# # ------------------ Block 6: Merge VEP + ClinGen ------------------

# def merge_vep_clingen(df, df_clin_agg):
#     clingen_with_vcf = df.merge(
#         df_clin_agg,
#         left_on="SYMBOL",
#         right_on="Gene Symbol",
#         how="left",
#         suffixes=('_vcf','_clingen'),
#         validate='many_to_one'
#     ).drop(columns=['Gene ID','cytoBand','Genomic Location','Date Last Evaluated'], errors='ignore') \
#      .drop_duplicates()
#     return clingen_with_vcf

# # ------------------ Block 7: Merge with VCF ------------------

# def merge_with_vcf(clingen_with_vcf, vcf_simple):
#     merged_df = pd.merge(
#         clingen_with_vcf,
#         vcf_simple,
#         left_on='Adjusted_Location',
#         right_on='CHROM_POS',
#         how='left'
#     )
#     merged_df['ClinVar_CLNDISDB'] = merged_df['ClinVar_CLNDISDB'].replace(r'^\s*$', np.nan, regex=True)
#     merged_df['ClinVar_CLNDISDB'] = merged_df['ClinVar_CLNDISDB'].str.replace(
#         r'\b(MONDO):\1:', r'\1:', regex=True)
#     return merged_df

# # ------------------ Block 8: Compound Hets ------------------

# def add_compound_hets(merged_df, ped_data, gene_col='Gene'):
#     child_id  = ped_data['child']
#     mother_id = ped_data['mother']
#     father_id = ped_data['father']
#     if   'Child_GT' in merged_df.columns:
#         gt_col = 'Child_GT'
#     elif (child_id, 'GT') in merged_df.columns:
#         gt_col = (child_id, 'GT')
#     else:
#         raise KeyError(f"Could not find GT column for child {child_id}")
#     if   'Mother_GT' in merged_df.columns:
#         mom_gt_col = 'Mother_GT'
#     elif (mother_id, 'GT') in merged_df.columns:
#         mom_gt_col = (mother_id, 'GT')
#     else:
#         raise KeyError(f"Could not find GT column for mother {mother_id}")
#     if   'Father_GT' in merged_df.columns:
#         dad_gt_col = 'Father_GT'
#     elif (father_id, 'GT') in merged_df.columns:
#         dad_gt_col = (father_id, 'GT')
#     else:
#         raise KeyError(f"Could not find GT column for father {father_id}")
#     mat = merged_df[merged_df['INHERITANCE']=='ad_maternal']
#     pat = merged_df[merged_df['INHERITANCE']=='ad_paternal']
#     for gene, mat_sub in mat.groupby(gene_col):
#         pat_sub = pat[pat[gene_col]==gene]
#         if pat_sub.empty:
#             continue
#         for mi in mat_sub.index:
#             for pi in pat_sub.index:
#                 in_trans = True
#                 gt_m   = merged_df.at[mi, gt_col]
#                 gt_p   = merged_df.at[pi, gt_col]
#                 mom_gt = merged_df.at[mi, mom_gt_col]
#                 dad_gt = merged_df.at[pi, dad_gt_col]
#                 if mom_gt == './.' or dad_gt == './.':
#                     continue
#                 if ('|' in gt_m and '|' in gt_p and '|' in mom_gt and '|' in dad_gt):
#                     m0,  m1  = gt_m.split('|')
#                     mm0, mm1 = mom_gt.split('|')
#                     p0,  p1  = gt_p.split('|')
#                     fd0, fd1 = dad_gt.split('|')
#                     mat_from_mom = (m1  == '1') and (mm1 == '1')
#                     pat_from_dad = (p0  == '1') and (fd0 == '1')
#                     in_trans = mat_from_mom and pat_from_dad
#                 if in_trans:
#                     merged_df.at[mi, 'INHERITANCE'] = 'ar_comp_het'
#                     merged_df.at[pi, 'INHERITANCE'] = 'ar_comp_het'
#     return merged_df

# # ------------------ Block 9: Zygosity Classification ------------------

# def classify_zygosity(gt, vaf=None, chrom=None, pos=None, sex=None):
#     chrom_clean = None
#     if isinstance(chrom, str):
#         chrom_clean = chrom.lstrip("chr").upper()
#     try:
#         pos = int(pos)
#     except (TypeError, ValueError):
#         pos = None
#     PAR_RANGES = {
#         'PAR1': {'X': (10000, 2781479),      'Y': (10000, 2781479)},
#         'PAR2': {'X': (155701382, 156030895),'Y': (56887903, 57217415)},
#     }
#     in_par = False
#     if chrom_clean in {'X', 'Y'} and pos is not None:
#         for start, end in PAR_RANGES['PAR1'][chrom_clean], PAR_RANGES['PAR2'][chrom_clean]:
#             if start <= pos <= end:
#                 in_par = True
#                 break
#     if pd.isna(gt):
#         alleles = []
#     else:
#         alleles = re.split(r"[\/|]", str(gt))
#         alleles = [a for a in alleles if a != ""]
#     if (isinstance(sex, str) and sex.lower() == "male"
#         and chrom_clean in {"X", "Y"} and not in_par):
#         if len(alleles) == 1 and alleles[0] != ".":
#             return "Hemizygous"
#     if len(alleles) == 2 and "." not in alleles:
#         return "Homozygous" if alleles[0] == alleles[1] else "Heterozygous"
#     if vaf is not None and not pd.isna(vaf):
#         try:
#             v = float(vaf)
#         except (TypeError, ValueError):
#             return "Unclassified"
#         if v >= 0.85 or v <= 0.15:
#             return "Homozygous"
#         if 0.35 <= v <= 0.65:
#             return "Heterozygous"
#     return "Unclassified"

# def add_zygosity(merged_df):
#     for col in ["Child_DP", "Child_GQ", "Child_VAF"]:
#         merged_df[col] = pd.to_numeric(merged_df[col], errors="coerce")
#     merged_df["POS"] = pd.to_numeric(merged_df["POS"], errors="coerce")
#     # You may need to add a column for Child_sex from ped_data
#     merged_df["Zygosity"] = merged_df.apply(
#         lambda r: classify_zygosity(
#             gt    = r["Child_GT"],
#             vaf   = r.get("Child_VAF", None),
#             chrom = r.get("#CHROM", r.get("CHROM", None)),
#             pos   = r.get("POS", None),
#             sex   = r.get("Child_sex", None)
#         ),
#         axis=1
#     )
#     return merged_df

# # ------------------ Block 10: Final Labels and Output ------------------

# def finalize_labels(merged_df, zmap):
#     merged_df.insert(
#         merged_df.columns.get_loc("Zygosity")+1,
#         "Zygosity_label",
#         merged_df["Zygosity"].map(zmap)
#     )
#     clinvar_map = {
#         "likely_benign": "LB",
#         "benign": "Ben",
#         "likely_pathogenic": "LP",
#         "pathogenic": "PAT",
#         "uncertain_significance": "VUS",
#         "conflicting_interpretations_of_pathogenicity": "Conflicting"
#     }
#     def map_clinsig(c):
#         if pd.isna(c): return "NA"
#         pat = r'([^/,]+)'
#         def rep(m):
#             term = m.group(0).strip().lower()
#             return clinvar_map.get(term, term)
#         return re.sub(pat, rep, c)
#     merged_df.insert(
#         merged_df.columns.get_loc('ClinVar_CLNSIG')+1,
#         "CLINSIG_label",
#         merged_df['ClinVar_CLNSIG'].apply(map_clinsig)
#     )
#     def merge_prioritized(r):
#         if pd.notna(r.get('OMIM', None)):
#             return f"OMIM:{r.get('OMIM')}"
#         elif pd.notna(r.get('Orpha', None)):
#             return f"Orpha:{r.get('Orpha')}"
#         elif pd.notna(r.get('MONDO', None)):
#             return f"MONDO:{r.get('MONDO')}"
#         elif pd.notna(r.get('MedGen', None)):
#             return f"MedGen:{r.get('MedGen')}"
#         return ""
#     # def merge_prioritized(r):
#     #     if pd.notna(r['OMIM']):    return f"OMIM:{r['OMIM']}"
#     #     elif pd.notna(r['Orpha']): return f"Orpha:{r['Orpha']}"
#     #     elif pd.notna(r['MONDO']): return f"Orpha:{r['MONDO']}"
#     #     elif pd.notna(r['MedGen']):    return f"MedGen:{r['MedGen']}"
#     #     return ""
#     merged_df['ACMG DISEASE ID (OMIM/ORPHA ID)'] = merged_df.apply(merge_prioritized, axis=1)
# #    merged_df = merged_df.rename(columns={
# #        'Gene': 'GENE (GENE ID)',
# #        'SYMBOL': 'SYMBOL (Gene Name)',
# #        'Existing_variation': 'SNPs/Rsid',
# #        'ClinVar_CLNSIG': 'CLINVAR CLNSIG',
# #        'Adjusted_Location': 'LOCATION',
# #        'Consequence': 'EFFECT'
# #    })
#     return merged_df

# # ------------------ Main Automation Script ------------------

# def main():
#     parser = argparse.ArgumentParser(description="Automated trio variant annotation pipeline")
#     parser.add_argument('--vep', required=True, help='VEP-annotated file')
#     parser.add_argument('--clingen', required=True, help='ClinGen gene curation file')
#     parser.add_argument('--ped', required=True, help='PED file')
#     parser.add_argument('--vcf', required=True, help='Phased VCF file (gzipped)')
#     parser.add_argument('--output', required=True, help='Output file (TSV)')
#     args = parser.parse_args()

#     print("Loading VEP...")
#     df = load_vep(args.vep)
#     print("Loading ClinGen...")
#     df_clin_agg = load_clingen(args.clingen, df)
#     print("Loading PED...")
#     ped_data = load_ped_file(args.ped)
#     child_sex = ped_data['sex_map'][('child', ped_data['child'])]
#     with open("child_sex.txt", "w") as f:
#         f.write(child_sex)
#     print("Loading VCF...")
#     vcf = read_vcf(args.vcf)
#     print("Tagging Mendelian inheritance...")
#     vcf = add_mendelian_inheritance(vcf, ped_data)
#     vcf_simple = simplify_columns(vcf, ped_data)
#     print("Merging VEP and ClinGen...")
#     clingen_with_vcf = merge_vep_clingen(df, df_clin_agg)
#     print("Merging with VCF...")
#     merged_df = merge_with_vcf(clingen_with_vcf, vcf_simple)
#     print("Flagging compound hets...")
#     merged_df = add_compound_hets(merged_df, ped_data)
#     print("Computing zygosity...")
#     merged_df = add_zygosity(merged_df)
#     # You must define zmap for Zygosity_label mapping
#     zmap = {"Homozygous": "Hom", "Heterozygous": "Het", "Hemizygous": "Hem", "Unclassified": "Unc"}
#     print("Finalizing labels...")
#     merged_df = finalize_labels(merged_df, zmap)
#     print(f"Saving output to {args.output}")
#     merged_df.to_csv(args.output, sep="\t", index=False)
#     print("Done.")

# if __name__ == "__main__":
#     main()


################ ACMG integration ##############
### worked
## Usage:  python3 script2.py --vep cohort.phased.merged.annotatedvep.txt --clingen ClinGen_gene_curation_list_GRCh38_cleaned.tsv --ped trioPRJ8766.ped --vcf cohort.phased.merged.vcf.gz --output annotated_trio050725.tsv
import argparse
import pandas as pd
import numpy as np
import gzip
import sys
import re

from mendelian_inheritance3 import mendelian_tag

# ------------------ Block 1: VEP File ------------------

def load_vep(vep_file):
    headers = []
    with open(vep_file, "r") as f:
        for line in f:
            if line.startswith("##"):
                headers.append(line)
            elif line.startswith("#"):
                headers.append(line)
                break
    meta = headers[:-1]
    df = pd.read_csv(vep_file, sep="\t", skiprows=len(meta), low_memory=False)
    df.columns = [c.lstrip("#") for c in df.columns]
    exclude = {"Location", "Uploaded_variation", "Allele", "UPLOADED_ALLELE"}
    cols = [c for c in df.columns if c not in exclude]
    df.loc[:, cols] = df.loc[:, cols].replace("-", np.nan)

#    def adjust_location(r):
#        try:
#            allele = str(r.get("Allele","")).strip()
#            loc    = str(r["Location"]).strip()
#            uv     = str(r["Uploaded_variation"]).strip()
#            chrom, pos = loc.split(":")
#            left = int(pos.split("-")[0])
#            up   = int(uv.split("_")[1])
#            if len(allele)==1 and allele.upper() in "ACGT":
#                return loc if "-" not in pos else f"{chrom}:{pos.split('-')[0]}"
#            return f"{chrom}:{left-1}" if left != up-1 else f"{chrom}:{left}"
#        except:
#            return r.get("Location","")
#    df["Adjusted_Location"] = df.apply(adjust_location, axis=1)
    def adjust_location(r):
        allele = str(r.get("Allele", "")).strip().upper()
        loc = str(r["Location"]).strip()
        uv = str(r["Uploaded_variation"]).strip()  # e.g., "chr1_1005757_A_AT"
    
        # Extract VCF's POS from Uploaded_variation (authoritative)
        chrom_uv, pos_uv, ref_uv, alt_uv = uv.split("_")[:4]
        pos_uv = int(pos_uv)  # VCF's 1-based POS
    
        # Case 1: Deletions (use VCF's POS)
        if "DEL" in allele or (len(ref_uv) > len(alt_uv)):
            return f"{chrom_uv}:{pos_uv}"
    
        # Case 2: Insertions (use VCF's POS)
        elif len(ref_uv) < len(alt_uv):  # Insertion
            return f"{chrom_uv}:{pos_uv}"
    
        # Case 3: SNPs (trim ranges if needed)
        elif len(allele) == 1 and allele in "ACGT":
            return loc.split("-")[0] if "-" in loc else loc
    
        # Default: Keep unchanged
        else:
            return loc

    df["Adjusted_Location"] = df.apply(adjust_location, axis=1)
    return df


# ------------------ Block 2: ClinGen File ------------------

def load_clingen(path, df):
    df_clin = pd.read_csv(path, sep="\t", low_memory=False, dtype= str)
    haplo_pm = [c for c in df.columns if c.startswith("Haploinsufficiency PMID")]
    triplo_pm= [c for c in df.columns if c.startswith("Triplosensitivity PMID")]
    def collapse_pmids(row, cols):
        pmids = []
        for col in cols:
            val = row.get(col)
            if pd.notna(val):
                try:
                    pmids.append(str(int(float(val))))
                except (ValueError, TypeError):
                    pass
        return ';'.join(sorted(set(pmids))) if pmids else np.nan
    df_clin['Haploinsufficiency_PMIDs'] = df_clin.apply(lambda r: collapse_pmids(r, haplo_pm), axis=1)
    df_clin['Triplosensitivity_PMIDs']  = df_clin.apply(lambda r: collapse_pmids(r, triplo_pm), axis=1)
    def combine_series_semis(series):
        items = set()
        for cell in series.dropna():
            items.update(cell.split(';'))
        return ';'.join(sorted(items)) if items else np.nan
    df_clin_agg = df_clin.groupby('Gene Symbol', as_index=False).agg({
        'Haploinsufficiency Score':     lambda s: pd.to_numeric(s, errors='coerce').max(),
        'Triplosensitivity Score':      lambda s: pd.to_numeric(s, errors='coerce').max(),
        'Haploinsufficiency Disease ID':lambda s: ';'.join(sorted(set(s.dropna()))),
        'Triplosensitivity Disease ID': lambda s: ';'.join(sorted(set(s.dropna()))),
        'Genomic Location':             lambda s: ';'.join(sorted(set(s.dropna()))),
        'Haploinsufficiency_PMIDs':     lambda s: combine_series_semis(s),
        'Triplosensitivity_PMIDs':      lambda s: combine_series_semis(s),
    })
    return df_clin_agg

# ------------------ Block 3: PED File ------------------

def load_ped_file(ped_path):
    ped = pd.read_csv(
        ped_path,
        sep=r'\s+',
        engine='python',
        comment='#',
        header=None,
        usecols=[1, 2, 3, 4],  # IID, Father_ID, Mother_ID, Sex
        names=['IID', 'Father', 'Mother', 'Sex']
    )
    parents = set(ped['Father']).union(set(ped['Mother']))
    child_ids = [iid for iid in ped['IID'] if iid not in parents]
    if len(child_ids) != 1:
        raise ValueError(f"Expected exactly one proband, found {len(child_ids)}")
    child_id = child_ids[0]
    mother_id = ped.loc[ped['IID'] == child_id, 'Mother'].iloc[0]
    father_id = ped.loc[ped['IID'] == child_id, 'Father'].iloc[0]
    sex_map = {}
    for iid, sex in zip(ped['IID'], ped['Sex']):
        if iid == child_id:
            role = 'child'
        elif iid == mother_id:
            role = 'mother'
        elif iid == father_id:
            role = 'father'
        else:
            continue
        if sex in (1, '1'):
            sex_map[(role, iid)] = 'male'
        elif sex in (2, '2'):
            sex_map[(role, iid)] = 'female'
        else:
            sex_map[(role, iid)] = None
    if sex_map.get(('mother', mother_id)) == 'male':
        raise ValueError(f"Mother {mother_id} is marked as male")
    if sex_map.get(('father', father_id)) == 'female':
        raise ValueError(f"Father {father_id} is marked as female")
    return {
        'child': child_id,
        'mother': mother_id,
        'father': father_id,
        'sex_map': sex_map
    }

# ------------------ Block 4: VCF File ------------------

def read_vcf(vcf):
    vcf_headers = []
    vcf_lines   = []
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("##"):
                vcf_headers.append(line)
            elif line.startswith("#CHROM"):
                vcf_headers.append(line)
                vcf_columns = line.rstrip().split("\t")
            else:
                vcf_lines.append(line.rstrip("\n").split("\t"))
    fixed_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    sample_cols = vcf_columns[9:]
    df_vcf = pd.DataFrame(vcf_lines, columns=fixed_columns + sample_cols)
    required_cols = ['#CHROM', 'POS', 'REF', 'ALT', 'FORMAT']
    if not all(col in df_vcf.columns for col in required_cols):
        print("Error: VCF missing required columns.", file=sys.stderr)
        sys.exit(1)
    df_vcf["CHROM_POS"] = df_vcf["#CHROM"] + ":" + df_vcf["POS"].astype(str)
    format_fields = df_vcf["FORMAT"].str.split(":", expand=False)
    sample_dfs = []
    for sample in sample_cols:
        geno_splits = df_vcf[sample].str.split(":", expand=False)
        keys = [k for sub in format_fields for k in sub]
        vals = [v for sub in geno_splits for v in sub]
        lengths = format_fields.str.len().to_numpy()
        rows_i = np.repeat(df_vcf.index.values, lengths)
        long = pd.DataFrame({"row": rows_i, "key": keys, "val": vals})
        sample_df = (
            long.pivot(index="row", columns="key", values="val")
            .reindex(index=df_vcf.index)
            .fillna(".")
        )
        sample_df.columns = [(sample, col) for col in sample_df.columns]
        sample_dfs.append(sample_df)
    sample_combined = pd.concat(sample_dfs, axis=1)
    vcf = pd.concat([df_vcf.drop(columns=["FORMAT"] + sample_cols), sample_combined], axis=1)
    samples = sample_cols
    for s in samples:
        if ('DP', s) not in vcf.columns:
            if ('DPI', s) in vcf.columns:
                vcf[(s, 'DP')] = vcf[(s, 'DPI')]
            elif ('AD', s) in vcf.columns:
                a0a1 = vcf[(s, 'AD')].str.split(',', expand=True).astype(float)
                vcf[(s, 'DP')] = a0a1.sum(axis=1)
        if ('DPI', s) in vcf.columns:
            vcf.drop(columns=[('DPI', s)], inplace=True)
        if ('GQ', s) not in vcf.columns:
            def calc_gq(r, sample):
                pl = None
                if ('PL', sample) in r and pd.notna(r[('PL', sample)]):
                    raw = list(map(int, str(r[('PL', sample)]).split(",")))
                    base = min(raw)
                    pl = [p - base for p in raw]
                elif ('GL', sample) in r and pd.notna(r[('GL', sample)]):
                    raw = [round(-10 * float(x)) for x in str(r[('GL', sample)]).split(",")]
                    base = min(raw)
                    pl = [p - base for p in raw]
                if pl is None:
                    return np.nan
                gt = str(r[('GT', sample)]).replace("|", "/")
                idx = {'0/0': 0, '0/1': 1, '1/1': 2}.get(gt)
                if idx is None:
                    return np.nan
                return min(p for i, p in enumerate(pl) if i != idx)
            vcf[(s, 'GQ')] = vcf.apply(lambda r: calc_gq(r, s), axis=1)
        if ('PL', s) not in vcf.columns and ('GL', s) in vcf.columns:
            def gl_to_pl(gl):
                raw = [round(-10 * float(x)) for x in str(gl).split(",")]
                base = min(raw)
                return ",".join(str(p - base) for p in raw)
            vcf[(s, 'PL')] = vcf[(s, 'GL')].apply(gl_to_pl)
        if (s,'VAF') not in vcf.columns and (s,'AD') in vcf.columns and (s,'DP') in vcf.columns:
            def calc_vaf(ad, dp):
                try:
                    a0, a1 = map(int, ad.split(","))
                    return a1 / (a0 + a1) if (a0 + a1) > 0 else np.nan
                except:
                    return np.nan
            vcf[(s, 'VAF')] = vcf.apply(lambda r: calc_vaf(r[(s, 'AD')], r[(s, 'DP')]), axis=1)
        for field in ('DP', 'GQ', 'VAF'):
            if (field, s) in vcf.columns:
                vcf[(field, s)] = pd.to_numeric(vcf[(field, s)], errors="coerce")
    return vcf

# ------------------ Block 5: Mendelian Inheritance ------------------

def add_mendelian_inheritance(vcf, ped_data, gq_thresh=20, dp_thresh=10, vaf_mosaic=(0.05, 0.35)):
    child_id = ped_data['child']
    mother_id = ped_data['mother']
    father_id = ped_data['father']
    sex_map = ped_data['sex_map']
    child_sex = sex_map[('child', child_id)]
    mother_sex = sex_map.get(('mother', mother_id), 'female')
    father_sex = sex_map.get(('father', father_id), 'male')
    vcf['INHERITANCE'] = vcf.apply(
        lambda row: mendelian_tag(
            chrom=row['#CHROM'].replace('chr', '').upper(),
            child_gt=row[(child_id, 'GT')],
            mother_gt=row[(mother_id, 'GT')],
            father_gt=row[(father_id, 'GT')],
            child_sex=child_sex,
            mother_sex=mother_sex,
            father_sex=father_sex,
            child_gq   =  pd.to_numeric(row.get((child_id,  'GQ')),   errors='coerce'),
            mother_gq  =  pd.to_numeric(row.get((mother_id, 'GQ')),   errors='coerce'),
            father_gq  =  pd.to_numeric(row.get((father_id, 'GQ')),   errors='coerce'),
            child_dp   =  pd.to_numeric(row.get((child_id,  'DP')),   errors='coerce'),
            mother_dp  =  pd.to_numeric(row.get((mother_id, 'DP')),   errors='coerce'),
            father_dp  =  pd.to_numeric(row.get((father_id, 'DP')),   errors='coerce'),
            child_vaf =  pd.to_numeric(row.get((child_id, 'VAF')),   errors='coerce')
        ),
        axis=1
    )
    vcf['label'] = vcf['INHERITANCE']
    return vcf

# ------------------ Block : ACMG ------------------
def read_acmg(acmg_file):
    df = pd.read_csv(acmg_file, sep="\t", low_memory=False, dtype=str)
    df = df[['chromosome', 'position', 'refAllele', 'altAllele', 'acmgClassification', 'rationale']]
    df = df.rename(columns={'chromosome': '#CHROM', 'position': 'POS', 'refAllele': 'REF', 'altAllele': 'ALT', 'acmgClassification': 'ACMGC_Classification', 'rationale': 'ACMG_Criteria'})
    #print(df)
    return df


def simplify_columns(vcf, ped_data, acmg_f):
    column_map = {}
    role_names = {'mother': 'Mother', 'father': 'Father', 'child': 'Child'}
    relationships = {'child': ped_data['child'], 'mother': ped_data['mother'], 'father': ped_data['father']}
    id_to_role = {v: k for k, v in relationships.items()}
    for col in vcf.columns:
        if isinstance(col, str):
            column_map[col] = col
        else:
            sample_id, field = col
            role = id_to_role.get(sample_id)
            if role:
                new_name = f"{role_names[role]}_{field.upper()}"
                column_map[col] = new_name
            else:
                column_map[col] = col
    # return vcf.rename(columns=column_map)
    # print(vcf)


    # Merge DataFrames on the key columns
    vcf = pd.merge(
        vcf,
        acmg_f[['#CHROM', 'POS', 'REF', 'ALT', 'ACMGC_Classification', 'ACMG_Criteria']],
        on=['#CHROM', 'POS', 'REF', 'ALT'],
        how='left'  # Use 'inner' if you only want rows with matches
    )
    return vcf.rename(columns=column_map)
    print(vcf)
# ------------------ Block 6: Merge VEP + ClinGen ------------------

def merge_vep_clingen(df, df_clin_agg):
    clingen_with_vcf = df.merge(
        df_clin_agg,
        left_on="SYMBOL",
        right_on="Gene Symbol",
        how="left",
        suffixes=('_vcf','_clingen'),
        validate='many_to_one'
    ).drop(columns=['Gene ID','cytoBand','Genomic Location','Date Last Evaluated'], errors='ignore') \
     .drop_duplicates()
    return clingen_with_vcf

# ------------------ Block 7: Merge with VCF ------------------

def merge_with_vcf(clingen_with_vcf, vcf_simple):
    merged_df = pd.merge(
        clingen_with_vcf,
        vcf_simple,
        left_on='Adjusted_Location',
        right_on='CHROM_POS',
        how='left'
    )
    merged_df['ClinVar_CLNDISDB'] = merged_df['ClinVar_CLNDISDB'].replace(r'^\s*$', np.nan, regex=True)
    merged_df['ClinVar_CLNDISDB'] = merged_df['ClinVar_CLNDISDB'].str.replace(
        r'\b(MONDO):\1:', r'\1:', regex=True)
    return merged_df

# ------------------ Block 8: Compound Hets ------------------

def add_compound_hets(merged_df, ped_data, gene_col='Gene'):
    child_id  = ped_data['child']
    mother_id = ped_data['mother']
    father_id = ped_data['father']
    if   'Child_GT' in merged_df.columns:
        gt_col = 'Child_GT'
    elif (child_id, 'GT') in merged_df.columns:
        gt_col = (child_id, 'GT')
    else:
        raise KeyError(f"Could not find GT column for child {child_id}")
    if   'Mother_GT' in merged_df.columns:
        mom_gt_col = 'Mother_GT'
    elif (mother_id, 'GT') in merged_df.columns:
        mom_gt_col = (mother_id, 'GT')
    else:
        raise KeyError(f"Could not find GT column for mother {mother_id}")
    if   'Father_GT' in merged_df.columns:
        dad_gt_col = 'Father_GT'
    elif (father_id, 'GT') in merged_df.columns:
        dad_gt_col = (father_id, 'GT')
    else:
        raise KeyError(f"Could not find GT column for father {father_id}")
    mat = merged_df[merged_df['INHERITANCE']=='ad_maternal']
    pat = merged_df[merged_df['INHERITANCE']=='ad_paternal']
    for gene, mat_sub in mat.groupby(gene_col):
        pat_sub = pat[pat[gene_col]==gene]
        if pat_sub.empty:
            continue
        for mi in mat_sub.index:
            for pi in pat_sub.index:
                in_trans = True
                gt_m   = merged_df.at[mi, gt_col]
                gt_p   = merged_df.at[pi, gt_col]
                mom_gt = merged_df.at[mi, mom_gt_col]
                dad_gt = merged_df.at[pi, dad_gt_col]
                if mom_gt == './.' or dad_gt == './.':
                    continue
                if ('|' in gt_m and '|' in gt_p and '|' in mom_gt and '|' in dad_gt):
                    m0,  m1  = gt_m.split('|')
                    mm0, mm1 = mom_gt.split('|')
                    p0,  p1  = gt_p.split('|')
                    fd0, fd1 = dad_gt.split('|')
                    mat_from_mom = (m1  == '1') and (mm1 == '1')
                    pat_from_dad = (p0  == '1') and (fd0 == '1')
                    in_trans = mat_from_mom and pat_from_dad
                if in_trans:
                    merged_df.at[mi, 'INHERITANCE'] = 'ar_comp_het'
                    merged_df.at[pi, 'INHERITANCE'] = 'ar_comp_het'
    return merged_df

# ------------------ Block 9: Zygosity Classification ------------------

def classify_zygosity(gt, vaf=None, chrom=None, pos=None, sex=None):
    chrom_clean = None
    if isinstance(chrom, str):
        chrom_clean = chrom.lstrip("chr").upper()
    try:
        pos = int(pos)
    except (TypeError, ValueError):
        pos = None
    PAR_RANGES = {
        'PAR1': {'X': (10000, 2781479),      'Y': (10000, 2781479)},
        'PAR2': {'X': (155701382, 156030895),'Y': (56887903, 57217415)},
    }
    in_par = False
    if chrom_clean in {'X', 'Y'} and pos is not None:
        for start, end in PAR_RANGES['PAR1'][chrom_clean], PAR_RANGES['PAR2'][chrom_clean]:
            if start <= pos <= end:
                in_par = True
                break
    if pd.isna(gt):
        alleles = []
    else:
        alleles = re.split(r"[\/|]", str(gt))
        alleles = [a for a in alleles if a != ""]
    if (isinstance(sex, str) and sex.lower() == "male"
        and chrom_clean in {"X", "Y"} and not in_par):
        if len(alleles) == 1 and alleles[0] != ".":
            return "Hemizygous"
    if len(alleles) == 2 and "." not in alleles:
        return "Homozygous" if alleles[0] == alleles[1] else "Heterozygous"
    if vaf is not None and not pd.isna(vaf):
        try:
            v = float(vaf)
        except (TypeError, ValueError):
            return "Unclassified"
        if v >= 0.85 or v <= 0.15:
            return "Homozygous"
        if 0.35 <= v <= 0.65:
            return "Heterozygous"
    return "Unclassified"

def add_zygosity(merged_df):
    for col in ["Child_DP", "Child_GQ", "Child_VAF"]:
        merged_df[col] = pd.to_numeric(merged_df[col], errors="coerce")
    merged_df["POS"] = pd.to_numeric(merged_df["POS"], errors="coerce")
    # You may need to add a column for Child_sex from ped_data
    merged_df["Zygosity"] = merged_df.apply(
        lambda r: classify_zygosity(
            gt    = r["Child_GT"],
            vaf   = r.get("Child_VAF", None),
            chrom = r.get("#CHROM", r.get("CHROM", None)),
            pos   = r.get("POS", None),
            sex   = r.get("Child_sex", None)
        ),
        axis=1
    )
    return merged_df

# ------------------ Block 10: Final Labels and Output ------------------

def finalize_labels(merged_df, zmap):
    merged_df.insert(
        merged_df.columns.get_loc("Zygosity")+1,
        "Zygosity_label",
        merged_df["Zygosity"].map(zmap)
    )
    clinvar_map = {
        "likely_benign": "LB",
        "benign": "Ben",
        "likely_pathogenic": "LP",
        "pathogenic": "PAT",
        "uncertain_significance": "VUS",
        "conflicting_interpretations_of_pathogenicity": "Conflicting"
    }
    def map_clinsig(c):
        if pd.isna(c): return "NA"
        pat = r'([^/,]+)'
        def rep(m):
            term = m.group(0).strip().lower()
            return clinvar_map.get(term, term)
        return re.sub(pat, rep, c)
    merged_df.insert(
        merged_df.columns.get_loc('ClinVar_CLNSIG')+1,
        "CLINSIG_label",
        merged_df['ClinVar_CLNSIG'].apply(map_clinsig)
    )
    def merge_prioritized(r):
        if pd.notna(r.get('OMIM', None)):
            return f"OMIM:{r.get('OMIM')}"
        elif pd.notna(r.get('Orpha', None)):
            return f"Orpha:{r.get('Orpha')}"
        elif pd.notna(r.get('MONDO', None)):
            return f"MONDO:{r.get('MONDO')}"
        elif pd.notna(r.get('MedGen', None)):
            return f"MedGen:{r.get('MedGen')}"
        return ""
    # def merge_prioritized(r):
    #     if pd.notna(r['OMIM']):    return f"OMIM:{r['OMIM']}"
    #     elif pd.notna(r['Orpha']): return f"Orpha:{r['Orpha']}"
    #     elif pd.notna(r['MONDO']): return f"Orpha:{r['MONDO']}"
    #     elif pd.notna(r['MedGen']):    return f"MedGen:{r['MedGen']}"
    #     return ""
    merged_df['ACMG DISEASE ID (OMIM/ORPHA ID)'] = merged_df.apply(merge_prioritized, axis=1)
#    merged_df = merged_df.rename(columns={
#        'Gene': 'GENE (GENE ID)',
#        'SYMBOL': 'SYMBOL (Gene Name)',
#        'Existing_variation': 'SNPs/Rsid',
#        'ClinVar_CLNSIG': 'CLINVAR CLNSIG',
#        'Adjusted_Location': 'LOCATION',
#        'Consequence': 'EFFECT'
#    })
    return merged_df

# ------------------ Main Automation Script ------------------

def main():
    parser = argparse.ArgumentParser(description="Automated trio variant annotation pipeline")
    parser.add_argument('--vep', required=True, help='VEP-annotated file')
    parser.add_argument('--clingen', required=True, help='ClinGen gene curation file')
    parser.add_argument('--ped', required=True, help='PED file')
    parser.add_argument('--vcf', required=True, help='Phased VCF file (gzipped)')
    parser.add_argument('--acmg', required=True, help='ACMG TSV file')
    parser.add_argument('--output', required=True, help='Output file (TSV)')
    args = parser.parse_args()

    print("Loading VEP...")
    df = load_vep(args.vep)
    print("Loading ClinGen...")
    df_clin_agg = load_clingen(args.clingen, df)
    print("Loading PED...")
    ped_data = load_ped_file(args.ped)
    child_sex = ped_data['sex_map'][('child', ped_data['child'])]
    with open("child_sex.txt", "w") as f:
        f.write(child_sex)
    print("Loading VCF...")
    vcf = read_vcf(args.vcf)
    print("Tagging Mendelian inheritance...")
    vcf = add_mendelian_inheritance(vcf, ped_data)

    acmg_f = read_acmg(args.acmg)
    print("Loading ACMG file")
    vcf_simple = simplify_columns(vcf, ped_data, acmg_f)
    print(vcf_simple)
    print("Merging VEP and ClinGen...")
    clingen_with_vcf = merge_vep_clingen(df, df_clin_agg)
    print("Merging with VCF...")
    merged_df = merge_with_vcf(clingen_with_vcf, vcf_simple)
    print("Flagging compound hets...")
    merged_df = add_compound_hets(merged_df, ped_data)
    print("Computing zygosity...")
    merged_df = add_zygosity(merged_df)
    # You must define zmap for Zygosity_label mapping
    zmap = {"Homozygous": "Hom", "Heterozygous": "Het", "Hemizygous": "Hem", "Unclassified": "Unc"}
    print("Finalizing labels...")
    merged_df = finalize_labels(merged_df, zmap)
    print(f"Saving output to {args.output}")
    merged_df.to_csv(args.output, sep="\t", index=False)
    print("Done.")

if __name__ == "__main__":
    main()

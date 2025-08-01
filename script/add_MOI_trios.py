#!/usr/bin/env python3
"""
Standalone script for annotating variants with Mode of Inheritance (MOI) information
by combining ClinVar, OMIM, and MONDO ontology data.
"""

import sys
import pandas as pd
import re
from owlready2 import get_ontology
from typing import List, Optional, Dict, Set

def normalize_phenotype(s: str) -> Optional[str]:
    """Normalize phenotype strings for consistent matching."""
    if pd.isna(s) or s in ("", "-"):
        return None
    t = s.replace("_", " ").replace("-", " ")
    t = t.lower()
    t = re.sub(r"[^a-z0-9\s]", " ", t)
    t = re.sub(r"\s+", " ", t).strip()
    return t or None

def split_terms(field_str: str) -> List[str]:
    """Split multi-term fields into individual terms."""
    if pd.isna(field_str):
        return []
    return [t.strip() for t in re.split(r"[|,;/]+", str(field_str)) if t.strip()]

def load_phenotype_to_moi_map(omim_file: str) -> Dict[str, List[str]]:
    """Create phenotype to MOI mapping from OMIM data."""
    omim = pd.read_csv(omim_file, sep="\t", dtype=str)
    omim["norm_pheno"] = omim["Phenotype"].apply(normalize_phenotype)
    omim = omim[omim["norm_pheno"].notna()]
    return (
        omim
        .dropna(subset=["MOI", "norm_pheno"])
        .groupby("norm_pheno")["MOI"]
        .unique()
        .to_dict()
    )

def fetch_mondo_synonyms(obo_file: str, curie: str) -> List[str]:
    """Retrieve synonyms for a MONDO term from ontology."""
    onto = get_ontology(f"file://{obo_file}").load()
    suf = curie.replace(":", "_")
    cls = onto.search_one(iri=f"*{suf}")
    if not cls:
        return []
    
    syns = list(getattr(cls, "label", []))
    for attr in ("hasExactSynonym", "hasRelatedSynonym", "hasBroadSynonym"):
        syns += list(getattr(cls, attr, []))
    return syns

def find_moi_for_variant(
    row: pd.Series,
    pheno2moi: Dict[str, List[str]],
    obo_file: str
) -> Optional[str]:
    """Determine MOI for a single variant by checking multiple evidence sources."""
    mois = set()
    
    # Check ClinVar disease names
    clin_terms = split_terms(row.get("ClinVar_CLNDN", ""))
    for term in clin_terms:
        n = normalize_phenotype(term)
        if n and n in pheno2moi:
            mois.update(pheno2moi[n])
    
    # Check VEP PHENOTYPES if no MOI found yet
    if not mois:
        pheno_terms = split_terms(row.get("PHENOTYPES", ""))
        for term in pheno_terms:
            n = normalize_phenotype(term)
            if n and n in pheno2moi:
                mois.update(pheno2moi[n])
    
    # Check MONDO terms if still no MOI found
    if not mois:
        mondo_terms = split_terms(row.get("MONDO", ""))
        for curie in mondo_terms:
            for syn in fetch_mondo_synonyms(obo_file, curie):
                n = normalize_phenotype(syn)
                if n and n in pheno2moi:
                    mois.update(pheno2moi[n])
    
    return ";".join(sorted(mois)) if mois else None

def main(
    vep_file: str = "all_zyg_g2p.csv",
    omim_file: str = "omim_clingene.tsv",
    obo_file: str = "mondo.owl",
    output_file: str = "trioPVS1_with_MOI_g2p.csv"
):
    """Main execution function for the script."""
    try:
        # Load and prepare data
        vep = pd.read_csv(vep_file, sep=",", dtype=str)
        pheno2moi = load_phenotype_to_moi_map(omim_file)
        
        # Annotate MOI
        vep["MOI"] = vep.apply(
            lambda row: find_moi_for_variant(row, pheno2moi, obo_file),
            axis=1
        )
        
        # Save results
        vep.to_csv(output_file, sep=",", index=False)
        print(f"Successfully annotated MOI for {vep['MOI'].notna().sum():,} variants")
        print(f"Results saved to: {output_file}")
        return 0
    
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Allow custom file paths via command line
        sys.exit(main(*sys.argv[1:5]))
    else:
        # Run with default file paths
        sys.exit(main())
# Clinical WES Trio Analysis Pipeline

![DeepVariant](https://img.shields.io/badge/DeepVariant-1.6.0-green)
![VEP](https://img.shields.io/badge/VEP-109-blue)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A production-ready pipeline for clinical exome trio analysis, from raw FASTQs to prioritized variant lists with ACMG classification.

## Key Features
- **Trio-Optimized**:
  - DeepVariant joint calling
  - BIAS for inheritance pattern checks
- **Structural Variants**:
  - Manta (SVs)
  - ExomeDepth (CNVs)
- **Clinical Annotation**:
  - VEP + AnnotSV (pathogenicity scores)
  - Custom Python prioritization (ACMG rules)
- **QC Metrics**:
  - Pedigree consistency checks
  - Coverage reports (20x/50x/100x)

## Quick Start
```bash
nextflow run trio-wes-pipeline/main.nf \
  --pedigree family.ped \
  --fastq_dir ./fastqs \
  --outdir results \
  -profile docker

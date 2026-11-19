# MATLAB Computational Workflow

## Overview

This repository contains a modular MATLAB workflow developed to investigate splice-site responses to the splice-modifying compounds risdiplam and branaplam.

The analysis integrates RNA splicing measurements, nucleotide sequence characteristics, regulatory features, statistical analysis, biological response classification, and machine learning.

The workflow consists of multiple MATLAB scripts executed as separate analytical steps rather than a single executable pipeline.

## Workflow Organization

The MATLAB scripts are organized into six computational modules.

| Module | Directory | Purpose |
|---|---|---|
| 1 | `01_sequence_preprocessing` | Sequence extraction and preprocessing |
| 2 | `02_feature_engineering` | PSI feature extraction and feature integration |
| 3 | `03_statistical_analysis` | Statistical feature computation |
| 4 | `04_biological_classification` | Drug-response classification and group comparisons |
| 5 | `05_sequence_functional_analysis` | Sequence-level functional analysis |
| 6 | `06_machine_learning` | Machine learning model training and evaluation |

## MATLAB Script Execution Order

The following sequence reflects the numbered research scripts.

### Stage 1: Sequence Preprocessing

Directory: `code/01_sequence_preprocessing/`

1. `00_extract_smn2_exon7.m`
2. `01_extract_splice_regulatory_features.m`

**Purpose:** Prepare nucleotide sequence information and extract regulatory features for downstream analysis.

### Stage 2: Feature Engineering

Directory: `code/02_feature_engineering/`

3. `02_extract_psi_features.m`
4. `03_merge_psi_and_regulatory_features.m`

**Purpose:** Extract PSI-based measurements and integrate them with sequence regulatory features.

### Stage 3: Statistical Feature Computation

Directory: `code/03_statistical_analysis/`

5. `04_compute_statistical_features.m`

**Purpose:** Compute statistical features used in subsequent analyses.

Return to `code/02_feature_engineering/`:

6. `05_merge_statistical_features.m`

**Purpose:** Integrate the computed statistical features with the existing feature dataset.

### Stage 4: Biological Response Classification

Directory: `code/04_biological_classification/`

7. `06_define_biological_classes.m`
8. `07_filter_branaplam_classes.m`
9. `07_filter_risdiplam_classes.m`
10. `08_compare_smn2_elp1_branaplam.m`
11. `08_compare_smn2_elp1_risdiplam.m`
12. `09_create_contrast_groups_branaplam.m`
13. `09_create_contrast_groups_risdiplam.m`

**Purpose:** Define drug-response categories, compare splice-site groups, and construct contrast groups for further analysis.

The branaplam and risdiplam scripts represent separate treatment-specific analyses.

### Stage 5: Sequence Functional Analysis

Directory: `code/05_sequence_functional_analysis/`

14. `10_sequence_functional_analysis.m`

**Purpose:** Investigate sequence-associated functional characteristics of the analyzed splice sites.

### Stage 6: Machine Learning

Directory: `code/06_machine_learning/`

15. `12_train_evaluate_models.m`

**Purpose:** Train and evaluate machine learning models using the prepared analytical features.

## Required Input Data

The primary source dataset is GEO accession GSE221868.

The locally prepared data directory contains:

```text
data/
├── raw/
│   ├── GSE221868_smn2_psi.csv
│   ├── GSE221868_elp1_psi.csv
│   └── SMN2_ensembl.fasta
└── intermediate/
```

Additional information about the data sources is provided in `docs/DATA.md`.

## Execution Notes

- MATLAB scripts are organized by computational purpose.
- Script numbering indicates the intended broad execution sequence.
- Some analytical stages contain separate scripts for risdiplam and branaplam.
- Intermediate outputs may be required by subsequent scripts.
- Input and output paths should be configured for the local computing environment.
- MATLAB dependencies and required toolboxes should be checked before execution.

## Reproducibility

This repository preserves the modular organization of the research workflow.

Researchers wishing to reproduce the analysis should obtain the original datasets, prepare the required input files, and execute the MATLAB scripts in the documented order.

Exact intermediate-file dependencies, software versions, and numerical reproducibility should be verified against the MATLAB source scripts and original experimental environment.

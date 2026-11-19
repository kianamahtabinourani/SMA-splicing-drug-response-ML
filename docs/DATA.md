# Dataset Documentation

## Data Source

This research uses publicly available RNA splicing measurements from the Gene Expression Omnibus (GEO).

**GEO Accession:** GSE221868

**Dataset:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221868

The dataset provides splice-site measurements under three experimental conditions:

- DMSO (control)
- Risdiplam (RG)
- Branaplam (NVS)

## Primary Input Files

| Filename | Description |
|---|---|
| `GSE221868_smn2_psi.csv` | SMN2 splice-site PSI measurements |
| `GSE221868_elp1_psi.csv` | ELP1 splice-site PSI measurements |
| `SMN2_ensembl.fasta` | Reference nucleotide sequences used for SMN2 exon 7 extraction |

The SMN2 and ELP1 datasets contain splice-site sequence identifiers and experimental measurements across drug-treatment conditions.

The reference FASTA file was obtained separately from Ensembl and is used during sequence preprocessing.

## Local Data Organization

Store the input datasets in the following directory:

```text
data/
├── raw/
│   ├── GSE221868_smn2_psi.csv
│   ├── GSE221868_elp1_psi.csv
│   └── SMN2_ensembl.fasta
└── intermediate/
    └── Generated preprocessing outputs
```

## Data Processing

The MATLAB research workflow includes:

1. Sequence preprocessing
2. PSI-based feature extraction
3. Statistical analysis
4. Biological response classification
5. Sequence functional analysis
6. Machine learning and model evaluation

Intermediate outputs are generated during the execution of individual MATLAB scripts.

## Data Availability

The original GEO datasets are publicly accessible through accession GSE221868.

Users should retrieve the required datasets from their original sources and preserve the filenames expected by the MATLAB scripts.

The processed outputs depend on the preprocessing steps and analysis parameters described in the repository.

## Reproducibility

The repository organizes MATLAB scripts according to their analytical purpose. The scripts were developed as separate research steps rather than a single executable pipeline.

Users should follow the documented execution order and ensure that required input files are available before running each script.

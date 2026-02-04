# MPLID: Membrane Protein-Lipid Interface Dataset

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.XXXXXXX-blue)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Data: CC-BY 4.0](https://img.shields.io/badge/Data-CC--BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

A large-scale, experimentally-validated dataset of membrane protein-lipid contact residues derived from crystallographic structures in the Protein Data Bank.

## Overview

MPLID provides residue-level annotations of lipid contacts for 2,792 membrane proteins, representing the largest experimentally-validated dataset of its kind. All labels are derived from crystallized lipids in PDB structures using a 4.0 Angstrom distance cutoff.

### Key Statistics

| Metric | Value |
|--------|-------|
| Proteins | 3,192 |
| Total residues | 5,134,242 |
| Contact residues | 38,435 |
| Contact rate | 0.75% |
| Label source | Experimental (crystallized lipids) |
| Sequence clusters (30% identity) | 594 |

### Dataset Splits

| Split | Proteins | Residues |
|-------|----------|----------|
| Training | 1,840 | 2,634,209 |
| Validation | 811 | 1,632,603 |
| Test | 541 | 867,430 |

Splits are performed at the protein level using sequence clustering to prevent data leakage.

## Installation

```bash
git clone https://github.com/omagebright/MPLID.git
cd MPLID
pip install -r requirements.txt
```

### Dependencies

- Python >= 3.8
- pandas >= 1.3.0
- numpy >= 1.20.0
- biopython >= 1.79
- scikit-learn >= 1.0.0
- matplotlib >= 3.4.0
- seaborn >= 0.11.0

## Quick Start

### Loading the Dataset

```python
import pandas as pd

# Load training data (compressed format)
train = pd.read_csv('data/processed/train_residues.csv.gz', compression='gzip')

print(f"Training samples: {len(train):,}")
print(f"Positive samples: {train['is_contact'].sum():,}")
print(f"Positive rate: {train['is_contact'].mean():.2%}")

# Filter for contact residues
contacts = train[train['is_contact'] == True]
print(f"\nContact residues: {len(contacts):,}")
print(f"Unique proteins with contacts: {contacts['pdb_id'].nunique():,}")
```

### Data Format

Each compressed CSV file (`.csv.gz`) contains the following columns:

| Column | Type | Description |
|--------|------|-------------|
| `pdb_id` | str | PDB identifier (4 characters) |
| `chain_id` | str | Chain identifier |
| `residue_number` | int | Residue position in chain |
| `residue_name` | str | Three-letter amino acid code |
| `is_contact` | bool | Contact label (True/False) |
| `label_source` | str | Label origin (EXPERIMENTAL) |
| `confidence` | str | Label confidence level |
| `min_distance` | float | Minimum distance to lipid (Angstroms) |
| `cluster_id` | int | Sequence cluster assignment |
| `split` | str | Dataset split (train/val/test) |

## Repository Structure

```
MPLID/
├── data/
│   ├── processed/           # Ready-to-use dataset files
│   │   ├── train_residues.csv
│   │   ├── val_residues.csv
│   │   ├── test_residues.csv
│   │   └── protein_metadata.csv
│   └── statistics/          # Dataset statistics
│       ├── lipid_distribution.csv
│       └── dataset_summary.json
├── scripts/                 # Data processing pipeline
│   ├── 01_query_rcsb_lipids.py
│   ├── 02_download_structures.py
│   ├── 03_calculate_contacts.py
│   ├── 04_cluster_sequences.py
│   ├── 05_create_splits.py
│   └── utils/
├── notebooks/               # Analysis notebooks
│   ├── 01_dataset_exploration.ipynb
│   └── 02_statistical_analysis.ipynb
└── docs/                    # Documentation
    ├── METHODOLOGY.md
    └── DATA_DICTIONARY.md
```

## Methodology

### Data Sources

1. **OPM Database**: Membrane protein identification and positioning
2. **RCSB PDB**: Structure coordinates and lipid ligand information

### Contact Definition

A residue is labeled as a lipid contact if any of its heavy atoms are within **4.0 Angstroms** of any heavy atom of a crystallized lipid molecule (HETATM records).

### Supported Lipid Types

The dataset recognizes 150+ lipid codes including:

- **Phospholipids**: CDL, POV, PCW, PEE, PGV, PLM
- **Sterols**: CLR, CHD, Y01, BCL
- **Sphingolipids**: SPH, S1P
- **Fatty acids**: MYR, OLA, STE, ARA
- **Detergents**: LDA, LMT, BOG, OLC

See `docs/DATA_DICTIONARY.md` for the complete list.

### Quality Control

1. **Sequence clustering**: 30% identity threshold using MMseqs2
2. **Stratified splitting**: Balanced contact rates across splits
3. **No data leakage**: Proteins from same cluster in same split

## Reproducing the Dataset

To regenerate the dataset from scratch:

```bash
# Step 1: Query RCSB for lipid-containing structures
python scripts/01_query_rcsb_lipids.py --output data/interim/

# Step 2: Download PDB structures
python scripts/02_download_structures.py --input data/interim/rcsb_lipids.json

# Step 3: Calculate contacts
python scripts/03_calculate_contacts.py --cutoff 4.0

# Step 4: Cluster sequences
python scripts/04_cluster_sequences.py --identity 0.3

# Step 5: Create train/val/test splits
python scripts/05_create_splits.py --train 0.6 --val 0.2 --test 0.2
```

## Comparison with Existing Datasets

| Dataset | Proteins | Label Source | Public |
|---------|----------|--------------|--------|
| DREAMM | 54 | Experimental | Limited |
| OPM-derived | 15,096 | Computational | Yes |
| **MPLID** | **2,792** | **Experimental** | **Yes** |

MPLID provides 52x more proteins than DREAMM while maintaining experimental validation.

## Citation

If you use this dataset in your research, please cite:

```bibtex
@article{omage2026mplid,
  author = {Omage, Folorunsho Bright},
  title = {{MPLID}: A Large-Scale Experimentally-Validated Dataset of
           Membrane Protein-Lipid Interface Residues},
  journal = {GigaScience},
  year = {2026},
  volume = {XX},
  pages = {XXX},
  doi = {10.1093/gigascience/XXXX}
}
```

## License

- **Code**: MIT License (see [LICENSE](LICENSE))
- **Data**: [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/)

## Acknowledgments

This research was funded by the São Paulo Research Foundation (FAPESP) (grant nos. 2023/02691-2, 2025/23708-6).

## Contact

**Folorunsho Bright Omage**
- Email: omagefolorunsho@gmail.com
- ORCID: [0000-0002-9750-5034](https://orcid.org/0000-0002-9750-5034)

---

*For questions or issues, please open a GitHub issue or contact the author directly.*

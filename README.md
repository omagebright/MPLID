# MPLID: Membrane Protein-Lipid Interface Dataset

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18487584.svg)](https://doi.org/10.5281/zenodo.18487584)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Data: CC0](https://img.shields.io/badge/Data-CC0%201.0-lightgrey.svg)](https://creativecommons.org/publicdomain/zero/1.0/)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

A large-scale dataset of experimentally validated lipid contact residues derived from experimentally determined structures in the Protein Data Bank. **100% experimental labels with no computational database dependencies.**

## Overview

MPLID provides residue-level annotations of lipid contacts for **4,704 proteins** with experimentally resolved lipid molecules, representing the largest experimentally validated dataset of its kind. All labels are derived exclusively from lipid molecules resolved in PDB structures (X-ray crystallography, cryo-EM) using a 4.0 Å all-atom heavy-atom distance cutoff.

### Key Statistics

| Metric | Value |
|--------|-------|
| Proteins | 4,704 |
| Total residues | 8,055,325 |
| Contact residues | 80,439 |
| Contact rate | 1.00% |
| Label source | Experimental (crystallized lipids) |
| Sequence clusters (30% identity) | 813 |
| Lipid codes recognized | 117 |

### Dataset Splits

| Split | Proteins | Residues | Contacts | Rate |
|-------|----------|----------|----------|------|
| Training | 2,578 | 4,907,696 | 56,976 | 1.16% |
| Validation | 1,051 | 1,403,838 | 9,626 | 0.69% |
| Test | 1,075 | 1,743,791 | 13,837 | 0.79% |

Splits are performed at the cluster level using 30% sequence identity to prevent data leakage.

## What Makes MPLID Different

### 100% Experimental Labels

MPLID derives labels **directly from experimentally determined PDB coordinates** (X-ray crystallography, cryo-EM). A residue is labeled positive if any of its heavy atoms (non-hydrogen) is within 4.0 Å of any heavy atom of an experimentally resolved lipid molecule. No computational predictions, no membrane plane algorithms, no literature curation required.

### Comparison with Other Datasets

MPLID and other protein-lipid datasets answer **different biological questions**:

| Aspect | DREAMM | MPLID |
|--------|--------|-------|
| Label source | Literature curation (EPR, fluorescence, mutagenesis) | Experimental (PDB structures) |
| Biological question | Functional membrane contact | Structural lipid proximity |
| Protein scope | Peripheral membrane proteins | All proteins with crystallized lipids |
| Proteins | 54 | 4,704 |
| Reproducibility | Requires literature access | Fully automated from PDB |

**Important**: DREAMM identifies residues that functionally interact with membranes through biophysical experiments; MPLID identifies residues structurally proximate to crystallized lipids. These are complementary approaches, not directly comparable benchmarks. Models trained on one may not transfer to the other.

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
- requests >= 2.25.0

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

Each CSV file contains the following columns:

| Column | Type | Description |
|--------|------|-------------|
| `pdb_id` | str | PDB identifier (4 characters) |
| `chain_id` | str | Chain identifier |
| `residue_number` | int | Residue position in chain |
| `residue_name` | str | Three-letter amino acid code |
| `is_contact` | bool | Contact label (True if any heavy atom ≤4.0 Å from lipid) |
| `label_source` | str | Always "EXPERIMENTAL" |
| `confidence` | str | Label confidence level |
| `min_distance` | float | Minimum heavy-atom-to-lipid distance (Å) |
| `cluster_id` | int | Sequence cluster assignment |
| `split` | str | Dataset split (train/val/test) |

## Repository Structure

```
MPLID/
├── data/
│   ├── processed/           # Ready-to-use dataset files
│   │   ├── train_residues.csv.gz
│   │   ├── val_residues.csv.gz
│   │   ├── test_residues.csv.gz
│   │   └── protein_metadata.csv
│   └── statistics/          # Dataset statistics
│       └── pipeline_stats.json
├── scripts/                 # Data processing pipeline
│   └── run_rcsb_only_pipeline.py
├── notebooks/               # Analysis notebooks
│   └── dataset_exploration.ipynb
└── docs/                    # Documentation
    ├── METHODOLOGY.md
    └── DATA_DICTIONARY.md
```

## Methodology

### Data Source

MPLID uses a single data source: the **RCSB Protein Data Bank**. We query RCSB for all structures containing recognized lipid chemical components (117 codes), download the original PDB files preserving HETATM records, and compute contacts directly from atomic coordinates.

**No OPM database, no membrane plane calculations, no external computational dependencies.**

### Contact Definition

A residue is labeled as a lipid contact if:

```
d_min(R, L) = min_{a ∈ R_heavy, b ∈ L_heavy} ||r_a - r_b|| ≤ 4.0 Å
```

Where `R_heavy` is the set of all heavy atoms (non-hydrogen) of residue R and `L_heavy` is the set of all heavy atoms of lipid molecule L.

### Supported Lipid Types (117 codes)

- **Phospholipids** (24): PC1, PCW, POV, PEE, PGV, LHG, PIP, etc.
- **Cardiolipin** (4): CDL, C9V, 18W, LCL
- **Sphingolipids** (7): SPH, S1P, HXJ, CRT, GCR, GSP, GLF
- **Sterols** (9): CLR, CHD, Y01, ERG, SIT, STI, etc.
- **Fatty acids** (14): PLM, MYR, OLA, STE, ARA, DHA, etc.
- **Glycerolipids** (6): TGL, DAG, MAG, GMO, etc.
- **Detergents** (22): LDA, LMT, BOG, OLC, DPC, SDS, etc.
- **Lipid A** (3): LPA, KDO, 6LP
- **CHARMM simulation** (27): POPC, POPE, POPG, etc.
- **Generic** (1): LIP

See `docs/DATA_DICTIONARY.md` for the complete list.

### Quality Control

1. **Sequence clustering**: 30% identity threshold using MMseqs2
2. **Cluster-level splitting**: Entire clusters assigned to same split
3. **No data leakage**: No protein in test shares >30% identity with training

## Reproducing the Dataset

```bash
# Single command to regenerate the entire dataset from RCSB
python scripts/run_rcsb_only_pipeline.py

# With options
python scripts/run_rcsb_only_pipeline.py --skip-rcsb-query  # Use cached RCSB results
python scripts/run_rcsb_only_pipeline.py --max-proteins 100  # Test with subset
```

## Limitations

1. **Crystallization artifacts**: Some lipid contacts may reflect crystallization conditions rather than native interactions
2. **Detergents included**: Dataset includes detergent molecules used in crystallization
3. **Static structures**: Does not capture dynamic or transient lipid interactions
4. **Class imbalance**: 1.00% positive rate (100:1 ratio) requires appropriate ML techniques

## Citation

If you use this dataset, please cite:

```bibtex
@article{omage2026mplid,
  author = {Omage, Folorunsho Bright and Neshich, Goran},
  title = {{MPLID} ({Membrane Protein-Lipid Interaction Database}): A Large-Scale Experimental Resource of Residue-Level Protein-Lipid Contacts},
  journal = {GigaScience},
  year = {2026},
  doi = {10.5281/zenodo.18487584}
}
```

## License

- **Code**: MIT License (see [LICENSE](LICENSE))
- **Data**: [CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/) (Public Domain)

## Acknowledgments

This research was funded by the São Paulo Research Foundation (FAPESP) (grant nos. 2023/02691-2, 2025/23708-6).

## Contact

**Folorunsho Bright Omage**
- Email: omagefolorunsho@gmail.com
- ORCID: [0000-0002-9750-5034](https://orcid.org/0000-0002-9750-5034)

---

*For questions or issues, please open a GitHub issue or contact the author directly.*

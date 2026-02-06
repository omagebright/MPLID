# MPLID: Membrane Protein-Lipid Interface Dataset

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18487585.svg)](https://doi.org/10.5281/zenodo.18487585)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Data: CC0](https://img.shields.io/badge/Data-CC0%201.0-lightgrey.svg)](https://creativecommons.org/publicdomain/zero/1.0/)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

A large-scale dataset of experimentally validated lipid contact residues derived from crystallographic structures in the Protein Data Bank. **100% crystallographic labels with no computational database dependencies.**

## Overview

MPLID provides residue-level annotations of lipid contacts for **7,541 proteins** with crystallized lipid molecules, representing the largest experimentally validated dataset of its kind. All labels are derived exclusively from crystallized lipids in PDB structures using a 4.0 Å Cα-to-lipid distance cutoff.

### Key Statistics

| Metric | Value |
|--------|-------|
| Proteins | 7,541 |
| Total residues | 10,486,965 |
| Contact residues | 98,101 |
| Contact rate | 0.94% |
| Label source | Experimental (crystallized lipids) |
| Sequence clusters (30% identity) | 1,848 |
| Lipid codes recognized | 121 |

### Dataset Splits

| Split | Proteins | Residues | Contacts | Rate |
|-------|----------|----------|----------|------|
| Training | 4,644 | 7,521,766 | 73,432 | 0.98% |
| Validation | 1,436 | 1,507,429 | 10,879 | 0.72% |
| Test | 1,413 | 1,391,166 | 12,980 | 0.93% |

Splits are performed at the cluster level using 30% sequence identity to prevent data leakage.

## What Makes MPLID Different

### 100% Crystallographic Labels

MPLID derives labels **directly from PDB crystallographic coordinates**. A residue is labeled positive if its Cα atom is within 4.0 Å of any heavy atom of a crystallized lipid molecule. No computational predictions, no membrane plane algorithms, no literature curation required.

### Comparison with DREAMM

| Aspect | DREAMM | MPLID |
|--------|--------|-------|
| Label source | Literature curation (EPR, fluorescence, mutagenesis) | PDB crystallography |
| Biological question | Functional membrane contact | Structural lipid proximity |
| Proteins | 54 | 7,541 (140× more) |
| Residues | ~15,000 | 10,486,965 (699× more) |
| Reproducibility | Requires literature access | Fully automated from PDB |

Both datasets use experimental data but answer **different biological questions**. DREAMM identifies residues that functionally interact with membranes; MPLID identifies residues structurally proximate to crystallized lipids.

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
| `is_contact` | bool | Contact label (True if ≤4.0 Å from lipid) |
| `label_source` | str | Always "EXPERIMENTAL" |
| `confidence` | str | Label confidence level |
| `min_distance` | float | Minimum Cα-to-lipid distance (Å) |
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

MPLID uses a single data source: the **RCSB Protein Data Bank**. We query RCSB for all structures containing recognized lipid chemical components (121 codes), download the original PDB files preserving HETATM records, and compute contacts directly from atomic coordinates.

**No OPM database, no membrane plane calculations, no external computational dependencies.**

### Contact Definition

A residue is labeled as a lipid contact if:

```
d_min(R, L) = min_{atom ∈ L} ||r_Cα - r_atom|| ≤ 4.0 Å
```

Where `r_Cα` is the Cα atom coordinate and `r_atom` is any heavy atom of lipid molecule L.

### Supported Lipid Types (121 codes)

- **Phospholipids**: CDL, POV, PCW, PEE, PGV, PLM, etc.
- **Sterols**: CLR, CHD, Y01, BCL
- **Sphingolipids**: SPH, S1P, HXJ
- **Fatty acids**: MYR, OLA, STE, ARA, DHA
- **Detergents**: LDA, LMT, BOG, OLC, DPC

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
4. **Class imbalance**: 0.94% positive rate requires specialized ML techniques

## Citation

If you use this dataset, please cite:

```bibtex
@article{omage2026mplid,
  author = {Omage, Folorunsho Bright and Mazoni, Ivan and Yano, Inácio Henrique and Neshich, Goran},
  title = {{MPLID}: A Large-Scale Dataset of Experimentally Validated Membrane Protein-Lipid Contact Residues for Machine Learning},
  journal = {GigaScience},
  year = {2026},
  doi = {10.5281/zenodo.18487585}
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

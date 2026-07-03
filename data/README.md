# MPLID Dataset

This directory contains the processed MPLID (Membrane Protein-Lipid Interface Dataset) files.

## Directory Structure

```
data/
├── processed/
│   ├── train_residues.csv.gz    # Training split (compressed)
│   ├── val_residues.csv.gz      # Validation split (compressed)
│   ├── test_residues.csv.gz     # Test split (compressed)
│   ├── protein_metadata.csv     # Per-protein statistics
│   └── dataset_summary.json     # Overall dataset statistics
└── statistics/
    └── (statistical summaries)
```

## File Descriptions

### Residue Files (train/val/test_residues.csv.gz)

Gzip-compressed CSV files containing residue-level annotations.

**Columns:**
| Column | Type | Description |
|--------|------|-------------|
| `pdb_id` | str | PDB structure identifier (4 characters) |
| `chain_id` | str | Protein chain identifier |
| `residue_number` | int | Residue position in chain |
| `residue_name` | str | Three-letter amino acid code |
| `is_contact` | bool | Lipid contact label (True/False) |
| `label_source` | str | Label origin (EXPERIMENTAL) |
| `confidence` | str | Label confidence level |
| `min_distance` | float | Minimum distance to lipid (Angstroms) |
| `cluster_id` | int | Sequence cluster assignment |
| `split` | str | Dataset split (train/val/test) |

### protein_metadata.csv

Per-protein summary statistics.

**Columns:**
| Column | Type | Description |
|--------|------|-------------|
| `pdb_id` | str | PDB identifier |
| `chain_id` | str | Primary chain |
| `cluster_id` | int | Sequence cluster |
| `split` | str | Dataset split |
| `n_residues` | int | Total residues |
| `n_contacts` | int | Contact residues |
| `contact_rate` | float | Fraction of contacts |

### dataset_summary.json

Overall dataset statistics in JSON format.

## Usage

### Loading Compressed Files

```python
import pandas as pd

# Load training data
train = pd.read_csv('data/processed/train_residues.csv.gz', compression='gzip')

print(f"Training samples: {len(train):,}")
print(f"Positive samples: {train['is_contact'].sum():,}")
print(f"Positive rate: {train['is_contact'].mean():.2%}")
```

### Quick Statistics

```python
import json

with open('data/processed/dataset_summary.json') as f:
    stats = json.load(f)

print(f"Total proteins: {stats['total_proteins']:,}")
print(f"Total residues: {stats['total_residues']:,}")
print(f"Contact rate: {stats['contact_rate']:.2%}")
```

## Dataset Statistics

| Metric | Value |
|--------|-------|
| Total proteins | 4,704 |
| Total residues | 8,055,325 |
| Contact residues | 80,439 |
| Contact rate | 1.00% |
| Sequence clusters | 813 |

### Split Distribution

| Split | Proteins | Residues |
|-------|----------|----------|
| Train | 2,578 | 4,907,696 |
| Validation | 1,051 | 1,403,838 |
| Test | 1,075 | 1,743,791 |

## License

Data: CC0 1.0 Universal (Public Domain Dedication) — https://creativecommons.org/publicdomain/zero/1.0/

This matches the CC0 dedication stated in the associated manuscript and the Zenodo deposit. Code in this repository is licensed separately under the MIT License (see `LICENSE`).

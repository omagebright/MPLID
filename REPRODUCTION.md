# MPLID Reproduction Guide

This document provides step-by-step instructions to reproduce the MPLID dataset from scratch using only publicly available data from RCSB PDB.

## System Requirements

### Hardware
- **Minimum**: 16 GB RAM, 100 GB disk space
- **Recommended**: 32 GB RAM, 200 GB disk space
- **GPU**: Not required for data generation

### Software
- Python 3.8 or higher
- Internet connection (for RCSB API queries and PDB downloads)

## Quick Start (Using Pre-built Dataset)

If you just want to use the dataset without regenerating it:

```bash
# Clone the repository
git clone https://github.com/omagebright/MPLID.git
cd MPLID

# Install dependencies
pip install -r requirements.txt

# Load and explore
python -c "import pandas as pd; df = pd.read_csv('data/processed/train_residues.csv.gz'); print(df.head())"
```

## Full Reproduction Pipeline

### Step 1: Setup Environment

```bash
# Clone repository
git clone https://github.com/omagebright/MPLID.git
cd MPLID

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # Linux/Mac
# or: venv\Scripts\activate  # Windows

# Install dependencies
pip install -r requirements.txt
```

### Step 2: Run the Pipeline

The main pipeline script handles everything automatically:

```bash
# Full pipeline (downloads ~8,200 PDB structures, takes 4-8 hours)
python scripts/run_rcsb_only_pipeline.py

# Test run with limited proteins
python scripts/run_rcsb_only_pipeline.py --max-proteins 100
```

### Step 3: Verify Output

After completion, verify the generated files:

```bash
# Check output files
ls -la data/processed/

# Expected files:
# - train_residues.csv.gz (~41 MB)
# - val_residues.csv.gz (~18 MB)
# - test_residues.csv.gz (~16 MB)
# - protein_metadata.csv (~340 KB)

# Check statistics
cat data/statistics/pipeline_stats.json
```

## Pipeline Stages in Detail

### Stage 1: Query RCSB for Lipid-Containing Structures

```python
# Or run individually:
python scripts/01_query_rcsb_lipids.py
```

This queries RCSB Search API for structures containing any of 28 primary lipid codes (a focused subset used for efficient querying; the full 117-code recognition set is applied during contact calculation).

**Output**: `data/interim/rcsb_lipid_structures.json`

### Stage 2: Download PDB Structures

```python
python scripts/02_download_structures.py
```

Downloads original PDB files preserving HETATM records with lipid coordinates.

**Output**: `data/raw/pdb_original/*.pdb`

### Stage 3: Calculate Lipid Contacts

```python
python scripts/03_calculate_contacts.py
```

For each protein:
1. Parse structure using BioPython
2. Identify protein residues (standard amino acids)
3. Identify lipid molecules (117 recognized codes)
4. Calculate minimum heavy-atom distances
5. Label residues as contact (≤4.0 Å) or non-contact

**Output**: `data/interim/contacts/*.csv`

### Stage 4: Cluster Sequences

```python
python scripts/04_cluster_sequences.py
```

Clusters proteins by sequence identity using MMseqs2 (30% threshold).

**Output**: `data/interim/sequence_clusters.tsv`

### Stage 5: Create Train/Val/Test Splits

```python
python scripts/05_create_splits.py
```

Performs cluster-level stratified splitting:
- Train: 60% of clusters
- Validation: 20% of clusters  
- Test: 20% of clusters

**Output**: Final CSV files in `data/processed/`

## Expected Statistics

Your regenerated dataset should match these statistics (±1% due to PDB updates):

| Metric | Expected Value |
|--------|----------------|
| Total proteins | ~4,700 |
| Total residues | ~8,000,000 |
| Contact residues | ~80,000 |
| Contact rate | ~1.0% |
| Sequence clusters | ~800 |

## Troubleshooting

### Network Proxy Issues

If behind a corporate proxy:

```python
# Add to script or set environment variables
import os
os.environ['HTTP_PROXY'] = 'http://your.proxy:port'
os.environ['HTTPS_PROXY'] = 'http://your.proxy:port'
```

### Memory Issues

For systems with limited RAM:

```bash
# Process in smaller batches
python scripts/run_rcsb_only_pipeline.py --batch-size 100
```

### PDB Download Failures

Some structures may fail to download. The pipeline logs these and continues.
Typical failure rate: <1% of structures.

```bash
# Retry failed downloads
python scripts/02_download_structures.py --retry-failed
```

### MMseqs2 Not Found

If MMseqs2 is not installed, the pipeline uses CD-HIT as fallback:

```bash
# Install MMseqs2 (recommended)
conda install -c bioconda mmseqs2

# Or use CD-HIT fallback
pip install cd-hit
```

## Reproducibility Notes

### Random Seeds

All stochastic operations use fixed seeds for reproducibility:
- Dataset splitting: seed=42
- Stratification: seed=42

### Version Pinning

Key package versions used in original dataset generation:
- Python 3.12
- pandas 2.3.3
- numpy 2.3.4
- biopython 1.79+
- scikit-learn 1.7.2

### PDB Database Updates

The RCSB PDB database is continuously updated. Results may vary slightly
if structures are added/removed/updated between generation runs.

To reproduce exact results from the paper, use the archived dataset from Zenodo:
https://doi.org/10.5281/zenodo.18487585

## Output File Formats

### train/val/test_residues.csv.gz

| Column | Type | Description |
|--------|------|-------------|
| pdb_id | str | 4-character PDB ID |
| chain_id | str | Chain identifier |
| residue_number | int | Residue position |
| residue_name | str | 3-letter AA code |
| is_contact | bool | Contact label |
| min_distance | float | Min distance to lipid (Å) |
| cluster_id | int | Sequence cluster |
| split | str | train/val/test |

### protein_metadata.csv

| Column | Type | Description |
|--------|------|-------------|
| pdb_id | str | PDB identifier |
| resolution | float | Structure resolution (Å) |
| method | str | X-ray/cryo-EM/NMR |
| lipid_codes | str | Lipids found (comma-separated) |
| n_residues | int | Total residues |
| n_contacts | int | Contact residues |

## Citation

If you use this dataset or reproduction pipeline, please cite:

```bibtex
@article{omage2026mplid,
  author = {Omage, Folorunsho Bright and Mazoni, Ivan and Yano, Inácio Henrique and Neshich, Goran},
  title = {{MPLID}: A Large-Scale Dataset of Experimentally Validated Membrane Protein-Lipid Contact Residues for Machine Learning},
  journal = {GigaScience},
  year = {2026},
  doi = {10.5281/zenodo.18487585}
}
```

## Support

For issues or questions:
1. Open an issue on GitHub: https://github.com/omagebright/MPLID/issues
2. Email: omagefolorunsho@gmail.com

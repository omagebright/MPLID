# MPLID Data Dictionary

This document provides detailed descriptions of all data fields in the MPLID dataset.

## Residue-Level Data Files

Files: `train_residues.csv`, `val_residues.csv`, `test_residues.csv`

### Column Descriptions

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `pdb_id` | string | PDB structure identifier (4 characters) | `4HYC` |
| `chain_id` | string | Protein chain identifier | `A` |
| `residue_number` | integer | Residue position in the chain | `142` |
| `residue_name` | string | Three-letter amino acid code | `LEU` |
| `is_contact` | integer | Binary contact label (1=contact, 0=non-contact) | `1` |
| `lipid_type` | string | Lipid code if contact, null otherwise | `CLR` |
| `min_distance` | float | Minimum distance to nearest lipid (Å), null if no lipid | `3.45` |
| `cluster_id` | integer | Sequence cluster assignment | `127` |
| `split` | string | Dataset split assignment | `train` |

### Value Constraints

- `pdb_id`: Uppercase alphanumeric, exactly 4 characters
- `chain_id`: Single character (A-Z, 0-9)
- `residue_number`: Positive integer
- `residue_name`: Standard three-letter amino acid codes (see below)
- `is_contact`: Binary (0 or 1)
- `lipid_type`: Valid PDB lipid code or null
- `min_distance`: Positive float or null
- `cluster_id`: Positive integer
- `split`: One of `train`, `val`, `test`

### Standard Amino Acid Codes

| Code | Full Name |
|------|-----------|
| ALA | Alanine |
| ARG | Arginine |
| ASN | Asparagine |
| ASP | Aspartic acid |
| CYS | Cysteine |
| GLN | Glutamine |
| GLU | Glutamic acid |
| GLY | Glycine |
| HIS | Histidine |
| ILE | Isoleucine |
| LEU | Leucine |
| LYS | Lysine |
| MET | Methionine |
| PHE | Phenylalanine |
| PRO | Proline |
| SER | Serine |
| THR | Threonine |
| TRP | Tryptophan |
| TYR | Tyrosine |
| VAL | Valine |

## Protein Metadata

File: `protein_metadata.csv`

### Column Descriptions

| Column | Type | Description |
|--------|------|-------------|
| `pdb_id` | string | PDB structure identifier |
| `chain_id` | string | Primary chain identifier |
| `cluster_id` | integer | Sequence cluster assignment |
| `split` | string | Dataset split assignment |
| `n_residues` | integer | Total number of residues |
| `n_contacts` | integer | Number of lipid contact residues |
| `contact_rate` | float | Proportion of contact residues |

## Lipid Distribution

File: `data/statistics/lipid_distribution.csv`

### Column Descriptions

| Column | Type | Description |
|--------|------|-------------|
| `lipid_code` | string | PDB chemical component code |
| `lipid_name` | string | Full lipid name |
| `category` | string | Lipid category |
| `count` | integer | Number of contact residues |
| `n_proteins` | integer | Number of proteins with this lipid |

## Recognized Lipid Codes

### Phospholipids

| Code | Name |
|------|------|
| CDL | Cardiolipin |
| POV | 1-palmitoyl-2-oleoyl-phosphatidylcholine |
| PCW | Phosphatidylcholine |
| PEE | Phosphatidylethanolamine |
| PGV | Phosphatidylglycerol |
| PLM | Palmitic acid |
| PLC | Dipalmitoyl phosphatidylcholine |
| POP | 1-palmitoyl-2-oleoyl-phosphatidylethanolamine |
| PPE | Palmitoyl-2-oleoyl-sn-glycero-3-phosphoethanolamine |
| PGP | Phosphatidylglycerol phosphate |
| PIO | Phosphatidylinositol |
| PIP | Phosphatidylinositol-4-phosphate |
| EPE | 1,2-dioleoyl-sn-glycero-3-phosphoethanolamine |
| LPE | Lysophosphatidylethanolamine |

### Sphingolipids

| Code | Name |
|------|------|
| SPH | Sphingosine |
| S1P | Sphingosine-1-phosphate |
| HXJ | Hexosylceramide |

### Sterols

| Code | Name |
|------|------|
| CLR | Cholesterol |
| CHD | Cholesteryl hemisuccinate |
| Y01 | Cholesterol hemisuccinate |
| BCL | Beta-sitosterol |
| CHL | Chlorophyll |
| LHG | Ergosterol |

### Fatty Acids

| Code | Name |
|------|------|
| MYR | Myristic acid |
| OLA | Oleic acid |
| STE | Stearic acid |
| ARA | Arachidonic acid |
| DHA | Docosahexaenoic acid |
| MYS | Myristoyl |
| PAM | Palmitate |
| LNL | Linoleic acid |

### Detergents and Membrane Mimetics

| Code | Name |
|------|------|
| LDA | Lauryl dimethylamine-N-oxide |
| LMT | n-Dodecyl-beta-D-maltoside |
| BOG | Beta-octyl glucoside |
| OLC | n-Octyl-beta-D-glucopyranoside |
| DPC | n-Dodecylphosphocholine |
| DMU | n-Decyl-beta-D-maltoside |
| UNL | Unknown ligand (lipid-like) |

## Dataset Summary Statistics

File: `data/statistics/dataset_summary.json`

### Structure

```json
{
  "total_proteins": 2792,
  "total_residues": 4717703,
  "total_contacts": 39224,
  "total_clusters": 594,
  "contact_rate": 0.0083,
  "random_seed": 42,
  "splits": {
    "train": {
      "proteins": 1840,
      "residues": 2634209,
      "contacts": 25891,
      "clusters": 392
    },
    "val": {
      "proteins": 811,
      "residues": 1632603,
      "contacts": 10112,
      "clusters": 121
    },
    "test": {
      "proteins": 541,
      "residues": 867430,
      "contacts": 5221,
      "clusters": 81
    }
  }
}
```

## File Formats

### CSV Files

- Encoding: UTF-8
- Delimiter: Comma (,)
- Header: First row contains column names
- Null values: Empty string
- Quoting: Double quotes for strings containing commas

### JSON Files

- Encoding: UTF-8
- Formatting: Indented with 2 spaces
- Null values: JSON null

## Data Quality Notes

1. **Residue numbering**: Matches original PDB numbering, which may include insertion codes or gaps
2. **Missing lipids**: Proteins with no recognized lipids have all `is_contact = 0`
3. **Multiple chains**: Each chain is processed independently
4. **Cluster assignment**: Based on representative chain sequence

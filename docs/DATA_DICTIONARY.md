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
| `is_contact` | boolean | True if any heavy atom within 4.0 A of lipid heavy atom | `True` |
| `label_source` | string | Always "EXPERIMENTAL" (from crystallized lipids) | `EXPERIMENTAL` |
| `confidence` | string | Label confidence level | `high` |
| `min_distance` | float | Minimum heavy-atom-to-lipid distance (A) | `3.45` |
| `cluster_id` | integer | Sequence cluster assignment | `127` |
| `split` | string | Dataset split assignment | `train` |

### Value Constraints

- `pdb_id`: Uppercase alphanumeric, exactly 4 characters
- `chain_id`: Single character (A-Z, 0-9)
- `residue_number`: Positive integer
- `residue_name`: Standard three-letter amino acid codes (see below)
- `is_contact`: Boolean (True or False)
- `label_source`: Always "EXPERIMENTAL"
- `confidence`: Label confidence level
- `min_distance`: Positive float
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

## Recognized Lipid Codes (117 total)

### Phospholipids (24 codes)

| Code | Name |
|------|------|
| PC1 | 1,2-diacyl-sn-glycero-3-phosphocholine |
| PCW | 1,2-dioleoyl-sn-glycero-3-phosphocholine (DOPC) |
| POV | Palmitoyl-oleoyl-phosphatidylcholine (POPC) |
| PLC | Diundecyl phosphatidylcholine |
| 1PL | 1-palmitoyl-sn-glycero-3-phosphocholine |
| 2PL | 1,2-dipalmitoyl-sn-glycero-3-phosphocholine |
| XPC | Phosphatidylcholine variant |
| PEE | 1,2-dioleoyl-sn-glycero-3-phosphoethanolamine (DOPE) |
| PTY | Phosphatidylethanolamine |
| MPE | Myristic-palmitoyl PE |
| DPE | Dipalmitoyl PE |
| PGV | Phosphatidylglycerol |
| EPG | Egg PG |
| LHG | 1,2-dipalmitoyl-phosphatidylglycerol |
| PGR | PG variant |
| DPG | Diphosphatidylglycerol |
| LPS | Lysophosphatidylserine |
| DPS | Dipalmitoyl PS |
| PTI | Phosphatidylinositol |
| IPC | Inositol phosphoceramide |
| PIP | PI phosphate |
| PI3 | PI 3-phosphate |
| PI4 | PI 4-phosphate |
| P34 | PI 3,4-bisphosphate |

### Cardiolipin (4 codes)

| Code | Name |
|------|------|
| CDL | Cardiolipin |
| C9V | Cardiolipin-boron complex |
| 18W | Tetraoleoyl cardiolipin |
| LCL | Lyso-cardiolipin |

### Sphingolipids (7 codes)

| Code | Name |
|------|------|
| HXJ | Sphingomyelin/ceramide variant |
| SPH | Sphingosine |
| S1P | Sphingosine-1-phosphate |
| CRT | Ceramide |
| GCR | Glucosylceramide |
| GSP | Ganglioside |
| GLF | Galactosylceramide |

### Sterols (9 codes)

| Code | Name |
|------|------|
| CLR | Cholesterol |
| CHD | Cholesterol derivative |
| Y01 | Cholesterol hemisuccinate |
| OHC | 25-hydroxycholesterol |
| OCL | Cholesteryl oleate |
| LNS | Lanosterol |
| ERG | Ergosterol |
| SIT | Beta-sitosterol |
| STI | Stigmasterol |

### Fatty Acids (14 codes)

| Code | Name |
|------|------|
| PLM | Palmitic acid (C16:0) |
| MYR | Myristic acid (C14:0) |
| OLA | Oleic acid (C18:1) |
| STE | Stearic acid (C18:0) |
| ARA | Arachidonic acid (C20:4) |
| DHA | Docosahexaenoic acid (C22:6) |
| LNL | Linoleic acid (C18:2) |
| LNN | Alpha-linolenic acid (C18:3) |
| EPA | Eicosapentaenoic acid (C20:5) |
| PAM | Palmitate ion |
| DCR | Decanoic acid |
| UND | Undecanoic acid |
| LAU | Lauric acid (C12:0) |
| CAP | Capric acid (C10:0) |

### Glycerolipids (6 codes)

| Code | Name |
|------|------|
| TGL | Triglyceride |
| TAG | Triacylglycerol |
| DGA | Diglyceride |
| DAG | Diacylglycerol |
| MAG | Monoacylglycerol |
| GMO | Glyceryl monooleate |

### Detergents and Membrane Mimetics (22 codes)

| Code | Name |
|------|------|
| LDA | Lauryl dimethylamine-N-oxide (LDAO) |
| LMT | Dodecyl-beta-D-maltoside (DDM) |
| OLC | Monoolein (glyceryl monooleate) |
| BOG | n-Octyl-beta-D-glucopyranoside |
| BGC | Beta-D-glucopyranosyl |
| C8E | Octyl-PEG |
| C10 | Decyl maltoside |
| C12 | Dodecyl maltoside |
| UDM | n-Undecyl-beta-D-maltoside |
| NM | n-Nonyl-beta-D-maltoside |
| DM | n-Decyl-beta-D-maltoside |
| HTG | Heptyl thioglucoside |
| OG | Octyl glucoside |
| NG | Nonyl glucoside |
| SDS | Sodium dodecyl sulfate |
| DPC | Dodecylphosphocholine |
| FC6 | Fos-choline 6 |
| FC8 | Fos-choline 8 |
| FC10 | Fos-choline 10 |
| FC12 | Fos-choline 12 |
| FC14 | Fos-choline 14 |
| FC16 | Fos-choline 16 |

### Lipid A Components (3 codes)

| Code | Name |
|------|------|
| LPA | Lipid A |
| KDO | 3-deoxy-D-manno-octulosonic acid |
| 6LP | Lipid A disaccharide |

### CHARMM Simulation Nomenclature (27 codes)

These codes are commonly found in cryo-EM depositions where lipid bilayers are modeled using molecular dynamics force field naming conventions.

| Code | Name |
|------|------|
| POPC, POPE, POPS, POPG, POPA | Palmitoyl-oleoyl variants |
| DPPC, DPPE, DPPS, DPPG | Dipalmitoyl variants |
| DOPC, DOPE, DOPS, DOPG | Dioleoyl variants |
| DMPC, DMPE, DMPG | Dimyristoyl variants |
| DLPC, DLPE, DLPG | Dilauroyl variants |
| SOPC, SDPC, SLPC | Stearoyl-containing variants |
| PLPC, PAPC, SAPC | Mixed-chain variants |
| CHL1, CHOL | Cholesterol (CHARMM names) |

### Generic (1 code)

| Code | Name |
|------|------|
| LIP | Generic lipid |

## Dataset Summary Statistics

File: `data/statistics/dataset_summary.json`

### Structure

```json
{
  "total_proteins": 4704,
  "total_residues": 8055325,
  "total_contacts": 80439,
  "total_clusters": 813,
  "contact_rate": 0.0100,
  "random_seed": 42,
  "splits": {
    "train": {
      "proteins": 2578,
      "residues": 4907696,
      "contacts": 56976
    },
    "val": {
      "proteins": 1051,
      "residues": 1403838,
      "contacts": 9626
    },
    "test": {
      "proteins": 1075,
      "residues": 1743791,
      "contacts": 13837
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

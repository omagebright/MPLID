# MPLID Methodology

This document describes the methodology used to create the Membrane Protein-Lipid Interface Dataset (MPLID).

## Overview

MPLID provides experimentally-validated residue-level annotations of lipid contacts for membrane proteins. Unlike computational approaches that predict membrane boundaries, MPLID uses crystallized lipid molecules in PDB structures as the ground truth for lipid contact identification.

## Data Sources

### Primary Source

**RCSB Protein Data Bank (PDB)**
- Source of atomic coordinates for all proteins
- Provides crystallized lipid molecule positions (HETATM records)
- URL: https://www.rcsb.org/

### Data Collection

MPLID uses a single data source: the RCSB PDB. We query RCSB for all structures containing recognized lipid chemical components (117 codes), download the original PDB files preserving HETATM records, and compute contacts directly from atomic coordinates. No OPM database, no membrane plane calculations, no external computational dependencies.

## Contact Definition

### Distance Criterion

A residue is labeled as a **lipid contact** if:

$$d_{min}(R, L) \leq 4.0 \, \text{Å}$$

where:
- $R$ represents a protein residue
- $L$ represents a lipid molecule
- $d_{min}$ is the minimum distance between any heavy atom of $R$ and any heavy atom of $L$

### Heavy Atoms

Only heavy atoms (non-hydrogen atoms) are considered for distance calculations. This includes:
- **Protein**: C, N, O, S atoms in amino acid residues
- **Lipid**: All non-hydrogen atoms in lipid molecules

### Cutoff Selection

The 4.0 Å cutoff was chosen based on:
1. Standard van der Waals contact distances in protein-ligand interactions
2. Consistency with prior studies on membrane protein-lipid interfaces
3. Balance between sensitivity (detecting true contacts) and specificity (avoiding false positives)

## Lipid Recognition

### Supported Lipid Types

MPLID recognizes 117 lipid codes from the PDB Chemical Component Dictionary and CHARMM simulation nomenclature. Major categories include:

| Category | Count | Examples | Description |
|----------|-------|----------|-------------|
| Phospholipids | 24 | PC1, PCW, POV, PEE, PGV | Membrane bilayer components |
| Cardiolipin | 4 | CDL, C9V, 18W, LCL | Mitochondrial lipids |
| Sphingolipids | 7 | SPH, S1P, HXJ, CRT | Sphingosine-based lipids |
| Sterols | 9 | CLR, CHD, Y01, ERG | Cholesterol and derivatives |
| Fatty acids | 14 | PLM, MYR, OLA, STE, ARA | Free fatty acid chains |
| Glycerolipids | 6 | TGL, DAG, MAG, GMO | Glycerol-based lipids |
| Detergents | 22 | LDA, LMT, BOG, OLC, DPC | Membrane-mimetic molecules |
| Lipid A | 3 | LPA, KDO, 6LP | LPS components |
| CHARMM simulation | 27 | POPC, POPE, POPG | Cryo-EM deposition names |
| Generic | 1 | LIP | Unspecified lipid |

### Detergent Handling

Detergent molecules used in crystallization are included as they occupy similar positions to native lipids and provide valid contact information for machine learning applications.

## Quality Control

### Sequence Clustering

To prevent data leakage between dataset splits, proteins are clustered by sequence similarity:

1. **Method**: MMseqs2 (or CD-HIT as fallback)
2. **Identity threshold**: 30%
3. **Coverage threshold**: 80%

Proteins within the same cluster are always assigned to the same split (train/val/test).

### Structure Validation

Each structure undergoes validation:
1. Parse PDB file using BioPython
2. Verify presence of protein chains (standard amino acids)
3. Verify presence of recognized lipid molecules
4. Calculate contact distances for all protein-lipid pairs

### Data Integrity

- Residue numbering matches PDB SEQRES records
- Chain identifiers are preserved
- All calculations are reproducible with provided scripts

## Dataset Splitting

### Split Strategy

The dataset is split at the **protein level** using stratified sampling:

| Split | Proportion | Purpose |
|-------|------------|---------|
| Training | 60% | Model training |
| Validation | 20% | Hyperparameter tuning |
| Test | 20% | Final evaluation |

### Stratification

Clusters are stratified by size (number of residues) to ensure balanced representation across splits. This prevents the test set from being dominated by either very large or very small proteins.

### Reproducibility

All splits use a fixed random seed (42) for reproducibility. The exact split assignments are provided in the dataset files.

## Limitations

### Coverage Limitations

1. **Crystallization bias**: Only proteins that crystallize with lipids are included
2. **Lipid representation**: Native membrane lipids may differ from crystallized lipids
3. **Dynamic interactions**: Crystal structures capture static snapshots

### Label Noise Sources

1. **Partial lipid occupancy**: Some lipid positions may have low occupancy
2. **Crystal packing artifacts**: Lipids may occupy non-native positions
3. **Resolution effects**: Lower resolution structures may have less accurate lipid positions

## Comparison with Related Approaches

### vs. OPM-derived Labels

| Aspect | MPLID | OPM-derived |
|--------|-------|-------------|
| Label source | Crystallized lipids | Computed membrane plane |
| Precision | Residue-level (4.0 Å) | Zone-based (±15 Å) |
| Validation | Experimental | Computational |
| Coverage | 4,704 proteins | 15,096 proteins |

### vs. DREAMM Dataset

| Aspect | MPLID | DREAMM |
|--------|-------|--------|
| Proteins | 4,704 | 54 |
| Label source | Crystallized lipids | Literature curation (EPR, fluorescence, mutagenesis) |
| Availability | Public (Zenodo + GitHub) | Distributed within prediction tool |

## References

1. Lomize, M. A., et al. (2012). OPM database and PPM web server: resources for positioning of proteins in membranes. Nucleic Acids Research, 40(D1), D370-D376.

2. Berman, H. M., et al. (2000). The Protein Data Bank. Nucleic Acids Research, 28(1), 235-242.

3. Steinegger, M., & Söding, J. (2017). MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology, 35(11), 1026-1028.

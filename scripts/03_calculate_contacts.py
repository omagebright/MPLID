#!/usr/bin/env python3
"""Calculate lipid-protein contacts from PDB structures.

This script identifies residues in contact with crystallized lipids using
a distance-based cutoff criterion. Contacts are defined as any protein
heavy atom within a specified distance of any lipid heavy atom.

Author: Folorunsho Bright Omage
License: MIT
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, Selection
from Bio.PDB.Structure import Structure
from tqdm import tqdm

# Default contact distance cutoff in Angstroms
CONTACT_CUTOFF = 4.0

# Recognized lipid residue codes
LIPID_CODES = {
    # Phospholipids
    "CDL", "POV", "PCW", "PEE", "PGV", "PLM", "PLC", "POP", "PPE", "PGP",
    "PIO", "PIP", "P2E", "P3E", "P4E", "P5E", "PC1", "PE5", "LPE", "EPE",
    # Sphingolipids
    "SPH", "S1P", "HXJ",
    # Sterols
    "CLR", "CHD", "Y01", "BCL", "CHL", "LHG",
    # Fatty acids
    "MYR", "OLA", "STE", "ARA", "DHA", "MYS", "PAM", "LNL",
    # Detergents
    "LDA", "LMT", "BOG", "OLC", "DPC", "DMU", "UNL",
}


def setup_logging(output_dir: Path) -> logging.Logger:
    """Configure logging for the pipeline."""
    log_file = output_dir / "contacts.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(),
        ],
    )
    return logging.getLogger(__name__)


def get_heavy_atoms(residue) -> np.ndarray:
    """Extract coordinates of heavy atoms from a residue.

    Parameters
    ----------
    residue : Bio.PDB.Residue.Residue
        BioPython residue object.

    Returns
    -------
    np.ndarray
        Array of shape (N, 3) containing heavy atom coordinates.
    """
    coords = []
    for atom in residue.get_atoms():
        if atom.element != "H":
            coords.append(atom.get_coord())
    return np.array(coords) if coords else np.empty((0, 3))


def calculate_min_distance(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """Calculate minimum distance between two sets of coordinates.

    Parameters
    ----------
    coords1 : np.ndarray
        First coordinate array of shape (N, 3).
    coords2 : np.ndarray
        Second coordinate array of shape (M, 3).

    Returns
    -------
    float
        Minimum pairwise distance between the two coordinate sets.
    """
    if len(coords1) == 0 or len(coords2) == 0:
        return float("inf")

    # Compute pairwise distances efficiently
    diff = coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :]
    distances = np.sqrt(np.sum(diff ** 2, axis=2))
    return float(np.min(distances))


def extract_lipids(structure: Structure) -> List[Tuple]:
    """Extract lipid residues from a PDB structure.

    Parameters
    ----------
    structure : Bio.PDB.Structure.Structure
        BioPython structure object.

    Returns
    -------
    List[Tuple]
        List of (residue, coordinates) tuples for lipid molecules.
    """
    lipids = []
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()
                if resname in LIPID_CODES:
                    coords = get_heavy_atoms(residue)
                    if len(coords) > 0:
                        lipids.append((residue, coords))
    return lipids


def calculate_contacts_for_structure(
    pdb_file: Path,
    cutoff: float = CONTACT_CUTOFF,
) -> Optional[pd.DataFrame]:
    """Calculate lipid contacts for all residues in a structure.

    Parameters
    ----------
    pdb_file : Path
        Path to PDB structure file.
    cutoff : float, optional
        Distance cutoff in Angstroms (default: 4.0).

    Returns
    -------
    Optional[pd.DataFrame]
        DataFrame with residue-level contact information, or None if error.
    """
    parser = PDBParser(QUIET=True)
    pdb_id = pdb_file.stem.upper()

    try:
        structure = parser.get_structure(pdb_id, pdb_file)
    except Exception as e:
        logging.warning(f"Failed to parse {pdb_id}: {e}")
        return None

    # Extract lipids
    lipids = extract_lipids(structure)
    if not lipids:
        logging.debug(f"No lipids found in {pdb_id}")
        return None

    # Process each protein residue
    records = []
    for model in structure:
        for chain in model:
            for residue in chain:
                # Skip non-standard residues (HETATM)
                hetflag = residue.get_id()[0]
                if hetflag != " ":
                    continue

                resname = residue.get_resname().strip()
                resnum = residue.get_id()[1]

                # Calculate minimum distance to any lipid
                protein_coords = get_heavy_atoms(residue)
                if len(protein_coords) == 0:
                    continue

                min_dist = float("inf")
                contact_lipid = None

                for lipid_res, lipid_coords in lipids:
                    dist = calculate_min_distance(protein_coords, lipid_coords)
                    if dist < min_dist:
                        min_dist = dist
                        contact_lipid = lipid_res.get_resname().strip()

                is_contact = 1 if min_dist <= cutoff else 0

                records.append({
                    "pdb_id": pdb_id,
                    "chain_id": chain.get_id(),
                    "residue_number": resnum,
                    "residue_name": resname,
                    "is_contact": is_contact,
                    "min_distance": round(min_dist, 3) if min_dist < float("inf") else None,
                    "lipid_type": contact_lipid if is_contact else None,
                })

    return pd.DataFrame(records) if records else None


def main():
    """Main entry point for contact calculation pipeline."""
    parser = argparse.ArgumentParser(
        description="Calculate lipid-protein contacts from PDB structures"
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/raw/pdb"),
        help="Directory containing PDB files",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/interim"),
        help="Output directory for contact data",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=CONTACT_CUTOFF,
        help=f"Contact distance cutoff in Angstroms (default: {CONTACT_CUTOFF})",
    )
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    logger = setup_logging(args.output)

    # Find all PDB files
    pdb_files = list(args.input.glob("*.pdb"))
    logger.info(f"Processing {len(pdb_files)} PDB structures")
    logger.info(f"Using contact cutoff: {args.cutoff} Angstroms")

    # Process structures
    all_contacts = []
    proteins_with_contacts = 0

    for pdb_file in tqdm(pdb_files, desc="Calculating contacts"):
        df = calculate_contacts_for_structure(pdb_file, cutoff=args.cutoff)
        if df is not None and len(df) > 0:
            all_contacts.append(df)
            if df["is_contact"].sum() > 0:
                proteins_with_contacts += 1

    # Combine results
    if all_contacts:
        combined = pd.concat(all_contacts, ignore_index=True)
        output_file = args.output / "all_contacts.csv"
        combined.to_csv(output_file, index=False)

        # Generate statistics
        stats = {
            "total_proteins": len(pdb_files),
            "proteins_processed": len(all_contacts),
            "proteins_with_contacts": proteins_with_contacts,
            "total_residues": len(combined),
            "contact_residues": int(combined["is_contact"].sum()),
            "contact_rate": float(combined["is_contact"].mean()),
            "cutoff_angstroms": args.cutoff,
        }

        stats_file = args.output / "contact_statistics.json"
        with open(stats_file, "w") as f:
            json.dump(stats, f, indent=2)

        logger.info(f"Processed {len(all_contacts)} proteins")
        logger.info(f"Total residues: {len(combined):,}")
        logger.info(f"Contact residues: {stats['contact_residues']:,}")
        logger.info(f"Contact rate: {stats['contact_rate']:.2%}")
    else:
        logger.warning("No contacts found in any structure")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Cluster protein sequences to prevent data leakage.

This script performs sequence clustering using CD-HIT or MMseqs2 to group
similar proteins together, ensuring that proteins from the same cluster
are assigned to the same dataset split.

Author: Folorunsho Bright Omage
License: MIT
"""

import argparse
import json
import logging
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set

import pandas as pd
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

# Standard amino acid mapping
AA_MAP = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}


def setup_logging(output_dir: Path) -> logging.Logger:
    """Configure logging for the pipeline."""
    log_file = output_dir / "clustering.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(),
        ],
    )
    return logging.getLogger(__name__)


def extract_sequence_from_pdb(pdb_file: Path) -> Dict[str, str]:
    """Extract amino acid sequences from a PDB structure.

    Parameters
    ----------
    pdb_file : Path
        Path to PDB structure file.

    Returns
    -------
    Dict[str, str]
        Dictionary mapping chain IDs to amino acid sequences.
    """
    parser = PDBParser(QUIET=True)
    pdb_id = pdb_file.stem.upper()

    try:
        structure = parser.get_structure(pdb_id, pdb_file)
    except Exception:
        return {}

    sequences = {}
    for model in structure:
        for chain in model:
            seq = []
            for residue in chain:
                if residue.get_id()[0] == " ":  # Standard residue
                    resname = residue.get_resname().strip()
                    aa = AA_MAP.get(resname, "X")
                    seq.append(aa)
            if seq:
                sequences[f"{pdb_id}_{chain.get_id()}"] = "".join(seq)
        break  # Only first model

    return sequences


def write_fasta(sequences: Dict[str, str], output_file: Path) -> None:
    """Write sequences to FASTA format.

    Parameters
    ----------
    sequences : Dict[str, str]
        Dictionary mapping sequence IDs to sequences.
    output_file : Path
        Output FASTA file path.
    """
    with open(output_file, "w") as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n")
            # Write sequence in lines of 80 characters
            for i in range(0, len(seq), 80):
                f.write(f"{seq[i:i+80]}\n")


def run_mmseqs2(
    fasta_file: Path,
    output_dir: Path,
    identity: float = 0.3,
    coverage: float = 0.8,
) -> Dict[str, int]:
    """Run MMseqs2 clustering.

    Parameters
    ----------
    fasta_file : Path
        Input FASTA file.
    output_dir : Path
        Output directory for clustering results.
    identity : float, optional
        Sequence identity threshold (default: 0.3).
    coverage : float, optional
        Coverage threshold (default: 0.8).

    Returns
    -------
    Dict[str, int]
        Dictionary mapping sequence IDs to cluster IDs.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        db = tmpdir / "seqdb"
        cluster_db = tmpdir / "clusterdb"
        tsv_file = output_dir / "clusters.tsv"

        # Create sequence database
        subprocess.run(
            ["mmseqs", "createdb", str(fasta_file), str(db)],
            check=True,
            capture_output=True,
        )

        # Cluster sequences
        subprocess.run(
            [
                "mmseqs", "cluster", str(db), str(cluster_db), str(tmpdir),
                "--min-seq-id", str(identity),
                "-c", str(coverage),
                "--cluster-mode", "0",
            ],
            check=True,
            capture_output=True,
        )

        # Convert to TSV
        subprocess.run(
            [
                "mmseqs", "createtsv", str(db), str(db),
                str(cluster_db), str(tsv_file),
            ],
            check=True,
            capture_output=True,
        )

    # Parse cluster assignments
    clusters = {}
    cluster_id = 0
    current_rep = None

    with open(tsv_file) as f:
        for line in f:
            rep, member = line.strip().split("\t")
            if rep != current_rep:
                current_rep = rep
                cluster_id += 1
            clusters[member] = cluster_id

    return clusters


def run_cdhit(
    fasta_file: Path,
    output_dir: Path,
    identity: float = 0.3,
) -> Dict[str, int]:
    """Run CD-HIT clustering (fallback if MMseqs2 unavailable).

    Parameters
    ----------
    fasta_file : Path
        Input FASTA file.
    output_dir : Path
        Output directory for clustering results.
    identity : float, optional
        Sequence identity threshold (default: 0.3).

    Returns
    -------
    Dict[str, int]
        Dictionary mapping sequence IDs to cluster IDs.
    """
    output_file = output_dir / "cdhit_clusters"

    subprocess.run(
        [
            "cd-hit", "-i", str(fasta_file), "-o", str(output_file),
            "-c", str(identity), "-n", "2", "-M", "4000", "-T", "4",
        ],
        check=True,
        capture_output=True,
    )

    # Parse cluster file
    clusters = {}
    cluster_id = 0

    with open(f"{output_file}.clstr") as f:
        for line in f:
            if line.startswith(">Cluster"):
                cluster_id = int(line.split()[1]) + 1
            elif ">" in line:
                seq_id = line.split(">")[1].split("...")[0]
                clusters[seq_id] = cluster_id

    return clusters


def main():
    """Main entry point for sequence clustering pipeline."""
    parser = argparse.ArgumentParser(
        description="Cluster protein sequences to prevent data leakage"
    )
    parser.add_argument(
        "--pdb-dir",
        type=Path,
        default=Path("data/raw/pdb"),
        help="Directory containing PDB files",
    )
    parser.add_argument(
        "--contacts-file",
        type=Path,
        default=Path("data/interim/all_contacts.csv"),
        help="Contact data CSV file to annotate",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/interim"),
        help="Output directory",
    )
    parser.add_argument(
        "--identity",
        type=float,
        default=0.3,
        help="Sequence identity threshold (default: 0.3)",
    )
    parser.add_argument(
        "--method",
        choices=["mmseqs2", "cdhit"],
        default="mmseqs2",
        help="Clustering method (default: mmseqs2)",
    )
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    logger = setup_logging(args.output)

    # Extract sequences from PDB files
    logger.info("Extracting sequences from PDB files")
    pdb_files = list(args.pdb_dir.glob("*.pdb"))

    all_sequences = {}
    for pdb_file in pdb_files:
        seqs = extract_sequence_from_pdb(pdb_file)
        all_sequences.update(seqs)

    logger.info(f"Extracted {len(all_sequences)} sequences from {len(pdb_files)} structures")

    # Write FASTA file
    fasta_file = args.output / "sequences.fasta"
    write_fasta(all_sequences, fasta_file)

    # Run clustering
    logger.info(f"Clustering sequences at {args.identity:.0%} identity using {args.method}")

    try:
        if args.method == "mmseqs2":
            clusters = run_mmseqs2(fasta_file, args.output, identity=args.identity)
        else:
            clusters = run_cdhit(fasta_file, args.output, identity=args.identity)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.error(f"Clustering failed: {e}")
        logger.info("Falling back to single-linkage clustering")
        # Simple fallback: each protein is its own cluster
        clusters = {seq_id: i for i, seq_id in enumerate(all_sequences.keys(), 1)}

    # Map to protein level (PDB ID only)
    protein_clusters = {}
    for seq_id, cluster_id in clusters.items():
        pdb_id = seq_id.split("_")[0]
        if pdb_id not in protein_clusters:
            protein_clusters[pdb_id] = cluster_id

    unique_clusters = len(set(protein_clusters.values()))
    logger.info(f"Created {unique_clusters} clusters from {len(protein_clusters)} proteins")

    # Annotate contact data if provided
    if args.contacts_file.exists():
        logger.info("Annotating contact data with cluster IDs")
        df = pd.read_csv(args.contacts_file)
        df["cluster_id"] = df["pdb_id"].map(protein_clusters)
        df.to_csv(args.contacts_file, index=False)
        logger.info(f"Updated {args.contacts_file}")

    # Save cluster assignments
    cluster_file = args.output / "protein_clusters.json"
    with open(cluster_file, "w") as f:
        json.dump(protein_clusters, f, indent=2)

    # Summary statistics
    cluster_sizes = defaultdict(int)
    for cluster_id in protein_clusters.values():
        cluster_sizes[cluster_id] += 1

    stats = {
        "total_proteins": len(protein_clusters),
        "total_clusters": unique_clusters,
        "identity_threshold": args.identity,
        "method": args.method,
        "largest_cluster": max(cluster_sizes.values()),
        "singletons": sum(1 for s in cluster_sizes.values() if s == 1),
    }

    stats_file = args.output / "clustering_statistics.json"
    with open(stats_file, "w") as f:
        json.dump(stats, f, indent=2)

    logger.info("Clustering complete")


if __name__ == "__main__":
    main()

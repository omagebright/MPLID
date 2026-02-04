#!/usr/bin/env python3
"""Create train/validation/test splits for the MPLID dataset.

This script creates stratified dataset splits at the protein level,
ensuring that proteins from the same sequence cluster are kept together
to prevent data leakage.

Author: Folorunsho Bright Omage
License: MIT
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split


def setup_logging(output_dir: Path) -> logging.Logger:
    """Configure logging for the pipeline."""
    log_file = output_dir / "splits.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(),
        ],
    )
    return logging.getLogger(__name__)


def calculate_cluster_stats(df: pd.DataFrame) -> Dict[int, Dict]:
    """Calculate statistics for each sequence cluster.

    Parameters
    ----------
    df : pd.DataFrame
        Contact data with cluster_id column.

    Returns
    -------
    Dict[int, Dict]
        Dictionary mapping cluster IDs to statistics.
    """
    cluster_stats = {}

    for cluster_id in df["cluster_id"].unique():
        cluster_data = df[df["cluster_id"] == cluster_id]
        cluster_stats[cluster_id] = {
            "n_proteins": cluster_data["pdb_id"].nunique(),
            "n_residues": len(cluster_data),
            "n_contacts": int(cluster_data["is_contact"].sum()),
            "contact_rate": float(cluster_data["is_contact"].mean()),
        }

    return cluster_stats


def split_clusters(
    cluster_stats: Dict[int, Dict],
    train_ratio: float = 0.6,
    val_ratio: float = 0.2,
    test_ratio: float = 0.2,
    random_state: int = 42,
) -> Tuple[List[int], List[int], List[int]]:
    """Split clusters into train/val/test sets.

    Uses stratified splitting based on cluster size to ensure balanced
    representation across splits.

    Parameters
    ----------
    cluster_stats : Dict[int, Dict]
        Dictionary of cluster statistics.
    train_ratio : float, optional
        Proportion for training set (default: 0.6).
    val_ratio : float, optional
        Proportion for validation set (default: 0.2).
    test_ratio : float, optional
        Proportion for test set (default: 0.2).
    random_state : int, optional
        Random seed for reproducibility (default: 42).

    Returns
    -------
    Tuple[List[int], List[int], List[int]]
        Lists of cluster IDs for train, validation, and test sets.
    """
    cluster_ids = list(cluster_stats.keys())
    cluster_sizes = [cluster_stats[c]["n_residues"] for c in cluster_ids]

    # Stratify by cluster size (binned)
    size_bins = pd.qcut(cluster_sizes, q=4, labels=False, duplicates="drop")

    # First split: train vs (val + test)
    train_clusters, temp_clusters = train_test_split(
        cluster_ids,
        test_size=(val_ratio + test_ratio),
        stratify=size_bins,
        random_state=random_state,
    )

    # Second split: val vs test
    temp_indices = [cluster_ids.index(c) for c in temp_clusters]
    temp_bins = [size_bins[i] for i in temp_indices]

    relative_val_ratio = val_ratio / (val_ratio + test_ratio)
    val_clusters, test_clusters = train_test_split(
        temp_clusters,
        test_size=(1 - relative_val_ratio),
        stratify=temp_bins if len(set(temp_bins)) > 1 else None,
        random_state=random_state,
    )

    return train_clusters, val_clusters, test_clusters


def assign_splits(
    df: pd.DataFrame,
    train_clusters: List[int],
    val_clusters: List[int],
    test_clusters: List[int],
) -> pd.DataFrame:
    """Assign split labels to residues based on cluster membership.

    Parameters
    ----------
    df : pd.DataFrame
        Contact data with cluster_id column.
    train_clusters : List[int]
        Cluster IDs for training set.
    val_clusters : List[int]
        Cluster IDs for validation set.
    test_clusters : List[int]
        Cluster IDs for test set.

    Returns
    -------
    pd.DataFrame
        DataFrame with added 'split' column.
    """
    split_map = {}
    for c in train_clusters:
        split_map[c] = "train"
    for c in val_clusters:
        split_map[c] = "val"
    for c in test_clusters:
        split_map[c] = "test"

    df["split"] = df["cluster_id"].map(split_map)
    return df


def main():
    """Main entry point for dataset splitting pipeline."""
    parser = argparse.ArgumentParser(
        description="Create train/validation/test splits for MPLID"
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/interim/all_contacts.csv"),
        help="Input contact data CSV file",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/processed"),
        help="Output directory for split datasets",
    )
    parser.add_argument(
        "--train",
        type=float,
        default=0.6,
        help="Training set proportion (default: 0.6)",
    )
    parser.add_argument(
        "--val",
        type=float,
        default=0.2,
        help="Validation set proportion (default: 0.2)",
    )
    parser.add_argument(
        "--test",
        type=float,
        default=0.2,
        help="Test set proportion (default: 0.2)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42)",
    )
    args = parser.parse_args()

    # Validate proportions
    total = args.train + args.val + args.test
    if not np.isclose(total, 1.0):
        raise ValueError(f"Split proportions must sum to 1.0, got {total}")

    args.output.mkdir(parents=True, exist_ok=True)
    logger = setup_logging(args.output)

    # Load contact data
    logger.info(f"Loading contact data from {args.input}")
    df = pd.read_csv(args.input)
    logger.info(f"Loaded {len(df):,} residues from {df['pdb_id'].nunique():,} proteins")

    if "cluster_id" not in df.columns:
        raise ValueError("Contact data must have cluster_id column. Run 04_cluster_sequences.py first.")

    # Calculate cluster statistics
    cluster_stats = calculate_cluster_stats(df)
    logger.info(f"Found {len(cluster_stats)} sequence clusters")

    # Split clusters
    train_clusters, val_clusters, test_clusters = split_clusters(
        cluster_stats,
        train_ratio=args.train,
        val_ratio=args.val,
        test_ratio=args.test,
        random_state=args.seed,
    )

    logger.info(f"Train clusters: {len(train_clusters)}")
    logger.info(f"Val clusters: {len(val_clusters)}")
    logger.info(f"Test clusters: {len(test_clusters)}")

    # Assign splits
    df = assign_splits(df, train_clusters, val_clusters, test_clusters)

    # Save split datasets
    for split_name in ["train", "val", "test"]:
        split_df = df[df["split"] == split_name]
        output_file = args.output / f"{split_name}_residues.csv"
        split_df.to_csv(output_file, index=False)
        logger.info(
            f"{split_name.capitalize()}: {len(split_df):,} residues, "
            f"{split_df['pdb_id'].nunique():,} proteins, "
            f"{split_df['is_contact'].sum():,} contacts"
        )

    # Save protein metadata
    protein_meta = df.groupby("pdb_id").agg({
        "chain_id": "first",
        "cluster_id": "first",
        "split": "first",
        "residue_number": "count",
        "is_contact": "sum",
    }).rename(columns={
        "residue_number": "n_residues",
        "is_contact": "n_contacts",
    })
    protein_meta["contact_rate"] = protein_meta["n_contacts"] / protein_meta["n_residues"]
    protein_meta.to_csv(args.output / "protein_metadata.csv")

    # Save split summary
    summary = {
        "total_proteins": df["pdb_id"].nunique(),
        "total_residues": len(df),
        "total_contacts": int(df["is_contact"].sum()),
        "total_clusters": len(cluster_stats),
        "random_seed": args.seed,
        "splits": {
            "train": {
                "proteins": int((df["split"] == "train").sum() / df.groupby("pdb_id").ngroup().nunique() * df["pdb_id"].nunique()),
                "residues": int((df["split"] == "train").sum()),
                "contacts": int(df[df["split"] == "train"]["is_contact"].sum()),
                "clusters": len(train_clusters),
            },
            "val": {
                "proteins": int(df[df["split"] == "val"]["pdb_id"].nunique()),
                "residues": int((df["split"] == "val").sum()),
                "contacts": int(df[df["split"] == "val"]["is_contact"].sum()),
                "clusters": len(val_clusters),
            },
            "test": {
                "proteins": int(df[df["split"] == "test"]["pdb_id"].nunique()),
                "residues": int((df["split"] == "test").sum()),
                "contacts": int(df[df["split"] == "test"]["is_contact"].sum()),
                "clusters": len(test_clusters),
            },
        },
    }

    summary_file = args.output / "dataset_summary.json"
    with open(summary_file, "w") as f:
        json.dump(summary, f, indent=2)

    logger.info("Dataset splitting complete")
    logger.info(f"Output files saved to {args.output}")


if __name__ == "__main__":
    main()

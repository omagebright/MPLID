#!/usr/bin/env python3
"""Run MPLID pipeline using ONLY RCSB crystallographic lipid contacts.

This script uses 100% experimental data with zero dependency on OPM:
1. Query RCSB to find ALL PDB structures with lipid ligands (~8,221 structures)
2. Download original PDB files directly from RCSB
3. Extract experimental lipid contacts using 4.0Å Cα-to-lipid cutoff
4. Cluster by sequence identity and create train/val/test splits

Key difference from run_experimental_only_pipeline.py:
- NO OPM intersection filtering (captures 100% of RCSB lipid structures)
- NO OPM file downloads needed
- Purely crystallographic ground truth

Dataset framing:
- "First large-scale predictor trained on crystallographically-observed lipid
  contacts rather than literature-curated functional data (DREAMM's 54 proteins)"
- Objective, reproducible labels from PDB HETATM records
- Scale: ~7,000+ proteins vs DREAMM's 54

Usage:
    python scripts/run_rcsb_only_pipeline.py
    python scripts/run_rcsb_only_pipeline.py --max-proteins 100
    python scripts/run_rcsb_only_pipeline.py --skip-rcsb-query  # Use cached list
"""

import argparse
import logging
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Set, List, Dict, Tuple, Optional
from dataclasses import dataclass, field

import requests
import pandas as pd
import numpy as np
from tqdm import tqdm

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.data.config import DataConfig
from src.data.lipid_contacts import LipidContactExtractor, LIPID_RESIDUE_NAMES
from src.data.sequence_clustering import SequenceClusterer, stratified_cluster_split

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

# Proxy for network requests (required at Embrapa)
PROXY = {
    "http": "http://proxy.cnptia.embrapa.br:3128",
    "https": "http://proxy.cnptia.embrapa.br:3128"
}


def query_rcsb_for_lipid_structures() -> Dict[str, List[str]]:
    """Query RCSB to find ALL PDB structures containing lipid ligands.

    Uses RCSB Search API to find structures with known lipid chemical components.
    Returns both the PDB IDs and which lipid codes were found in each.

    Returns:
        Dict mapping PDB IDs (uppercase) to list of lipid codes found
    """
    logger.info("Querying RCSB for structures with lipid ligands...")
    logger.info(f"Searching for {len(LIPID_RESIDUE_NAMES)} lipid codes from our curated list")

    # Use all lipid codes from our comprehensive list
    # Filter to 3-letter codes (RCSB uses 3-letter component IDs)
    lipid_codes = [code for code in LIPID_RESIDUE_NAMES if len(code) <= 3]

    # Add common 4-letter codes that we should query individually
    # (4-letter codes like POPC, DPPC are CHARMM names, less common in PDB)

    pdb_to_lipids: Dict[str, List[str]] = {}

    # RCSB Search API
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"

    # Create a session for connection pooling
    session = requests.Session()
    session.proxies = PROXY

    # Query for each lipid code
    for lipid_code in tqdm(lipid_codes, desc="Querying RCSB"):
        query = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
                    "operator": "exact_match",
                    "value": lipid_code
                }
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {
                    "start": 0,
                    "rows": 50000  # Get all results
                }
            }
        }

        try:
            response = requests.post(search_url, json=query, proxies=PROXY, timeout=60)
            if response.status_code == 200:
                data = response.json()
                count = data.get("total_count", 0)
                if "result_set" in data and count > 0:
                    for result in data["result_set"]:
                        pdb_id = result["identifier"].upper()
                        if pdb_id not in pdb_to_lipids:
                            pdb_to_lipids[pdb_id] = []
                        pdb_to_lipids[pdb_id].append(lipid_code)
                    logger.debug(f"  {lipid_code}: {count} structures")
            else:
                logger.warning(f"RCSB returned status {response.status_code} for {lipid_code}")
        except Exception as e:
            logger.warning(f"Failed to query for {lipid_code}: {e}")

    logger.info(f"Found {len(pdb_to_lipids)} unique PDB structures with lipid ligands")

    # Log lipid type distribution
    all_lipids = [lip for lips in pdb_to_lipids.values() for lip in lips]
    from collections import Counter
    lipid_counts = Counter(all_lipids)
    top_10 = lipid_counts.most_common(10)
    logger.info("Top 10 lipid codes by structure count:")
    for code, count in top_10:
        logger.info(f"  {code}: {count} structures")

    return pdb_to_lipids


@dataclass
class RCSBPipelineStats:
    """Statistics from RCSB-only pipeline."""
    total_rcsb_lipid_structures: int = 0
    proteins_processed: int = 0
    proteins_with_contacts: int = 0
    proteins_failed: int = 0
    total_residues: int = 0
    positive_residues: int = 0
    n_clusters: int = 0
    train_proteins: int = 0
    train_residues: int = 0
    train_positives: int = 0
    val_proteins: int = 0
    val_residues: int = 0
    val_positives: int = 0
    test_proteins: int = 0
    test_residues: int = 0
    test_positives: int = 0
    lipid_types_found: Dict[str, int] = field(default_factory=dict)
    processing_errors: List[str] = field(default_factory=list)


def download_pdb_structure(pdb_id: str, output_dir: Path, session: requests.Session) -> Optional[Path]:
    """Download original PDB structure from RCSB.

    Args:
        pdb_id: 4-character PDB ID
        output_dir: Directory to save the file
        session: Requests session with proxy configured

    Returns:
        Path to downloaded file, or None if failed
    """
    pdb_path = output_dir / f"{pdb_id.lower()}.pdb"

    # Skip if already downloaded
    if pdb_path.exists() and pdb_path.stat().st_size > 0:
        return pdb_path

    # Try multiple RCSB endpoints
    urls = [
        f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb",
        f"https://files.rcsb.org/view/{pdb_id.upper()}.pdb",
    ]

    for url in urls:
        try:
            response = session.get(url, timeout=60)
            if response.status_code == 200 and len(response.text) > 100:
                pdb_path.write_text(response.text)
                return pdb_path
        except Exception as e:
            continue

    return None


def extract_sequence(pdb_path: Path) -> Optional[str]:
    """Extract protein sequence from PDB file."""
    from Bio.PDB import PDBParser
    from Bio.SeqUtils import seq1

    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", str(pdb_path))

        sequences = []
        for model in structure:
            for chain in model:
                seq = ""
                for residue in chain:
                    if residue.id[0] == " ":  # Standard residue
                        try:
                            seq += seq1(residue.get_resname())
                        except:
                            pass
                if len(seq) >= 20:  # Minimum length for meaningful sequence
                    sequences.append(seq)

        # Return longest sequence (main chain)
        if sequences:
            return max(sequences, key=len)
        return None
    except Exception as e:
        return None


def process_protein(
    pdb_id: str,
    pdb_dir: Path,
    session: requests.Session,
    lipid_extractor: LipidContactExtractor
) -> Tuple[List[Dict], Optional[str], Optional[str]]:
    """Process a single protein for experimental lipid contacts.

    Args:
        pdb_id: PDB identifier
        pdb_dir: Directory for PDB files
        session: Requests session with proxy
        lipid_extractor: LipidContactExtractor instance

    Returns:
        Tuple of (residue dicts, sequence string, error message if any)
    """
    residues = []
    sequence = None
    error = None

    try:
        # Download PDB structure
        pdb_path = download_pdb_structure(pdb_id, pdb_dir, session)

        if not pdb_path or not pdb_path.exists():
            return [], None, f"Failed to download PDB {pdb_id}"

        # Extract contacts using EXPERIMENTAL lipid detection only
        # use_ca_only=True: Use Cα atom for residue representation
        # use_opm_membrane_plane=False: Never fall back to OPM Z-coordinates
        contacts = lipid_extractor.extract_contacts_from_pdb(
            pdb_path,
            use_ca_only=True,
            use_opm_membrane_plane=False
        )

        if not contacts:
            return [], None, f"No residues extracted from {pdb_id}"

        for contact in contacts:
            residues.append({
                "pdb_id": pdb_id,
                "chain_id": contact.chain_id,
                "residue_number": contact.residue_number,
                "residue_name": contact.residue_name,
                "is_contact": contact.is_contact,
                "label_source": "CRYSTALLOGRAPHIC",
                "confidence": "high",
                "min_distance": contact.min_distance
            })

        # Extract sequence from PDB
        sequence = extract_sequence(pdb_path)

    except Exception as e:
        error = f"Error processing {pdb_id}: {str(e)}"

    return residues, sequence, error


def main():
    parser = argparse.ArgumentParser(
        description="Run MPLID pipeline using ONLY RCSB crystallographic lipid contacts (no OPM dependency)"
    )
    parser.add_argument(
        "--max-proteins", type=int, default=None,
        help="Limit number of proteins to process (for testing)"
    )
    parser.add_argument(
        "--skip-rcsb-query", action="store_true",
        help="Skip RCSB query, use cached lipid structure list"
    )
    parser.add_argument(
        "--resume", action="store_true",
        help="Resume from previous run, skip already processed proteins"
    )
    args = parser.parse_args()

    config = DataConfig()
    stats = RCSBPipelineStats()

    # Setup logging to file
    log_dir = config.project_root / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"rcsb_only_pipeline_{timestamp}.log"

    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logging.getLogger().addHandler(file_handler)

    logger.info("=" * 70)
    logger.info("MPLID RCSB-Only Pipeline (100% Crystallographic Labels)")
    logger.info("=" * 70)
    logger.info("NO OPM dependency - using all RCSB structures with lipid ligands")
    logger.info("")

    # Step 1: Query RCSB for ALL lipid structures (no OPM filtering!)
    cache_file = config.interim_data_dir / "rcsb_lipid_structures_full.json"
    config.interim_data_dir.mkdir(parents=True, exist_ok=True)

    # Check for cached results (supports both old list format and new dict format)
    old_cache_file = config.interim_data_dir / "rcsb_lipid_structures.json"

    if args.skip_rcsb_query and (cache_file.exists() or old_cache_file.exists()):
        # Try new cache first, then old cache
        actual_cache = cache_file if cache_file.exists() else old_cache_file
        logger.info(f"Loading cached lipid structure list from {actual_cache}")
        with open(actual_cache) as f:
            cached_data = json.load(f)

        # Handle both formats: old (list) or new (dict)
        if isinstance(cached_data, list):
            # Old format: list of PDB IDs
            pdb_to_lipids = {pdb_id.upper(): ["UNKNOWN"] for pdb_id in cached_data}
            logger.info(f"Converted old cache format: {len(pdb_to_lipids)} structures")
        else:
            # New format: dict mapping PDB ID to lipid codes
            pdb_to_lipids = cached_data
    else:
        pdb_to_lipids = query_rcsb_for_lipid_structures()
        # Cache for future use
        with open(cache_file, "w") as f:
            json.dump(pdb_to_lipids, f, indent=2)
        logger.info(f"Cached lipid structure list to {cache_file}")

    proteins_to_process = list(pdb_to_lipids.keys())
    stats.total_rcsb_lipid_structures = len(proteins_to_process)

    logger.info(f"\n{'='*70}")
    logger.info(f"RCSB structures with lipid ligands: {stats.total_rcsb_lipid_structures:,}")
    logger.info(f"(100% of available data - NO OPM filtering)")
    logger.info(f"{'='*70}\n")

    if stats.total_rcsb_lipid_structures == 0:
        logger.error("No proteins found with lipids!")
        return

    # Limit proteins if requested (for testing)
    if args.max_proteins:
        proteins_to_process = proteins_to_process[:args.max_proteins]
        logger.info(f"Limited to {len(proteins_to_process)} proteins for testing")

    # Setup directories
    pdb_dir = config.raw_data_dir / "pdb_original"
    pdb_dir.mkdir(parents=True, exist_ok=True)

    # Resume support
    processed_file = config.interim_data_dir / "rcsb_processed_proteins.json"
    already_processed = set()
    existing_residues = []

    if args.resume and processed_file.exists():
        with open(processed_file) as f:
            resume_data = json.load(f)
            already_processed = set(resume_data.get("processed", []))
            logger.info(f"Resuming: {len(already_processed)} proteins already processed")

            # Load existing residues
            existing_csv = config.processed_data_dir / "all_residues.csv"
            if existing_csv.exists():
                existing_df = pd.read_csv(existing_csv)
                existing_residues = existing_df.to_dict('records')
                logger.info(f"Loaded {len(existing_residues):,} existing residues")

    # Step 2: Process proteins
    session = requests.Session()
    session.proxies = PROXY

    lipid_extractor = LipidContactExtractor(config.distance_cutoff)

    all_residues = existing_residues.copy()
    sequences = {}

    # Track lipid types found
    lipid_type_counts: Dict[str, int] = {}

    proteins_remaining = [p for p in proteins_to_process if p not in already_processed]
    logger.info(f"Processing {len(proteins_remaining)} proteins...")

    for pdb_id in tqdm(proteins_remaining, desc="Processing proteins"):
        residues, sequence, error = process_protein(
            pdb_id, pdb_dir, session, lipid_extractor
        )

        if error:
            stats.processing_errors.append(error)
            stats.proteins_failed += 1
            if len(stats.processing_errors) <= 10:
                logger.debug(error)

        if residues:
            all_residues.extend(residues)
            stats.proteins_processed += 1

            if any(r["is_contact"] for r in residues):
                stats.proteins_with_contacts += 1

                # Track lipid types
                for lipid_code in pdb_to_lipids.get(pdb_id, []):
                    lipid_type_counts[lipid_code] = lipid_type_counts.get(lipid_code, 0) + 1

        if sequence:
            sequences[pdb_id] = sequence

        already_processed.add(pdb_id)

        # Checkpoint every 500 proteins
        if len(already_processed) % 500 == 0:
            with open(processed_file, "w") as f:
                json.dump({"processed": list(already_processed)}, f)
            logger.info(f"Checkpoint: {len(already_processed)} proteins processed")

    # Save final checkpoint
    with open(processed_file, "w") as f:
        json.dump({"processed": list(already_processed)}, f)

    stats.lipid_types_found = lipid_type_counts
    stats.proteins_processed = len(already_processed)

    logger.info(f"\nProcessed {stats.proteins_processed} proteins")
    logger.info(f"  With contacts: {stats.proteins_with_contacts}")
    logger.info(f"  Failed: {stats.proteins_failed}")
    logger.info(f"  Total residues: {len(all_residues):,}")

    if len(all_residues) == 0:
        logger.error("No residues extracted!")
        return

    # Step 3: Create DataFrame
    df = pd.DataFrame(all_residues)
    stats.total_residues = len(df)
    stats.positive_residues = int(df["is_contact"].sum())

    logger.info(f"Total residues: {stats.total_residues:,}")
    logger.info(f"Positive residues: {stats.positive_residues:,}")
    logger.info(f"Positive rate: {100*stats.positive_residues/stats.total_residues:.3f}%")

    # Step 4: Cluster sequences
    logger.info("\nClustering sequences by 30% identity...")
    clusterer = SequenceClusterer(config.sequence_identity_threshold, config)

    if sequences:
        clusters = clusterer.cluster_sequences(sequences)
        stats.n_clusters = len(set(clusters.values()))
        df["cluster_id"] = df["pdb_id"].map(clusters)
        logger.info(f"Created {stats.n_clusters} sequence clusters")
    else:
        logger.warning("No sequences available for clustering")

    # Step 5: Split into train/val/test
    logger.info("Creating train/val/test splits...")
    if "cluster_id" in df.columns and len(sequences) > 0:
        # Calculate per-protein contact rate for stratification
        protein_labels = df.groupby("pdb_id")["is_contact"].mean().apply(
            lambda x: 1 if x > 0 else 0
        ).to_dict()

        splits = stratified_cluster_split(
            clusters, protein_labels,
            config.train_ratio, config.val_ratio, config.test_ratio
        )
        df["split"] = df["pdb_id"].map(splits)

        # Calculate split statistics
        for split_name in ["train", "val", "test"]:
            split_df = df[df["split"] == split_name]
            n_proteins = split_df["pdb_id"].nunique()
            n_residues = len(split_df)
            n_positives = int(split_df["is_contact"].sum())

            if split_name == "train":
                stats.train_proteins = n_proteins
                stats.train_residues = n_residues
                stats.train_positives = n_positives
            elif split_name == "val":
                stats.val_proteins = n_proteins
                stats.val_residues = n_residues
                stats.val_positives = n_positives
            else:
                stats.test_proteins = n_proteins
                stats.test_residues = n_residues
                stats.test_positives = n_positives

            rate = 100 * n_positives / n_residues if n_residues > 0 else 0
            logger.info(f"  {split_name}: {n_proteins} proteins, {n_residues:,} residues, {n_positives:,} positives ({rate:.2f}%)")

    # Step 6: Save results
    output_dir = config.processed_data_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save full dataset
    df.to_csv(output_dir / "all_residues.csv", index=False)
    logger.info(f"Saved all residues to {output_dir / 'all_residues.csv'}")

    # Save split-specific files
    if "split" in df.columns:
        for split in ["train", "val", "test"]:
            split_df = df[df["split"] == split]
            split_df.to_csv(output_dir / f"{split}_residues.csv", index=False)

    # Save protein-level metadata
    protein_meta = df.groupby("pdb_id").agg({
        "is_contact": ["sum", "count"],
        "chain_id": "first",
        "cluster_id": "first",
        "split": "first"
    }).reset_index()
    protein_meta.columns = ["pdb_id", "contact_count", "residue_count", "chain_id", "cluster_id", "split"]
    protein_meta["contact_rate"] = protein_meta["contact_count"] / protein_meta["residue_count"]
    protein_meta.to_csv(output_dir / "protein_metadata.csv", index=False)

    # Save comprehensive statistics
    stats_dict = {
        "pipeline_version": "RCSB-only (no OPM dependency)",
        "label_source": "CRYSTALLOGRAPHIC",
        "distance_cutoff_angstrom": config.distance_cutoff,
        "residue_representation": "C-alpha atom",
        "timestamp": datetime.now().isoformat(),

        "total_rcsb_lipid_structures": stats.total_rcsb_lipid_structures,
        "proteins_processed": stats.proteins_processed,
        "proteins_with_contacts": stats.proteins_with_contacts,
        "proteins_failed": stats.proteins_failed,

        "total_residues": stats.total_residues,
        "positive_residues": stats.positive_residues,
        "positive_rate": stats.positive_residues / stats.total_residues if stats.total_residues > 0 else 0,

        "n_clusters": stats.n_clusters,
        "sequence_identity_threshold": config.sequence_identity_threshold,

        "train": {
            "proteins": stats.train_proteins,
            "residues": stats.train_residues,
            "positives": stats.train_positives,
            "positive_rate": stats.train_positives / stats.train_residues if stats.train_residues > 0 else 0
        },
        "val": {
            "proteins": stats.val_proteins,
            "residues": stats.val_residues,
            "positives": stats.val_positives,
            "positive_rate": stats.val_positives / stats.val_residues if stats.val_residues > 0 else 0
        },
        "test": {
            "proteins": stats.test_proteins,
            "residues": stats.test_residues,
            "positives": stats.test_positives,
            "positive_rate": stats.test_positives / stats.test_residues if stats.test_residues > 0 else 0
        },

        "lipid_types_found": stats.lipid_types_found,
        "total_lipid_codes_queried": len(LIPID_RESIDUE_NAMES),

        "comparison_to_dreamm": {
            "dreamm_proteins": 54,
            "dreamm_residues": 12805,
            "dreamm_label_source": "Literature curation (EPR, fluorescence, mutagenesis)",
            "mplid_proteins": stats.proteins_with_contacts,
            "mplid_residues": stats.total_residues,
            "mplid_label_source": "Crystallographic (4.0Å Cα-to-lipid)",
            "scale_increase_proteins": f"{stats.proteins_with_contacts / 54:.1f}x",
            "scale_increase_residues": f"{stats.total_residues / 12805:.1f}x"
        }
    }

    with open(output_dir / "pipeline_stats.json", "w") as f:
        json.dump(stats_dict, f, indent=2)

    # Print summary
    logger.info("\n" + "=" * 70)
    logger.info("RCSB-ONLY PIPELINE COMPLETE")
    logger.info("=" * 70)
    logger.info(f"Total RCSB structures queried: {stats.total_rcsb_lipid_structures:,}")
    logger.info(f"Proteins processed: {stats.proteins_processed:,}")
    logger.info(f"Proteins with contacts: {stats.proteins_with_contacts:,}")
    logger.info(f"Total residues: {stats.total_residues:,}")
    logger.info(f"Positive residues: {stats.positive_residues:,}")
    logger.info(f"Positive rate: {100*stats.positive_residues/stats.total_residues:.3f}%")
    logger.info(f"Sequence clusters: {stats.n_clusters}")
    logger.info(f"")
    logger.info(f"Train: {stats.train_proteins} proteins, {stats.train_residues:,} residues")
    logger.info(f"Val:   {stats.val_proteins} proteins, {stats.val_residues:,} residues")
    logger.info(f"Test:  {stats.test_proteins} proteins, {stats.test_residues:,} residues")
    logger.info(f"")
    logger.info(f"Comparison to DREAMM (54 proteins, 12,805 residues):")
    logger.info(f"  Scale increase: {stats.proteins_with_contacts/54:.1f}x proteins, {stats.total_residues/12805:.1f}x residues")
    logger.info(f"")
    logger.info(f"Output saved to: {output_dir}")
    logger.info(f"Log saved to: {log_file}")
    logger.info("=" * 70)


if __name__ == "__main__":
    main()

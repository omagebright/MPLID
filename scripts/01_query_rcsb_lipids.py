#!/usr/bin/env python3
"""Query RCSB PDB for membrane proteins with crystallized lipids.

This script queries the RCSB PDB database to identify membrane protein structures
that contain crystallized lipid molecules. It cross-references with the OPM
database to ensure only validated membrane proteins are included.

Author: Folorunsho Bright Omage
License: MIT
"""

import argparse
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Set

import requests

# Recognized lipid codes for RCSB query
LIPID_CODES = [
    # Phospholipids
    "CDL", "POV", "PCW", "PEE", "PGV", "PLC", "POP", "PPE", "PGP",
    "PIO", "PIP", "P2E", "P3E", "P4E", "P5E", "PC1", "PE5", "LPE", "LHG",
    # Sphingolipids
    "SPH", "S1P", "HXJ",
    # Sterols
    "CLR", "CHD", "Y01",
    # Fatty acids
    "MYR", "OLA", "STE", "ARA", "DHA", "PLM", "MYS", "PAM", "LNL",
    # Detergents (membrane mimetics)
    "LDA", "LMT", "BOG", "OLC", "DPC", "DMU",
]

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"


def setup_logging(output_dir: Path) -> logging.Logger:
    """Configure logging for the pipeline."""
    log_file = output_dir / "rcsb_query.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(),
        ],
    )
    return logging.getLogger(__name__)


def build_lipid_query(lipid_codes: List[str]) -> Dict:
    """Build RCSB search query for structures containing specified lipids.

    Parameters
    ----------
    lipid_codes : List[str]
        List of lipid residue codes to search for.

    Returns
    -------
    Dict
        RCSB search API query object.
    """
    lipid_groups = [
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_nonpolymer_instance_container_identifiers.comp_id",
                "operator": "exact_match",
                "value": code,
            },
        }
        for code in lipid_codes
    ]

    query = {
        "query": {
            "type": "group",
            "logical_operator": "or",
            "nodes": lipid_groups,
        },
        "return_type": "entry",
        "request_options": {
            "results_content_type": ["experimental"],
            "return_all_hits": True,
        },
    }
    return query


def query_rcsb(query: Dict, timeout: int = 120) -> List[str]:
    """Execute RCSB search query.

    Parameters
    ----------
    query : Dict
        RCSB search API query object.
    timeout : int, optional
        Request timeout in seconds (default: 120).

    Returns
    -------
    List[str]
        List of PDB identifiers matching the query.
    """
    response = requests.post(
        RCSB_SEARCH_URL,
        json=query,
        headers={"Content-Type": "application/json"},
        timeout=timeout,
    )
    response.raise_for_status()

    results = response.json()
    pdb_ids = [hit["identifier"] for hit in results.get("result_set", [])]
    return pdb_ids


def load_opm_proteins(opm_file: Path) -> Set[str]:
    """Load PDB identifiers from OPM database file.

    Parameters
    ----------
    opm_file : Path
        Path to OPM database file (JSON or text format).

    Returns
    -------
    Set[str]
        Set of PDB identifiers from OPM.
    """
    if opm_file.suffix == ".json":
        with open(opm_file) as f:
            data = json.load(f)
            return {entry["pdb_id"].upper() for entry in data}
    else:
        with open(opm_file) as f:
            return {line.strip().upper() for line in f if line.strip()}


def main():
    """Main entry point for RCSB lipid query pipeline."""
    parser = argparse.ArgumentParser(
        description="Query RCSB for membrane proteins with crystallized lipids"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/interim"),
        help="Output directory for query results",
    )
    parser.add_argument(
        "--opm-file",
        type=Path,
        help="Optional OPM database file for cross-referencing",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=120,
        help="RCSB query timeout in seconds",
    )
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    logger = setup_logging(args.output)

    logger.info("Starting RCSB lipid query")
    logger.info(f"Searching for {len(LIPID_CODES)} lipid codes")

    # Build and execute query
    query = build_lipid_query(LIPID_CODES)
    pdb_ids = query_rcsb(query, timeout=args.timeout)
    logger.info(f"Found {len(pdb_ids)} structures with lipids in RCSB")

    # Cross-reference with OPM if provided
    if args.opm_file and args.opm_file.exists():
        opm_ids = load_opm_proteins(args.opm_file)
        logger.info(f"Loaded {len(opm_ids)} proteins from OPM")

        intersection = set(pdb_ids) & opm_ids
        logger.info(f"Intersection (membrane proteins with lipids): {len(intersection)}")
        pdb_ids = sorted(intersection)

    # Save results
    output_file = args.output / "rcsb_lipid_structures.json"
    results = {
        "query_date": datetime.now().isoformat(),
        "lipid_codes_queried": LIPID_CODES,
        "total_structures": len(pdb_ids),
        "pdb_ids": pdb_ids,
    }

    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)

    logger.info(f"Results saved to {output_file}")
    logger.info("RCSB query complete")


if __name__ == "__main__":
    main()

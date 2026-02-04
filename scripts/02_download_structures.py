#!/usr/bin/env python3
"""Download PDB structures for lipid contact analysis.

This script downloads PDB structure files from the RCSB database for
proteins identified as containing crystallized lipids.

Author: Folorunsho Bright Omage
License: MIT
"""

import argparse
import gzip
import json
import logging
import time
from pathlib import Path
from typing import List, Optional

import requests
from tqdm import tqdm

RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb.gz"


def setup_logging(output_dir: Path) -> logging.Logger:
    """Configure logging for the pipeline."""
    log_file = output_dir / "download.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(),
        ],
    )
    return logging.getLogger(__name__)


def download_pdb(
    pdb_id: str,
    output_dir: Path,
    timeout: int = 30,
    retries: int = 3,
) -> Optional[Path]:
    """Download a single PDB structure file.

    Parameters
    ----------
    pdb_id : str
        PDB identifier (4 characters).
    output_dir : Path
        Directory to save downloaded files.
    timeout : int, optional
        Request timeout in seconds (default: 30).
    retries : int, optional
        Number of retry attempts (default: 3).

    Returns
    -------
    Optional[Path]
        Path to downloaded file, or None if download failed.
    """
    output_file = output_dir / f"{pdb_id.upper()}.pdb"

    if output_file.exists():
        return output_file

    url = RCSB_DOWNLOAD_URL.format(pdb_id=pdb_id.upper())

    for attempt in range(retries):
        try:
            response = requests.get(url, timeout=timeout)
            response.raise_for_status()

            # Decompress gzipped content
            content = gzip.decompress(response.content)
            output_file.write_bytes(content)
            return output_file

        except requests.RequestException as e:
            if attempt < retries - 1:
                time.sleep(2 ** attempt)  # Exponential backoff
            else:
                logging.warning(f"Failed to download {pdb_id}: {e}")
                return None


def load_pdb_list(input_file: Path) -> List[str]:
    """Load list of PDB identifiers from JSON file.

    Parameters
    ----------
    input_file : Path
        Path to JSON file containing PDB identifiers.

    Returns
    -------
    List[str]
        List of PDB identifiers.
    """
    with open(input_file) as f:
        data = json.load(f)
    return data.get("pdb_ids", [])


def main():
    """Main entry point for structure download pipeline."""
    parser = argparse.ArgumentParser(
        description="Download PDB structures for lipid contact analysis"
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="JSON file containing PDB identifiers",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/raw/pdb"),
        help="Output directory for PDB files",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=30,
        help="Download timeout in seconds per file",
    )
    parser.add_argument(
        "--max-structures",
        type=int,
        help="Maximum number of structures to download (for testing)",
    )
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    logger = setup_logging(args.output)

    # Load PDB list
    pdb_ids = load_pdb_list(args.input)
    if args.max_structures:
        pdb_ids = pdb_ids[:args.max_structures]

    logger.info(f"Downloading {len(pdb_ids)} PDB structures")

    # Download structures
    successful = 0
    failed = []

    for pdb_id in tqdm(pdb_ids, desc="Downloading"):
        result = download_pdb(pdb_id, args.output, timeout=args.timeout)
        if result:
            successful += 1
        else:
            failed.append(pdb_id)

    # Save summary
    summary = {
        "total_requested": len(pdb_ids),
        "successful": successful,
        "failed_count": len(failed),
        "failed_ids": failed,
    }

    summary_file = args.output / "download_summary.json"
    with open(summary_file, "w") as f:
        json.dump(summary, f, indent=2)

    logger.info(f"Download complete: {successful}/{len(pdb_ids)} successful")
    if failed:
        logger.warning(f"Failed downloads: {failed}")


if __name__ == "__main__":
    main()

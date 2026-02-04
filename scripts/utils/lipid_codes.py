"""Lipid code definitions for MPLID dataset.

This module contains the complete list of recognized lipid codes
from the PDB Chemical Component Dictionary.

Author: Folorunsho Bright Omage
License: MIT
"""

from typing import Dict, Set

# Primary lipid codes used for RCSB queries
LIPID_QUERY_CODES = [
    "CDL", "POV", "PCW", "PEE", "PGV", "PLM", "SPH", "S1P",
    "CLR", "CHD", "Y01", "BCL", "MYR", "OLA", "STE", "ARA",
    "DHA", "LDA", "LMT", "BOG", "OLC", "DPC", "DMU", "EPE",
    "LPE", "PLC", "POP", "PPE", "PGP", "UNL", "HXJ",
]

# Complete lipid code dictionary with categories
LIPID_CATEGORIES: Dict[str, Set[str]] = {
    "phospholipids": {
        "CDL",  # Cardiolipin
        "POV",  # 1-palmitoyl-2-oleoyl-phosphatidylcholine
        "PCW",  # Phosphatidylcholine
        "PEE",  # Phosphatidylethanolamine
        "PGV",  # Phosphatidylglycerol
        "PLM",  # Palmitic acid
        "PLC",  # Dipalmitoyl phosphatidylcholine
        "POP",  # 1-palmitoyl-2-oleoyl-phosphatidylethanolamine
        "PPE",  # Palmitoyl-2-oleoyl-sn-glycero-3-phosphoethanolamine
        "PGP",  # Phosphatidylglycerol phosphate
        "PIO",  # Phosphatidylinositol
        "PIP",  # Phosphatidylinositol-4-phosphate
        "P2E",  # Phosphatidylinositol-4,5-bisphosphate
        "P3E",  # Phosphatidylinositol-3,4-bisphosphate
        "P4E",  # Phosphatidylinositol-3,4,5-trisphosphate
        "P5E",  # Phosphatidylinositol-5-phosphate
        "PC1",  # Phosphatidylcholine variant
        "PE5",  # Phosphatidylethanolamine variant
        "LPE",  # Lysophosphatidylethanolamine
        "EPE",  # 1,2-dioleoyl-sn-glycero-3-phosphoethanolamine
        "LPC",  # Lysophosphatidylcholine
        "LPG",  # Lysophosphatidylglycerol
        "PGE",  # Phosphatidylglycerol
        "PSE",  # Phosphatidylserine
        "LPS",  # Lysophosphatidylserine
    },
    "sphingolipids": {
        "SPH",  # Sphingosine
        "S1P",  # Sphingosine-1-phosphate
        "HXJ",  # Hexosylceramide
        "CER",  # Ceramide
        "SM",   # Sphingomyelin
        "GLC",  # Glucosylceramide
        "GAL",  # Galactosylceramide
    },
    "sterols": {
        "CLR",  # Cholesterol
        "CHD",  # Cholesteryl hemisuccinate
        "Y01",  # Cholesterol hemisuccinate
        "BCL",  # Beta-sitosterol
        "CHL",  # Chlorophyll
        "LHG",  # Ergosterol
        "ERG",  # Ergosterol
        "STR",  # Stigmasterol
        "SIT",  # Sitosterol
        "CPS",  # Campesterol
    },
    "fatty_acids": {
        "MYR",  # Myristic acid (C14:0)
        "OLA",  # Oleic acid (C18:1)
        "STE",  # Stearic acid (C18:0)
        "ARA",  # Arachidonic acid (C20:4)
        "DHA",  # Docosahexaenoic acid (C22:6)
        "MYS",  # Myristoyl
        "PAM",  # Palmitate (C16:0)
        "LNL",  # Linoleic acid (C18:2)
        "LNA",  # Alpha-linolenic acid (C18:3)
        "EPA",  # Eicosapentaenoic acid (C20:5)
        "LAU",  # Lauric acid (C12:0)
        "CAP",  # Capric acid (C10:0)
    },
    "detergents": {
        "LDA",  # Lauryl dimethylamine-N-oxide
        "LMT",  # n-Dodecyl-beta-D-maltoside
        "BOG",  # Beta-octyl glucoside
        "OLC",  # n-Octyl-beta-D-glucopyranoside
        "DPC",  # n-Dodecylphosphocholine
        "DMU",  # n-Decyl-beta-D-maltoside
        "UNL",  # Unknown ligand (often lipid-like)
        "C8E",  # Octyl-polyoxyethylene
        "C10E", # Decyl-polyoxyethylene
        "C12E", # Dodecyl-polyoxyethylene
        "SDS",  # Sodium dodecyl sulfate
        "TX1",  # Triton X-100
        "DDM",  # n-Dodecyl-beta-D-maltoside
        "OG",   # n-Octyl-beta-D-glucoside
        "NG",   # n-Nonyl-beta-D-glucoside
        "LDAO", # Lauryldimethylamine oxide
        "CHAPS",# CHAPS
        "CHAPSO",# CHAPSO
    },
    "other": {
        "TGL",  # Triglyceride
        "DGD",  # Digalactosyldiacylglycerol
        "MGD",  # Monogalactosyldiacylglycerol
        "DAG",  # Diacylglycerol
        "LPL",  # Lipoprotein lipase
        "BCD",  # Beta-cyclodextrin (lipid mimetic)
    },
}

# All recognized lipid codes as a flat set
LIPID_CODES: Set[str] = set()
for category_codes in LIPID_CATEGORIES.values():
    LIPID_CODES.update(category_codes)


def get_lipid_category(code: str) -> str:
    """Get the category for a lipid code.

    Parameters
    ----------
    code : str
        Three-letter lipid code.

    Returns
    -------
    str
        Category name, or 'unknown' if not found.
    """
    for category, codes in LIPID_CATEGORIES.items():
        if code.upper() in codes:
            return category
    return "unknown"


def is_lipid(code: str) -> bool:
    """Check if a residue code is a recognized lipid.

    Parameters
    ----------
    code : str
        Three-letter residue code.

    Returns
    -------
    bool
        True if the code is a recognized lipid.
    """
    return code.upper() in LIPID_CODES

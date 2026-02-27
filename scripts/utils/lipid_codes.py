"""Lipid code definitions for MPLID dataset.

This module contains the complete list of 117 recognized lipid codes
used in the MPLID contact calculation pipeline. Codes are sourced from
the PDB Chemical Component Dictionary and CHARMM lipid simulation
nomenclature commonly found in cryo-EM depositions.

Author: Folorunsho Bright Omage
License: MIT
"""

from typing import Dict, Set

# Primary lipid codes used for RCSB structure queries.
# This is a focused subset for efficient RCSB Search API queries;
# the full recognition set (LIPID_CODES, 117 codes) is used during
# contact calculation.
LIPID_QUERY_CODES = [
    # Phospholipids
    "CDL", "PC1", "PCW", "PEE", "PGV", "POV", "PLM", "PLC", "LHG", "PTY",
    # Sterols
    "CLR", "CHD", "Y01",
    # Sphingolipids
    "SPH", "S1P", "HXJ",
    # Fatty acids
    "MYR", "OLA", "STE", "ARA", "DHA",
    # Detergents
    "LDA", "LMT", "BOG", "OLC", "DPC",
    # Additional common codes
    "C9V",
]

# Complete lipid code dictionary with categories.
# These 117 codes match the production pipeline exactly.
LIPID_CATEGORIES: Dict[str, Set[str]] = {
    "phospholipids": {
        # Phosphatidylcholine (PC) variants
        "PC1",   # 1,2-diacyl-sn-glycero-3-phosphocholine
        "PCW",   # 1,2-dioleoyl-sn-glycero-3-phosphocholine (DOPC)
        "POV",   # Palmitoyl-oleoyl-phosphatidylcholine (POPC)
        "PLC",   # Diundecyl phosphatidylcholine
        "1PL",   # 1-palmitoyl-sn-glycero-3-phosphocholine
        "2PL",   # 1,2-dipalmitoyl-sn-glycero-3-phosphocholine
        "XPC",   # Phosphatidylcholine variant
        # Phosphatidylethanolamine (PE) variants
        "PEE",   # 1,2-dioleoyl-sn-glycero-3-phosphoethanolamine (DOPE)
        "PTY",   # Phosphatidylethanolamine
        "MPE",   # Myristic-palmitoyl PE
        "DPE",   # Dipalmitoyl PE
        # Phosphatidylglycerol (PG) variants
        "PGV",   # Phosphatidylglycerol
        "EPG",   # Egg PG
        "LHG",   # 1,2-dipalmitoyl-phosphatidylglycerol
        "PGR",   # PG variant
        "DPG",   # Diphosphatidylglycerol
        # Phosphatidylserine (PS) variants
        "LPS",   # Lysophosphatidylserine
        "DPS",   # Dipalmitoyl PS
        # Phosphatidylinositol (PI) variants
        "PTI",   # Phosphatidylinositol
        "IPC",   # Inositol phosphoceramide
        "PIP",   # PI phosphate
        "PI3",   # PI 3-phosphate
        "PI4",   # PI 4-phosphate
        "P34",   # PI 3,4-bisphosphate
    },
    "cardiolipin": {
        "CDL",   # Cardiolipin (canonical PDB code)
        "C9V",   # Cardiolipin-boron complex
        "18W",   # Tetraoleoyl cardiolipin
        "LCL",   # Lyso-cardiolipin
    },
    "sphingolipids": {
        "HXJ",   # Sphingomyelin/ceramide variant
        "SPH",   # Sphingosine
        "S1P",   # Sphingosine-1-phosphate
        "CRT",   # Ceramide
        "GCR",   # Glucosylceramide
        "GSP",   # Ganglioside
        "GLF",   # Galactosylceramide
    },
    "sterols": {
        "CLR",   # Cholesterol
        "CHD",   # Cholesterol derivative
        "Y01",   # Cholesterol hemisuccinate
        "OHC",   # 25-hydroxycholesterol
        "OCL",   # Cholesteryl oleate
        "LNS",   # Lanosterol
        "ERG",   # Ergosterol
        "SIT",   # Beta-sitosterol
        "STI",   # Stigmasterol
    },
    "fatty_acids": {
        "PLM",   # Palmitic acid (C16:0)
        "MYR",   # Myristic acid (C14:0)
        "OLA",   # Oleic acid (C18:1)
        "STE",   # Stearic acid (C18:0)
        "ARA",   # Arachidonic acid (C20:4)
        "DHA",   # Docosahexaenoic acid (C22:6)
        "LNL",   # Linoleic acid (C18:2)
        "LNN",   # Alpha-linolenic acid (C18:3)
        "EPA",   # Eicosapentaenoic acid (C20:5)
        "PAM",   # Palmitate ion
        "DCR",   # Decanoic acid
        "UND",   # Undecanoic acid
        "LAU",   # Lauric acid (C12:0)
        "CAP",   # Capric acid (C10:0)
    },
    "glycerolipids": {
        "TGL",   # Triglyceride
        "TAG",   # Triacylglycerol
        "DGA",   # Diglyceride
        "DAG",   # Diacylglycerol
        "MAG",   # Monoacylglycerol
        "GMO",   # Glyceryl monooleate
    },
    "detergents": {
        "LDA",   # Lauryl dimethylamine-N-oxide (LDAO)
        "LMT",   # Dodecyl-beta-D-maltoside (DDM)
        "OLC",   # Monoolein (glyceryl monooleate)
        "BOG",   # n-Octyl-beta-D-glucopyranoside
        "BGC",   # Beta-D-glucopyranosyl
        "C8E",   # Octyl-PEG
        "C10",   # Decyl maltoside
        "C12",   # Dodecyl maltoside
        "UDM",   # n-Undecyl-beta-D-maltoside
        "NM",    # n-Nonyl-beta-D-maltoside
        "DM",    # n-Decyl-beta-D-maltoside
        "HTG",   # Heptyl thioglucoside
        "OG",    # Octyl glucoside
        "NG",    # Nonyl glucoside
        "SDS",   # Sodium dodecyl sulfate
        "DPC",   # Dodecylphosphocholine
        "FC6",   # Fos-choline 6
        "FC8",   # Fos-choline 8
        "FC10",  # Fos-choline 10
        "FC12",  # Fos-choline 12
        "FC14",  # Fos-choline 14
        "FC16",  # Fos-choline 16
    },
    "lipid_a": {
        "LPA",   # Lipid A
        "KDO",   # 3-deoxy-D-manno-octulosonic acid
        "6LP",   # Lipid A disaccharide
    },
    "charmm_simulation": {
        # CHARMM/simulation lipid names (commonly found in cryo-EM depositions)
        "POPC", "POPE", "POPS", "POPG", "POPA",
        "DPPC", "DPPE", "DPPS", "DPPG",
        "DOPC", "DOPE", "DOPS", "DOPG",
        "DMPC", "DMPE", "DMPG",
        "DLPC", "DLPE", "DLPG",
        "SOPC", "SDPC", "SLPC",
        "PLPC", "PAPC", "SAPC",
        "CHL1", "CHOL",
    },
    "generic": {
        "LIP",   # Generic lipid
    },
}

# All recognized lipid codes as a flat set (117 total)
LIPID_CODES: Set[str] = set()
for _category_codes in LIPID_CATEGORIES.values():
    LIPID_CODES.update(_category_codes)


def get_lipid_category(code: str) -> str:
    """Get the category for a lipid code.

    Parameters
    ----------
    code : str
        Lipid residue code from PDB.

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
        Residue code from PDB.

    Returns
    -------
    bool
        True if the code is a recognized lipid.
    """
    return code.upper() in LIPID_CODES

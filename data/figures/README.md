# Figure Source Data

Underlying numerical data for the figures in the MPLID manuscript (*GigaScience*,
GIGA-D-26-00059), provided so each figure can be reproduced from tabular values.
**Figure numbers below refer to the figure numbering in the published article.**
All values derive from the deposited MPLID v2.1.0 dataset (4,704 proteins;
8,055,325 residues) released under a CC0 public domain dedication.

| File | Article figure | Contents |
|------|----------------|----------|
| `figure2_cluster_size_distribution.csv` | Fig. 2 (Sequence cluster analysis) | Number of proteins per sequence cluster (`cluster_id`, `n_proteins`) at 30% identity. Panel A (size histogram) and Panel B (cumulative coverage) are both derived from these 813 cluster sizes. |
| `figure4_distance_distribution.csv` | Fig. 4 (Distance distribution) | Histogram counts of the minimum residue-to-lipid heavy-atom distance in 0.25 Å bins, given separately for contact and non-contact residues (`bin_lower_angstrom`, `bin_upper_angstrom`, `contact_residue_count`, `noncontact_residue_count`). |
| `figure5_protein_size_contact_density.csv` | Fig. 5 (Protein-level statistics) | Per-protein size (`n_residues`) and contact density (`contact_rate = n_contacts / n_residues`) for all proteins — the values plotted in the size and contact-density distributions. |
| `figure8_opm_comparison_summary.csv` | Fig. 8A | Residue-level MPLID-vs-OPM overlap: contact / membrane-zone counts and the percentages shown in Panel A. |
| `figure8_opm_comparison_per_protein.csv` | Fig. 8B–D | Per-protein counts (`n_residues`, `n_contacts`, `n_membrane`, `n_both`, `n_mplid_only`, `n_opm_only`) and derived `opm_precision_pct` and `contact_recall_by_opm_pct` plotted in Panels B, C, and D. |

## Notes

- **Fig. 4** — Contact residues are defined as minimum distance ≤ 4.0 Å, so their
  counts fall entirely in bins up to 4.0 Å; non-contact residues populate bins above
  4.0 Å. Both histograms use the full dataset (no subsampling).
- **Fig. 5** — `contact_rate` is a fraction (multiply by 100 for a percentage).
- **Fig. 8** — Panel B is filtered to proteins with ≥5 OPM membrane-zone residues;
  Panel C to proteins with ≥3 lipid contacts; Panel D to proteins with ≥1 contact,
  matching the figure-generation script.
- Note: the figure image filenames in `manuscript/figures/` (e.g. `fig6_cluster_analysis`)
  reflect an internal working order and do **not** all match the final article figure
  numbers; the table above uses the article numbering.

Regenerated from the deposited dataset with the pipeline in `scripts/`; see the
repository README and `docs/METHODOLOGY.md` for full details.

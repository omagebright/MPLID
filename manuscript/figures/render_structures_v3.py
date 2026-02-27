#!/usr/bin/env python3
"""
PyMOL rendering script for MPLID manuscript Figure 7 (v3).
Publication-quality structural visualizations with high-contrast
lipid contacts and bold annotations.

Key changes from v2:
- Surface transparency reduced (0.82 -> 0.70) for more visible protein body
- Contact residues rendered WITHOUT surface overlay (fully opaque sticks)
- Thicker lipid sticks (0.22 -> 0.35) with brighter colors
- Contact stick radius increased (0.18 -> 0.25)
- Contact cartoon uses a more saturated red
- Larger, bolder text labels in composite figure
- Higher resolution (3000x2500)

Usage:
    python3 render_structures_v3.py
"""

import sys
import os

try:
    import pymol
    from pymol import cmd, stored
except ImportError:
    print("PyMOL not available.")
    sys.exit(1)

pymol.finish_launching(['pymol', '-cq'])

# Paths
STRUCT_DIR = '/HOMENFS/c60238528073/Art/Exo_Database/Lipid_membrane_protein_interfaces/MPLID_paper/manuscript/figures/structures'
OUTPUT_DIR = '/HOMENFS/c60238528073/Art/Exo_Database/Lipid_membrane_protein_interfaces/MPLID_paper/manuscript/figures'

# Molecules to remove (non-lipid HETATM)
REMOVE_HETATM = {'HOH', 'SO4', 'CXT', 'NA', 'CL', 'CA', 'ZN', 'MG',
                  '1PE', 'PEG', 'GOL', 'EDO', 'ACT', 'HTO', 'ZMA'}

# Colors - higher contrast than v2
COLORS = {
    'protein_cartoon': [0.78, 0.82, 0.88],     # Steel blue-gray (slightly darker)
    'contact_cartoon': [0.90, 0.22, 0.15],      # Strong red (more saturated)
    'contact_carbon':  [0.85, 0.18, 0.12],      # Deep red
    'lipid_carbon':    [0.15, 0.72, 0.30],       # Vivid green (high contrast vs red)
    'lipid_oxygen':    [0.85, 0.20, 0.20],       # Red
    'lipid_nitrogen':  [0.20, 0.35, 0.85],       # Blue
    'lipid_phosphorus': [1.00, 0.55, 0.00],      # Orange
    'surface_color':   [0.85, 0.88, 0.93],       # Light blue-gray
}

PROTEINS = {
    '5NM2': {
        'title': r'A$_{2A}$ Adenosine Receptor',
        'subtitle': 'GPCR',
        'contacts': {
            'A': [3, 5, 6, 7, 9, 11, 19, 22, 23, 28, 43, 47, 57, 61, 65, 68,
                  71, 72, 73, 76, 97, 101, 118, 122, 126, 130, 134, 141, 244,
                  247, 248, 251, 252, 254, 255, 258, 259, 262, 263, 265, 268]
        },
        'lipids': ['OLA', 'OLB', 'OLC', 'CLR'],
        'n_contacts': 41,
        'n_residues': 388,
        'restrict_chain': 'A',
        'restrict_resi_range': (1, 310),
    },
    '4C9J': {
        'title': 'ADP/ATP Carrier (AAC3)',
        'subtitle': 'Transporter',
        'contacts': {
            'A': [58, 59, 60, 77, 78, 80, 161, 176, 180, 181, 182, 183, 184,
                  256, 258, 276, 278],
            'B': [58, 59, 77, 78, 80, 161, 162, 176, 181, 182, 183, 256, 258,
                  275, 276, 278]
        },
        'lipids': ['CDL'],
        'n_contacts': 33,
        'n_residues': 544,
        'restrict_chain': None,
    },
    '5T77': {
        'title': 'MurJ Flippase',
        'subtitle': 'Lipid II Flippase',
        'contacts': {
            'A': [11, 20, 23, 85, 91, 94, 95, 100, 101, 153, 157, 164, 168,
                  179, 183, 195, 206, 207, 208, 209, 216, 224, 268, 269, 303,
                  307, 332, 336, 368, 404, 408, 437, 442]
        },
        'lipids': ['OLC', 'OLB'],
        'n_contacts': 33,
        'n_residues': 467,
        'restrict_chain': None,
    },
    '4TSY': {
        'title': 'Fragaceatoxin C (FraC)',
        'subtitle': 'Pore-forming Toxin',
        'contacts': {
            'A': [84, 85, 107, 108, 113, 167, 168],
            'B': [78, 84, 85, 107, 108, 167, 168],
            'C': [78, 79, 84, 85, 107, 167, 168],
            'D': [78, 79, 84, 85, 107, 167, 168]
        },
        'lipids': ['HXJ'],
        'n_contacts': 28,
        'n_residues': 704,
        'restrict_chain': None,
    }
}


def set_color(name, rgb):
    cmd.set_color(name, rgb)


def setup_rendering():
    """Configure PyMOL for high-contrast publication output."""
    cmd.bg_color('white')
    cmd.set('ray_opaque_background', 1)

    # Ray tracing - sharper outlines
    cmd.set('ray_trace_mode', 1)
    cmd.set('ray_trace_gain', 0.012)
    cmd.set('ray_trace_color', 'black')
    cmd.set('ray_shadows', 'off')
    cmd.set('antialias', 2)

    # Orthographic projection
    cmd.set('orthoscopic', 1)
    cmd.set('depth_cue', 0)
    cmd.set('fog', 0)

    # Lighting - slightly more directional for depth
    cmd.set('ambient', 0.35)
    cmd.set('direct', 0.65)
    cmd.set('reflect', 0.4)
    cmd.set('spec_reflect', 0.25)
    cmd.set('spec_power', 200)
    cmd.set('spec_direct', 0)

    # Cartoon quality
    cmd.set('cartoon_fancy_helices', 1)
    cmd.set('cartoon_highlight_color', -1)
    cmd.set('cartoon_smooth_loops', 1)
    cmd.set('cartoon_oval_length', 1.2)
    cmd.set('cartoon_oval_width', 0.25)
    cmd.set('cartoon_loop_radius', 0.22)
    cmd.set('cartoon_tube_radius', 0.35)
    cmd.set('cartoon_discrete_colors', 1)

    # Stick quality
    cmd.set('stick_ball', 0)
    cmd.set('stick_quality', 20)

    # Define custom colors
    for name, rgb in COLORS.items():
        set_color(name, rgb)


def build_contact_sel(contacts_dict):
    parts = []
    for chain, residues in contacts_dict.items():
        resi_str = '+'.join(str(r) for r in residues)
        parts.append(f'(chain {chain} and resi {resi_str})')
    return ' or '.join(parts)


def render_protein(pdb_id, info):
    """Render a single protein with high-contrast lipid contacts."""
    cmd.reinitialize()
    setup_rendering()

    pdb_file = os.path.join(STRUCT_DIR, f'{pdb_id}.pdb')
    if not os.path.exists(pdb_file):
        print(f"ERROR: {pdb_file} not found!")
        return None

    cmd.load(pdb_file, pdb_id)
    cmd.hide('everything')

    # Remove water and non-lipid HETATM
    for het in REMOVE_HETATM:
        cmd.remove(f'resn {het}')
    cmd.remove('hydrogens')

    # Restrict to region if needed
    if info.get('restrict_resi_range'):
        rmin, rmax = info['restrict_resi_range']
        chain = info.get('restrict_chain', 'A')
        cmd.remove(f'{pdb_id} and polymer and chain {chain} and (resi 0-{rmin-1} or resi {rmax+1}-9999)')
        if chain:
            cmd.remove(f'{pdb_id} and polymer and not chain {chain}')

    # Build selections
    lipid_sel_parts = [f'resn {ln}' for ln in info['lipids']]
    lipid_sel = f'{pdb_id} and ({" or ".join(lipid_sel_parts)})'
    contact_sel_str = build_contact_sel(info['contacts'])
    cmd.select('contacts_sel', f'{pdb_id} and polymer and ({contact_sel_str})')
    cmd.select('noncontact_sel', f'{pdb_id} and polymer and not contacts_sel')
    cmd.select('lipids_sel', lipid_sel)

    # ---- NON-CONTACT PROTEIN: cartoon + transparent surface ----
    cmd.show('cartoon', 'noncontact_sel')
    cmd.color('protein_cartoon', 'noncontact_sel')

    # Surface ONLY on non-contact residues (so contacts show through)
    cmd.show('surface', 'noncontact_sel')
    cmd.color('surface_color', 'noncontact_sel')
    cmd.set('transparency', 0.70, 'noncontact_sel')

    # ---- CONTACT RESIDUES: bold cartoon + thick opaque sticks (NO surface) ----
    cmd.show('cartoon', 'contacts_sel')
    cmd.color('contact_cartoon', 'contacts_sel')

    cmd.show('sticks', 'contacts_sel')
    cmd.set('stick_radius', 0.25, 'contacts_sel')
    cmd.color('contact_carbon', 'contacts_sel and elem C')
    cmd.color('lipid_oxygen', 'contacts_sel and elem O')
    cmd.color('lipid_nitrogen', 'contacts_sel and elem N')
    cmd.color('tv_yellow', 'contacts_sel and elem S')

    # ---- LIPID MOLECULES: thick, vivid sticks ----
    cmd.show('sticks', 'lipids_sel')
    cmd.set('stick_radius', 0.35, 'lipids_sel')
    cmd.color('lipid_carbon', 'lipids_sel and elem C')
    cmd.color('lipid_oxygen', 'lipids_sel and elem O')
    cmd.color('lipid_nitrogen', 'lipids_sel and elem N')
    cmd.color('lipid_phosphorus', 'lipids_sel and elem P')
    cmd.color('tv_yellow', 'lipids_sel and elem S')

    # ---- Orient ----
    cmd.orient(f'{pdb_id} and polymer')
    cmd.zoom(f'{pdb_id}', buffer=4)

    # ---- Ray trace at high resolution ----
    width, height = 3000, 2500
    cmd.ray(width, height)

    out_png = os.path.join(OUTPUT_DIR, f'fig7_panel_{pdb_id}.png')
    cmd.png(out_png, dpi=300)
    print(f"  Saved: {out_png}")

    return out_png


def create_composite_figure(panel_data):
    """Create publication-quality 2x2 composite with legible labels."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    import matplotlib.patches as mpatches

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.subplots_adjust(wspace=0.04, hspace=0.14, left=0.01, right=0.99,
                        top=0.93, bottom=0.08)

    labels = ['A)', 'B)', 'C)', 'D)']

    for idx, (pdb_id, png_path, info) in enumerate(panel_data):
        ax = axes[idx // 2, idx % 2]

        if os.path.exists(png_path):
            img = mpimg.imread(png_path)
            ax.imshow(img)

        ax.axis('off')

        # Bold panel label - large
        ax.text(0.02, 0.97, labels[idx],
                transform=ax.transAxes, fontsize=26, fontweight='bold',
                va='top', ha='left', family='sans-serif',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                          edgecolor='none', alpha=0.9))

        # Protein info - larger, bolder text
        lipid_names = ', '.join(info['lipids'])
        title_str = f"PDB: {pdb_id}  |  {info['title']}"
        detail_str = f"{info['n_contacts']} contact residues  |  Lipids: {lipid_names}"

        ax.text(0.50, 0.01, detail_str,
                transform=ax.transAxes, fontsize=12,
                va='bottom', ha='center', family='sans-serif',
                color='#333333',
                bbox=dict(boxstyle='round,pad=0.35', facecolor='white',
                          edgecolor='#bbbbbb', alpha=0.95, linewidth=0.8))

        ax.text(0.50, 0.065, title_str,
                transform=ax.transAxes, fontsize=14, fontweight='bold',
                va='bottom', ha='center', family='sans-serif',
                color='#111111',
                bbox=dict(boxstyle='round,pad=0.35', facecolor='white',
                          edgecolor='#999999', alpha=0.95, linewidth=0.8))

    # Color legend at bottom
    legend_items = [
        mpatches.Patch(facecolor=[0.78, 0.82, 0.88], edgecolor='gray',
                       label='Non-contact residues'),
        mpatches.Patch(facecolor=[0.90, 0.22, 0.15], edgecolor='gray',
                       label='Lipid contact residues'),
        mpatches.Patch(facecolor=[0.15, 0.72, 0.30], edgecolor='gray',
                       label='Lipid molecules'),
    ]
    fig.legend(handles=legend_items, loc='lower center', ncol=3,
               fontsize=13, frameon=True, fancybox=True,
               edgecolor='#cccccc', facecolor='white',
               bbox_to_anchor=(0.5, 0.01))

    # Save
    for ext in ['png', 'pdf']:
        out_path = os.path.join(OUTPUT_DIR, f'fig7_structural_visualization.{ext}')
        fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"  Composite: {out_path}")

    plt.close()


def main():
    print("MPLID Figure 7 - Structural Visualization (v3)")
    print("=" * 60)

    panel_data = []
    for pdb_id, info in PROTEINS.items():
        print(f"\nRendering {pdb_id}: {info['title']} ({info['subtitle']})")
        print(f"  {info['n_contacts']} contacts, lipids: {', '.join(info['lipids'])}")
        out = render_protein(pdb_id, info)
        if out:
            panel_data.append((pdb_id, out, info))

    if panel_data:
        print("\nCreating composite figure...")
        create_composite_figure(panel_data)

    print("\nDone!")


if __name__ == '__main__':
    main()

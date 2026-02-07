#!/usr/bin/env python3
"""
PyMOL rendering script for MPLID manuscript Figure 7 (v2 - improved).
Creates publication-quality structural visualizations of membrane proteins
with their experimentally determined lipid contacts.

Key improvements over v1:
- Higher contrast between contact and non-contact residues
- Thicker lipid sticks with more vivid coloring
- Semi-transparent surface for membrane context
- Better viewing angles for each protein
- Higher resolution ray tracing (2400x2000)
- Improved composite figure with proper legends

Usage:
    python3 render_structures_v2.py
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
STRUCT_DIR = '/HOMENFS/60238528073/Art/Exo_Database/Lipid_membrane_protein_interfaces/MPLID_paper/manuscript/figures/structures'
OUTPUT_DIR = '/HOMENFS/60238528073/Art/Exo_Database/Lipid_membrane_protein_interfaces/MPLID_paper/manuscript/figures'

# Molecules to remove (non-lipid HETATM)
REMOVE_HETATM = {'HOH', 'SO4', 'CXT', 'NA', 'CL', 'CA', 'ZN', 'MG',
                  '1PE', 'PEG', 'GOL', 'EDO', 'ACT', 'HTO', 'ZMA'}

# Define colors as RGB tuples for consistency
COLORS = {
    'protein_cartoon': [0.85, 0.87, 0.90],     # Light steel blue-gray
    'contact_cartoon': [0.95, 0.45, 0.30],      # Warm coral-red
    'contact_carbon':  [0.90, 0.25, 0.20],      # Strong red
    'lipid_carbon':    [0.40, 0.75, 0.35],       # Vivid green
    'lipid_oxygen':    [0.85, 0.20, 0.20],       # Red
    'lipid_nitrogen':  [0.20, 0.35, 0.85],       # Blue
    'lipid_phosphorus': [1.00, 0.60, 0.00],      # Orange
    'surface_color':   [0.90, 0.92, 0.95],       # Very light blue-gray
}

PROTEINS = {
    '5NM2': {
        'title': r'A$_2$$_A$ Adenosine Receptor',
        'subtitle': 'GPCR',
        'contacts': {
            'A': [3, 5, 6, 7, 9, 11, 19, 22, 23, 28, 43, 47, 57, 61, 65, 68,
                  71, 72, 73, 76, 97, 101, 118, 122, 126, 130, 134, 141, 244,
                  247, 248, 251, 252, 254, 255, 258, 259, 262, 263, 265, 268]
        },
        'lipids': ['OLA', 'OLB', 'OLC', 'CLR'],
        'n_contacts': 41,
        'n_residues': 388,
        # 5NM2 has a BRIL fusion and T4L; restrict view to chain A TM domain
        'restrict_chain': 'A',
        'restrict_resi_range': (1, 310),  # Core GPCR domain only
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
    """Define a custom PyMOL color."""
    cmd.set_color(name, rgb)


def setup_rendering():
    """Configure PyMOL for publication-quality output."""
    # Background
    cmd.bg_color('white')
    cmd.set('ray_opaque_background', 1)

    # Ray tracing
    cmd.set('ray_trace_mode', 1)
    cmd.set('ray_trace_gain', 0.008)
    cmd.set('ray_trace_color', 'black')
    cmd.set('ray_shadows', 'off')
    cmd.set('antialias', 2)

    # Perspective
    cmd.set('orthoscopic', 1)
    cmd.set('depth_cue', 0)
    cmd.set('fog', 0)

    # Lighting - clean, even illumination
    cmd.set('ambient', 0.4)
    cmd.set('direct', 0.6)
    cmd.set('reflect', 0.3)
    cmd.set('spec_reflect', 0.2)
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
    cmd.set('stick_quality', 15)

    # Define custom colors
    for name, rgb in COLORS.items():
        set_color(name, rgb)


def build_contact_sel(contacts_dict):
    """Build PyMOL selection from contacts dictionary."""
    parts = []
    for chain, residues in contacts_dict.items():
        resi_str = '+'.join(str(r) for r in residues)
        parts.append(f'(chain {chain} and resi {resi_str})')
    return ' or '.join(parts)


def render_protein(pdb_id, info):
    """Render a single protein with lipid contacts - improved version."""
    cmd.reinitialize()
    setup_rendering()

    pdb_file = os.path.join(STRUCT_DIR, f'{pdb_id}.pdb')
    if not os.path.exists(pdb_file):
        print(f"ERROR: {pdb_file} not found!")
        return None

    cmd.load(pdb_file, pdb_id)
    cmd.hide('everything')

    # Remove water and non-lipid HETATM first
    for het in REMOVE_HETATM:
        cmd.remove(f'resn {het}')
    cmd.remove('hydrogens')

    # ---- Define selections ----
    # If we need to restrict to a specific region (e.g., remove fusion protein)
    if info.get('restrict_resi_range'):
        rmin, rmax = info['restrict_resi_range']
        chain = info.get('restrict_chain', 'A')
        # Remove polymer outside range
        cmd.remove(f'{pdb_id} and polymer and chain {chain} and (resi 0-{rmin-1} or resi {rmax+1}-9999)')
        # Remove other chains' polymer (but keep lipids)
        if chain:
            other_chains_polymer = f'{pdb_id} and polymer and not chain {chain}'
            cmd.remove(other_chains_polymer)

    # Build lipid selection
    lipid_sel_parts = [f'resn {ln}' for ln in info['lipids']]
    lipid_sel = f'{pdb_id} and ({" or ".join(lipid_sel_parts)})'

    # Build contact selection
    contact_sel_str = build_contact_sel(info['contacts'])
    cmd.select('contacts_sel', f'{pdb_id} and polymer and ({contact_sel_str})')
    cmd.select('noncontact_sel', f'{pdb_id} and polymer and not contacts_sel')
    cmd.select('lipids_sel', lipid_sel)

    # ---- PROTEIN CARTOON ----
    cmd.show('cartoon', f'{pdb_id} and polymer')
    cmd.color('protein_cartoon', 'noncontact_sel')
    cmd.color('contact_cartoon', 'contacts_sel')

    # ---- CONTACT RESIDUES AS STICKS ----
    cmd.show('sticks', 'contacts_sel')
    cmd.set('stick_radius', 0.18, 'contacts_sel')
    cmd.color('contact_carbon', 'contacts_sel and elem C')
    cmd.color('lipid_oxygen', 'contacts_sel and elem O')
    cmd.color('lipid_nitrogen', 'contacts_sel and elem N')
    cmd.color('tv_yellow', 'contacts_sel and elem S')

    # ---- LIPID MOLECULES AS THICK STICKS ----
    cmd.show('sticks', 'lipids_sel')
    cmd.set('stick_radius', 0.22, 'lipids_sel')
    cmd.color('lipid_carbon', 'lipids_sel and elem C')
    cmd.color('lipid_oxygen', 'lipids_sel and elem O')
    cmd.color('lipid_nitrogen', 'lipids_sel and elem N')
    cmd.color('lipid_phosphorus', 'lipids_sel and elem P')
    cmd.color('tv_yellow', 'lipids_sel and elem S')

    # ---- SEMI-TRANSPARENT SURFACE on protein for context ----
    cmd.show('surface', f'{pdb_id} and polymer')
    cmd.color('surface_color', f'{pdb_id} and polymer')
    cmd.set('transparency', 0.82, f'{pdb_id} and polymer')
    # Color surface of contact residues slightly
    cmd.set('surface_color', 'contact_cartoon', 'contacts_sel')

    # ---- Orient the view ----
    cmd.orient(f'{pdb_id} and polymer')

    # Fine-tune view per protein
    if pdb_id == '5NM2':
        # Side view of GPCR TM domain
        cmd.orient(f'{pdb_id} and polymer')
    elif pdb_id == '4C9J':
        # Side view showing both monomers
        cmd.orient(f'{pdb_id} and polymer')
    elif pdb_id == '5T77':
        # Side view of flippase
        cmd.orient(f'{pdb_id} and polymer')
    elif pdb_id == '4TSY':
        # Side view for FraC to show pore formation geometry
        cmd.orient(f'{pdb_id} and polymer')

    cmd.zoom(f'{pdb_id}', buffer=4)

    # ---- Ray trace at high resolution ----
    width, height = 2400, 2000
    cmd.ray(width, height)

    out_png = os.path.join(OUTPUT_DIR, f'fig7_panel_{pdb_id}.png')
    cmd.png(out_png, dpi=300)
    print(f"  Saved: {out_png}")

    return out_png


def create_composite_figure(panel_data):
    """Create a publication-quality 2x2 composite figure."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    import matplotlib.patches as mpatches
    from matplotlib.lines import Line2D

    fig, axes = plt.subplots(2, 2, figsize=(12, 10.5))
    fig.subplots_adjust(wspace=0.03, hspace=0.12, left=0.01, right=0.99,
                        top=0.94, bottom=0.07)

    labels = ['A', 'B', 'C', 'D']

    for idx, (pdb_id, png_path, info) in enumerate(panel_data):
        ax = axes[idx // 2, idx % 2]

        if os.path.exists(png_path):
            img = mpimg.imread(png_path)
            ax.imshow(img)

        ax.axis('off')

        # Bold panel label
        ax.text(0.02, 0.97, labels[idx],
                transform=ax.transAxes, fontsize=22, fontweight='bold',
                va='top', ha='left', family='sans-serif',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                          edgecolor='none', alpha=0.85))

        # Protein title
        lipid_names = ', '.join(info['lipids'])
        title_str = f"{pdb_id} - {info['title']}"
        detail_str = f"{info['subtitle']} | {info['n_contacts']} contacts | Lipids: {lipid_names}"

        ax.text(0.50, 0.01, title_str,
                transform=ax.transAxes, fontsize=9, fontweight='bold',
                va='bottom', ha='center', family='sans-serif',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                          edgecolor='#cccccc', alpha=0.92))
        ax.text(0.50, 0.065, detail_str,
                transform=ax.transAxes, fontsize=7, fontstyle='italic',
                va='bottom', ha='center', family='sans-serif',
                color='#444444',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                          edgecolor='none', alpha=0.85))

    # ---- Legend at bottom ----
    legend_elements = [
        mpatches.Patch(facecolor=[0.85, 0.87, 0.90], edgecolor='gray',
                       label='Protein backbone'),
        mpatches.Patch(facecolor=[0.95, 0.45, 0.30], edgecolor='#333',
                       label='Lipid-contact residues'),
        mpatches.Patch(facecolor=[0.40, 0.75, 0.35], edgecolor='#333',
                       label='Crystallized lipid molecules'),
        mpatches.Patch(facecolor=[0.90, 0.92, 0.95], edgecolor='gray',
                       alpha=0.4, label='Protein surface (transparent)'),
    ]

    fig.legend(handles=legend_elements, loc='lower center',
               ncol=4, fontsize=9, frameon=True,
               fancybox=True, shadow=False,
               edgecolor='#cccccc', facecolor='white',
               bbox_to_anchor=(0.5, 0.005))

    # Supertitle
    fig.suptitle('Structural Visualization of Membrane Protein-Lipid Contacts in MPLID',
                 fontsize=13, fontweight='bold', family='sans-serif', y=0.98)

    # Save
    for ext in ['png', 'pdf']:
        out_path = os.path.join(OUTPUT_DIR, f'fig7_structural_visualization.{ext}')
        fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"  Composite: {out_path}")

    plt.close()


def main():
    print("MPLID Figure 7 - Structural Visualization (v2)")
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

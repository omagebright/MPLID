#!/usr/bin/env python3
"""
PyMOL rendering script for MPLID manuscript Figure 7.
Creates publication-quality structural visualizations of membrane proteins
with their experimentally determined lipid contacts.

Produces individual panels and a combined multi-panel figure.

Usage:
    pymol -cq render_structures.py
    OR
    python3 render_structures.py  (uses PyMOL Python API)
"""

import sys
import os

# Add pymol to path if needed
try:
    import pymol
    from pymol import cmd, stored
except ImportError:
    print("PyMOL not available. Run with: pymol -cq render_structures.py")
    sys.exit(1)

# Initialize PyMOL in quiet mode
pymol.finish_launching(['pymol', '-cq'])

# Paths
STRUCT_DIR = '/HOMENFS/60238528073/Art/Exo_Database/Lipid_membrane_protein_interfaces/MPLID_paper/manuscript/figures/structures'
OUTPUT_DIR = '/HOMENFS/60238528073/Art/Exo_Database/Lipid_membrane_protein_interfaces/MPLID_paper/manuscript/figures'

# Known lipid residue names in our structures
LIPID_NAMES = {
    'OLA', 'OLB', 'OLC', 'CLR', 'CDL', 'HXJ', 'PLM', 'PEE', 'PGV',
    'PCW', 'POV', 'LDA', 'LMT', 'BOG', 'DPC', 'SPH', 'MYR', 'STE',
    'ARA', 'DHA', 'Y01', 'CHD', 'DMU', 'S1P', 'LHG'
}

# Water and buffer molecules to hide
NON_LIPID_HETATM = {'HOH', 'SO4', 'CXT', 'NA', 'CL', 'CA', 'ZN', 'MG',
                     '1PE', 'PEG', 'GOL', 'EDO', 'ACT', 'HTO', 'ZMA'}

# Protein selections with contact residues and lipid types
PROTEINS = {
    '5NM2': {
        'title': 'A2A Adenosine Receptor (GPCR)',
        'short_title': 'A$_2$$_A$AR (GPCR)',
        'contacts': {
            'A': [3, 5, 6, 7, 9, 11, 19, 22, 23, 28, 43, 47, 57, 61, 65, 68,
                  71, 72, 73, 76, 97, 101, 118, 122, 126, 130, 134, 141, 244,
                  247, 248, 251, 252, 254, 255, 258, 259, 262, 263, 265, 268]
        },
        'lipids': ['OLA', 'OLB', 'OLC', 'CLR'],
        'n_contacts': 41,
        'n_residues': 388,
        'view_preset': 'membrane_side'
    },
    '4C9J': {
        'title': 'Mitochondrial ADP/ATP Carrier',
        'short_title': 'AAC3 (Transporter)',
        'contacts': {
            'A': [58, 59, 60, 77, 78, 80, 161, 176, 180, 181, 182, 183, 184,
                  256, 258, 276, 278],
            'B': [58, 59, 77, 78, 80, 161, 162, 176, 181, 182, 183, 256, 258,
                  275, 276, 278]
        },
        'lipids': ['CDL'],
        'n_contacts': 33,
        'n_residues': 544,
        'view_preset': 'membrane_side'
    },
    '5T77': {
        'title': 'MurJ MOP Flippase',
        'short_title': 'MurJ (Flippase)',
        'contacts': {
            'A': [11, 20, 23, 85, 91, 94, 95, 100, 101, 153, 157, 164, 168,
                  179, 183, 195, 206, 207, 208, 209, 216, 224, 268, 269, 303,
                  307, 332, 336, 368, 404, 408, 437, 442]
        },
        'lipids': ['OLC', 'OLB'],
        'n_contacts': 33,
        'n_residues': 467,
        'view_preset': 'membrane_side'
    },
    '4TSY': {
        'title': 'Fragaceatoxin C (FraC)',
        'short_title': 'FraC (Pore-forming)',
        'contacts': {
            'A': [84, 85, 107, 108, 113, 167, 168],
            'B': [78, 84, 85, 107, 108, 167, 168],
            'C': [78, 79, 84, 85, 107, 167, 168],
            'D': [78, 79, 84, 85, 107, 167, 168]
        },
        'lipids': ['HXJ'],
        'n_contacts': 28,
        'n_residues': 704,
        'view_preset': 'membrane_top'
    }
}


def setup_global_settings():
    """Configure PyMOL for publication-quality rendering."""
    cmd.set('ray_opaque_background', 1)
    cmd.set('ray_shadows', 'off')
    cmd.set('antialias', 2)
    cmd.set('ray_trace_mode', 1)
    cmd.set('ray_trace_gain', 0.012)
    cmd.set('ray_trace_color', 'black')
    cmd.set('orthoscopic', 1)
    cmd.set('depth_cue', 0)
    cmd.set('fog', 0)
    cmd.set('spec_reflect', 0.3)
    cmd.set('spec_power', 200)
    cmd.set('ambient', 0.35)
    cmd.set('direct', 0.7)
    # Cartoon settings
    cmd.set('cartoon_fancy_helices', 1)
    cmd.set('cartoon_highlight_color', -1)
    cmd.set('cartoon_smooth_loops', 1)
    cmd.set('cartoon_oval_length', 1.2)
    cmd.set('cartoon_oval_width', 0.25)
    cmd.set('cartoon_loop_radius', 0.2)
    cmd.set('cartoon_tube_radius', 0.3)
    # Stick settings
    cmd.set('stick_radius', 0.15)
    cmd.set('stick_ball', 0)
    # Surface settings
    cmd.set('surface_quality', 1)
    cmd.set('transparency', 0.0)
    # Label settings
    cmd.set('label_size', -0.5)
    cmd.set('label_font_id', 7)
    cmd.set('label_color', 'black')


def build_contact_selection(pdb_id, contacts_dict):
    """Build a PyMOL selection string for contact residues."""
    parts = []
    for chain, residues in contacts_dict.items():
        resi_str = '+'.join(str(r) for r in residues)
        parts.append(f'(chain {chain} and resi {resi_str})')
    return ' or '.join(parts)


def render_protein(pdb_id, info, panel_label=''):
    """Render a single protein with its lipid contacts."""
    cmd.reinitialize()
    setup_global_settings()
    cmd.bg_color('white')

    pdb_file = os.path.join(STRUCT_DIR, f'{pdb_id}.pdb')
    if not os.path.exists(pdb_file):
        print(f"ERROR: {pdb_file} not found!")
        return

    cmd.load(pdb_file, pdb_id)

    # ---- Hide everything, then build up ----
    cmd.hide('everything')

    # ---- PROTEIN: cartoon representation ----
    cmd.show('cartoon', f'{pdb_id} and polymer')
    cmd.color('gray85', f'{pdb_id} and polymer')

    # ---- LIPIDS: stick representation ----
    # Build lipid selection from known lipid residue names present
    lipid_resnames = info['lipids']
    lipid_sel_parts = [f'resn {ln}' for ln in lipid_resnames]
    lipid_sel = f'{pdb_id} and ({" or ".join(lipid_sel_parts)})'

    cmd.show('sticks', lipid_sel)
    cmd.set('stick_radius', 0.18, lipid_sel)
    # Color lipid carbons in a distinctive green-gold
    cmd.color('palegreen', f'{lipid_sel} and elem C')
    cmd.color('red', f'{lipid_sel} and elem O')
    cmd.color('blue', f'{lipid_sel} and elem N')
    cmd.color('orange', f'{lipid_sel} and elem P')
    cmd.color('yellow', f'{lipid_sel} and elem S')

    # ---- CONTACT RESIDUES: sticks in warm red/orange ----
    contact_sel_str = build_contact_selection(pdb_id, info['contacts'])
    cmd.select('contact_residues', f'{pdb_id} and polymer and ({contact_sel_str})')

    cmd.show('sticks', 'contact_residues')
    cmd.set('stick_radius', 0.15, 'contact_residues')
    cmd.color('tv_red', 'contact_residues and elem C')
    cmd.color('red', 'contact_residues and elem O')
    cmd.color('tv_blue', 'contact_residues and elem N')
    cmd.color('yellow', 'contact_residues and elem S')

    # Also color the cartoon backbone of contact residues
    cmd.set('cartoon_color', 'salmon', 'contact_residues')

    # ---- Remove water and non-lipid HETATM ----
    cmd.remove('resn HOH')
    for het in NON_LIPID_HETATM:
        cmd.remove(f'resn {het}')

    # ---- Remove hydrogen atoms for clarity ----
    cmd.remove('hydrogens')

    # ---- Orient the view ----
    cmd.orient(f'{pdb_id}')

    if info.get('view_preset') == 'membrane_top' and pdb_id == '4TSY':
        # For FraC, show top-down view of the pore
        cmd.orient(f'{pdb_id}')
        cmd.turn('x', 90)

    cmd.zoom(f'{pdb_id}', buffer=5)

    # ---- Ray trace and save ----
    width, height = 1200, 1000
    cmd.ray(width, height)

    # Save individual panel
    out_png = os.path.join(OUTPUT_DIR, f'fig7_panel_{pdb_id}.png')
    cmd.png(out_png, dpi=300)
    print(f"Saved: {out_png}")

    return out_png


def render_all_proteins():
    """Render all selected proteins."""
    panel_files = []
    labels = ['A', 'B', 'C', 'D']

    for i, (pdb_id, info) in enumerate(PROTEINS.items()):
        print(f"\n{'='*60}")
        print(f"Rendering {pdb_id}: {info['title']}")
        print(f"  Contacts: {info['n_contacts']} / {info['n_residues']} residues")
        print(f"  Lipids: {', '.join(info['lipids'])}")
        print(f"{'='*60}")

        out_file = render_protein(pdb_id, info, panel_label=labels[i])
        if out_file:
            panel_files.append((pdb_id, out_file, info))

    return panel_files


def create_composite_figure(panel_files):
    """Create a 2x2 composite figure from individual panels using matplotlib."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    from matplotlib.offsetbox import AnchoredText

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.subplots_adjust(wspace=0.02, hspace=0.08, left=0.02, right=0.98,
                        top=0.95, bottom=0.02)

    labels = ['A', 'B', 'C', 'D']

    for idx, (pdb_id, png_file, info) in enumerate(panel_files):
        ax = axes[idx // 2, idx % 2]

        if os.path.exists(png_file):
            img = mpimg.imread(png_file)
            ax.imshow(img)

        ax.axis('off')

        # Panel label (A, B, C, D) in top-left
        ax.text(0.02, 0.98, f'{labels[idx]}',
                transform=ax.transAxes, fontsize=20, fontweight='bold',
                va='top', ha='left',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                          edgecolor='none', alpha=0.8))

        # Protein info below panel label
        lipid_str = ', '.join(info['lipids'])
        caption = (f"{pdb_id}: {info['short_title']}\n"
                   f"{info['n_contacts']} contacts / {info['n_residues']} residues\n"
                   f"Lipids: {lipid_str}")
        ax.text(0.02, 0.02, caption,
                transform=ax.transAxes, fontsize=8,
                va='bottom', ha='left',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                          edgecolor='gray', alpha=0.9),
                family='sans-serif')

    # Add color legend at bottom
    legend_text = (r'$\mathbf{Gray}$: protein backbone    '
                   r'$\mathbf{Red}$: lipid-contact residues    '
                   r'$\mathbf{Green}$: crystallized lipids')
    fig.text(0.5, 0.005, legend_text, ha='center', va='bottom', fontsize=10,
             family='sans-serif')

    # Save composite
    composite_png = os.path.join(OUTPUT_DIR, 'fig7_structural_visualization.png')
    composite_pdf = os.path.join(OUTPUT_DIR, 'fig7_structural_visualization.pdf')

    fig.savefig(composite_png, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(composite_pdf, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"\nComposite saved: {composite_png}")
    print(f"Composite saved: {composite_pdf}")


if __name__ == '__main__':
    print("MPLID Structural Visualization Rendering")
    print("=" * 60)

    # Render individual panels
    panel_files = render_all_proteins()

    # Create composite figure
    if panel_files:
        create_composite_figure(panel_files)

    print("\nDone!")

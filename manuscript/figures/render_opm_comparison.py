#!/usr/bin/env python3
"""
PyMOL rendering script for MPLID Figure 8:
Side-by-side comparison of MPLID experimental lipid contacts vs. OPM
computational membrane boundaries for the A2A adenosine receptor (5NM2).

Panel A: MPLID experimental contacts (red residues + green lipids)
Panel B: OPM membrane zone (blue = in membrane, gray = outside)
Panel C: Agreement/disagreement overlay

Usage:
    python3 render_opm_comparison.py
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

STRUCT_DIR = '/HOMENFS/60238528073/Art/Exo_Database/Lipid_membrane_protein_interfaces/MPLID_paper/manuscript/figures/structures'
OPM_DIR = '/HOMENFS/60238528073/Art/Exo_Database/Lipid_membrane_protein_interfaces/data/raw/opm'
OUTPUT_DIR = '/HOMENFS/60238528073/Art/Exo_Database/Lipid_membrane_protein_interfaces/MPLID_paper/manuscript/figures'

# 5NM2 contact residues from MPLID
CONTACTS_5NM2 = [3, 5, 6, 7, 9, 11, 19, 22, 23, 28, 43, 47, 57, 61, 65, 68,
                 71, 72, 73, 76, 97, 101, 118, 122, 126, 130, 134, 141, 244,
                 247, 248, 251, 252, 254, 255, 258, 259, 262, 263, 265, 268]

LIPID_RESNAMES = ['OLA', 'OLB', 'OLC', 'CLR']

REMOVE_HETATM = {'HOH', 'SO4', 'CXT', 'NA', 'CL', 'CA', 'ZN', 'MG',
                  '1PE', 'PEG', 'GOL', 'EDO', 'ACT', 'HTO', 'ZMA'}

# OPM membrane boundaries for 5nm2_opm.pdb
MEMBRANE_Z_UPPER = 15.70
MEMBRANE_Z_LOWER = -15.70


def setup_rendering():
    """Configure PyMOL for publication quality."""
    cmd.bg_color('white')
    cmd.set('ray_opaque_background', 1)
    cmd.set('ray_trace_mode', 1)
    cmd.set('ray_trace_gain', 0.008)
    cmd.set('ray_trace_color', 'black')
    cmd.set('ray_shadows', 'off')
    cmd.set('antialias', 2)
    cmd.set('orthoscopic', 1)
    cmd.set('depth_cue', 0)
    cmd.set('fog', 0)
    cmd.set('ambient', 0.4)
    cmd.set('direct', 0.6)
    cmd.set('reflect', 0.3)
    cmd.set('spec_reflect', 0.2)
    cmd.set('spec_power', 200)
    cmd.set('cartoon_fancy_helices', 1)
    cmd.set('cartoon_smooth_loops', 1)
    cmd.set('cartoon_oval_length', 1.2)
    cmd.set('cartoon_oval_width', 0.25)
    cmd.set('cartoon_discrete_colors', 1)
    cmd.set('stick_quality', 15)


def get_stored_view():
    """Get the current view matrix to reuse across panels."""
    return cmd.get_view()


def render_panel_a(saved_view=None):
    """Panel A: MPLID experimental lipid contacts."""
    cmd.reinitialize()
    setup_rendering()

    pdb_file = os.path.join(STRUCT_DIR, '5NM2.pdb')
    cmd.load(pdb_file, 'protein')
    cmd.hide('everything')

    # Clean up
    for het in REMOVE_HETATM:
        cmd.remove(f'resn {het}')
    cmd.remove('hydrogens')
    # Remove fusion protein region
    cmd.remove('protein and polymer and chain A and (resi 0 or resi 311-9999)')
    cmd.remove('protein and polymer and not chain A')

    # Protein cartoon
    cmd.show('cartoon', 'protein and polymer')
    cmd.color('gray85', 'protein and polymer')

    # Contact residues
    resi_str = '+'.join(str(r) for r in CONTACTS_5NM2)
    cmd.select('contacts', f'protein and polymer and chain A and resi {resi_str}')
    cmd.show('sticks', 'contacts')
    cmd.set('stick_radius', 0.18, 'contacts')

    # Color contacts
    cmd.set_color('contact_c', [0.90, 0.25, 0.20])
    cmd.set_color('contact_cartoon_c', [0.95, 0.45, 0.30])
    cmd.color('contact_c', 'contacts and elem C')
    cmd.color('red', 'contacts and elem O')
    cmd.color('tv_blue', 'contacts and elem N')
    cmd.set('cartoon_color', 'contact_cartoon_c', 'contacts')

    # Lipids
    lipid_sel = ' or '.join([f'resn {ln}' for ln in LIPID_RESNAMES])
    cmd.select('lipids', f'protein and ({lipid_sel})')
    cmd.show('sticks', 'lipids')
    cmd.set('stick_radius', 0.22, 'lipids')
    cmd.set_color('lipid_c', [0.40, 0.75, 0.35])
    cmd.color('lipid_c', 'lipids and elem C')
    cmd.color('red', 'lipids and elem O')
    cmd.color('tv_blue', 'lipids and elem N')
    cmd.set_color('lipid_p', [1.00, 0.60, 0.00])
    cmd.color('lipid_p', 'lipids and elem P')

    # Orient
    cmd.orient('protein and polymer')
    cmd.zoom('protein', buffer=5)

    if saved_view:
        cmd.set_view(saved_view)

    view = cmd.get_view()

    cmd.ray(1600, 1400)
    out = os.path.join(OUTPUT_DIR, 'fig8_panel_A_mplid.png')
    cmd.png(out, dpi=300)
    print(f"  Panel A: {out}")
    return view


def render_panel_b(saved_view):
    """Panel B: OPM membrane zone labeling."""
    cmd.reinitialize()
    setup_rendering()

    opm_file = os.path.join(OPM_DIR, '5nm2_opm.pdb')
    cmd.load(opm_file, 'opm_protein')
    cmd.hide('everything')

    # Remove HETATM except DUM
    for het in REMOVE_HETATM:
        cmd.remove(f'resn {het}')
    cmd.remove('hydrogens')

    # Show protein cartoon
    cmd.show('cartoon', 'opm_protein and polymer')

    # Color based on Z-coordinate (membrane zone)
    # We need to iterate over residues and color based on CA Z
    stored.membrane_resi = []
    stored.outside_resi = []

    cmd.iterate_state(1, 'opm_protein and polymer and name CA',
                      f'stored.membrane_resi.append(resi) if z >= {MEMBRANE_Z_LOWER} and z <= {MEMBRANE_Z_UPPER} else stored.outside_resi.append(resi)')

    membrane_resi = list(set(stored.membrane_resi))
    outside_resi = list(set(stored.outside_resi))

    if membrane_resi:
        mem_sel = '+'.join(membrane_resi)
        cmd.select('in_membrane', f'opm_protein and polymer and resi {mem_sel}')
        cmd.set_color('membrane_blue', [0.30, 0.55, 0.85])
        cmd.color('membrane_blue', 'in_membrane')
        cmd.show('sticks', 'in_membrane')
        cmd.set('stick_radius', 0.12, 'in_membrane')
        cmd.color('membrane_blue', 'in_membrane and elem C')
        cmd.set('cartoon_color', 'membrane_blue', 'in_membrane')

    if outside_resi:
        out_sel = '+'.join(outside_resi)
        cmd.select('outside_membrane', f'opm_protein and polymer and resi {out_sel}')
        cmd.color('gray85', 'outside_membrane')

    # Show DUM atoms as membrane boundary planes (thin dots)
    cmd.show('nb_spheres', 'resn DUM')
    cmd.set('sphere_scale', 0.3, 'resn DUM')
    cmd.set_color('dum_color', [0.75, 0.85, 1.0])
    cmd.color('dum_color', 'resn DUM')
    cmd.set('sphere_transparency', 0.6, 'resn DUM')

    # Apply same view as Panel A
    cmd.set_view(saved_view)

    cmd.ray(1600, 1400)
    out = os.path.join(OUTPUT_DIR, 'fig8_panel_B_opm.png')
    cmd.png(out, dpi=300)
    print(f"  Panel B: {out}")


def render_panel_c(saved_view):
    """Panel C: Agreement/disagreement overlay using OPM structure."""
    cmd.reinitialize()
    setup_rendering()

    opm_file = os.path.join(OPM_DIR, '5nm2_opm.pdb')
    cmd.load(opm_file, 'overlay')
    cmd.hide('everything')

    for het in REMOVE_HETATM:
        cmd.remove(f'resn {het}')
    cmd.remove('hydrogens')
    cmd.remove('resn DUM')

    cmd.show('cartoon', 'overlay and polymer')
    cmd.color('gray85', 'overlay and polymer')

    # Identify membrane residues from OPM Z coordinates
    stored.membrane_resi_c = []
    cmd.iterate_state(1, 'overlay and polymer and name CA',
                      f'stored.membrane_resi_c.append(resi) if z >= {MEMBRANE_Z_LOWER} and z <= {MEMBRANE_Z_UPPER} else None')
    opm_membrane_set = set(stored.membrane_resi_c)

    # MPLID contact set (as strings for comparison)
    mplid_contact_set = set(str(r) for r in CONTACTS_5NM2)

    # Categories:
    # True Positive (TP): MPLID contact AND OPM membrane
    # MPLID-only: MPLID contact but NOT in OPM membrane
    # OPM-only: In OPM membrane but NOT MPLID contact
    # Neither: not contact, not in membrane

    tp = mplid_contact_set & opm_membrane_set
    mplid_only = mplid_contact_set - opm_membrane_set
    opm_only = opm_membrane_set - mplid_contact_set

    print(f"\n  Agreement analysis:")
    print(f"    Both (TP): {len(tp)}")
    print(f"    MPLID-only: {len(mplid_only)} (contacts outside OPM membrane)")
    print(f"    OPM-only: {len(opm_only)} (in membrane but no lipid contact)")

    # Color categories
    cmd.set_color('agree_color', [0.20, 0.70, 0.30])      # Green = agreement
    cmd.set_color('mplid_only_color', [0.90, 0.25, 0.20]) # Red = MPLID-only
    cmd.set_color('opm_only_color', [0.30, 0.55, 0.85])   # Blue = OPM-only

    if tp:
        tp_sel = '+'.join(tp)
        cmd.select('agree_sel', f'overlay and polymer and resi {tp_sel}')
        cmd.color('agree_color', 'agree_sel')
        cmd.set('cartoon_color', 'agree_color', 'agree_sel')
        cmd.show('sticks', 'agree_sel')
        cmd.set('stick_radius', 0.15, 'agree_sel')
        cmd.color('agree_color', 'agree_sel and elem C')

    if mplid_only:
        mo_sel = '+'.join(mplid_only)
        cmd.select('mplid_only_sel', f'overlay and polymer and resi {mo_sel}')
        cmd.color('mplid_only_color', 'mplid_only_sel')
        cmd.set('cartoon_color', 'mplid_only_color', 'mplid_only_sel')
        cmd.show('sticks', 'mplid_only_sel')
        cmd.set('stick_radius', 0.15, 'mplid_only_sel')
        cmd.color('mplid_only_color', 'mplid_only_sel and elem C')

    if opm_only:
        oo_sel = '+'.join(opm_only)
        cmd.select('opm_only_sel', f'overlay and polymer and resi {oo_sel}')
        cmd.set('cartoon_color', 'opm_only_color', 'opm_only_sel')

    cmd.set_view(saved_view)

    cmd.ray(1600, 1400)
    out = os.path.join(OUTPUT_DIR, 'fig8_panel_C_overlay.png')
    cmd.png(out, dpi=300)
    print(f"  Panel C: {out}")


def create_composite():
    """Assemble 1x3 composite figure."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    import matplotlib.patches as mpatches

    panels = [
        ('A', 'fig8_panel_A_mplid.png', 'MPLID\n(Experimental Lipid Contacts)'),
        ('B', 'fig8_panel_B_opm.png', 'OPM\n(Computational Membrane Zone)'),
        ('C', 'fig8_panel_C_overlay.png', 'Agreement Overlay'),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(16, 6.5))
    fig.subplots_adjust(wspace=0.03, left=0.01, right=0.99, top=0.88, bottom=0.12)

    for idx, (label, fname, subtitle) in enumerate(panels):
        ax = axes[idx]
        img_path = os.path.join(OUTPUT_DIR, fname)
        if os.path.exists(img_path):
            img = mpimg.imread(img_path)
            ax.imshow(img)
        ax.axis('off')

        # Panel label
        ax.text(0.03, 0.97, label,
                transform=ax.transAxes, fontsize=20, fontweight='bold',
                va='top', ha='left', family='sans-serif',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                          edgecolor='none', alpha=0.85))
        # Subtitle
        ax.set_title(subtitle, fontsize=10, fontweight='bold',
                     family='sans-serif', pad=8)

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=[0.95, 0.45, 0.30], edgecolor='#333',
                       label='MPLID contact residues (Panel A)'),
        mpatches.Patch(facecolor=[0.40, 0.75, 0.35], edgecolor='#333',
                       label='Crystallized lipids (Panel A)'),
        mpatches.Patch(facecolor=[0.30, 0.55, 0.85], edgecolor='#333',
                       label='OPM membrane zone (Panel B)'),
        mpatches.Patch(facecolor=[0.20, 0.70, 0.30], edgecolor='#333',
                       label='Agreement: both methods (Panel C)'),
        mpatches.Patch(facecolor=[0.90, 0.25, 0.20], edgecolor='#333',
                       label='MPLID-only: contact outside OPM zone (Panel C)'),
        mpatches.Patch(facecolor=[0.30, 0.55, 0.85], edgecolor='#333',
                       alpha=0.6, label='OPM-only: in zone, no contact (Panel C)'),
    ]

    fig.legend(handles=legend_elements, loc='lower center',
               ncol=3, fontsize=7.5, frameon=True,
               fancybox=True, edgecolor='#cccccc', facecolor='white',
               bbox_to_anchor=(0.5, -0.01))

    fig.suptitle(r'Comparison: MPLID Experimental Contacts vs. OPM Membrane Zone (A$_2$$_A$AR, PDB: 5NM2)',
                 fontsize=12, fontweight='bold', family='sans-serif', y=0.97)

    for ext in ['png', 'pdf']:
        out = os.path.join(OUTPUT_DIR, f'fig8_opm_comparison.{ext}')
        fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"  Composite: {out}")

    plt.close()


def main():
    print("MPLID Figure 8 - OPM vs. MPLID Comparison")
    print("=" * 60)

    # Panel A: MPLID contacts (also captures view for reuse)
    print("\nRendering Panel A (MPLID contacts)...")
    view = render_panel_a()

    # Panel B: OPM membrane zone (same view)
    print("Rendering Panel B (OPM membrane zone)...")
    render_panel_b(view)

    # Panel C: Agreement overlay (same view)
    print("Rendering Panel C (Agreement overlay)...")
    render_panel_c(view)

    # Composite
    print("\nCreating composite figure...")
    create_composite()

    print("\nDone!")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
PyMOL rendering script for MPLID Figure 8 (v2):
Side-by-side comparison using the OPM-aligned structure for consistent orientation.

Panel A: MPLID experimental contacts (from original PDB lipids, mapped onto OPM coords)
Panel B: OPM membrane zone (blue = in membrane, gray = outside)
Panel C: Agreement/disagreement overlay

All panels use the OPM-aligned coordinate frame for consistent orientation.
"""

import sys
import os

try:
    import pymol
    from pymol import cmd, stored
except ImportError:
    sys.exit(1)

pymol.finish_launching(['pymol', '-cq'])

STRUCT_DIR = '/HOMENFS/60238528073/Art/Exo_Database/Lipid_membrane_protein_interfaces/MPLID_paper/manuscript/figures/structures'
OPM_DIR = '/HOMENFS/60238528073/Art/Exo_Database/Lipid_membrane_protein_interfaces/data/raw/opm'
OUTPUT_DIR = '/HOMENFS/60238528073/Art/Exo_Database/Lipid_membrane_protein_interfaces/MPLID_paper/manuscript/figures'

# MPLID contact residues
CONTACTS = [3, 5, 6, 7, 9, 11, 19, 22, 23, 28, 43, 47, 57, 61, 65, 68,
            71, 72, 73, 76, 97, 101, 118, 122, 126, 130, 134, 141, 244,
            247, 248, 251, 252, 254, 255, 258, 259, 262, 263, 265, 268]
CONTACTS_SET = set(str(r) for r in CONTACTS)

MEMBRANE_Z_UPPER = 15.70
MEMBRANE_Z_LOWER = -15.70

REMOVE_HETATM = {'HOH', 'SO4', 'CXT', 'NA', 'CL', 'CA', 'ZN', 'MG',
                  '1PE', 'PEG', 'GOL', 'EDO', 'ACT', 'HTO', 'ZMA'}


def setup():
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

    # Custom colors
    cmd.set_color('prot_gray', [0.85, 0.87, 0.90])
    cmd.set_color('contact_red', [0.90, 0.25, 0.20])
    cmd.set_color('contact_cartoon', [0.95, 0.45, 0.30])
    cmd.set_color('lipid_green', [0.40, 0.75, 0.35])
    cmd.set_color('lipid_p_orange', [1.00, 0.60, 0.00])
    cmd.set_color('membrane_blue', [0.30, 0.55, 0.85])
    cmd.set_color('membrane_cartoon', [0.35, 0.60, 0.88])
    cmd.set_color('agree_green', [0.20, 0.70, 0.30])
    cmd.set_color('mplid_only_red', [0.90, 0.25, 0.20])
    cmd.set_color('opm_only_blue', [0.50, 0.70, 0.90])
    cmd.set_color('dum_light', [0.75, 0.85, 1.0])


def get_membrane_residues():
    """Get the set of residues within the OPM membrane zone."""
    stored.mem_resi = []
    cmd.iterate_state(1, 'mol and polymer and name CA',
                      f'stored.mem_resi.append(resi) if z >= {MEMBRANE_Z_LOWER} and z <= {MEMBRANE_Z_UPPER} else None')
    return set(stored.mem_resi)


def standard_view():
    """Set a consistent side-on view of the membrane protein."""
    cmd.orient('mol and polymer')
    # OPM aligns membrane plane in XY, so side view is the default orient
    # Just ensure we're looking from a good angle
    cmd.turn('y', 0)
    cmd.zoom('mol and polymer', buffer=8)
    return cmd.get_view()


def render_panel_a():
    """Panel A: MPLID experimental lipid contacts on OPM-aligned structure."""
    cmd.reinitialize()
    setup()

    opm_file = os.path.join(OPM_DIR, '5nm2_opm.pdb')
    orig_file = os.path.join(STRUCT_DIR, '5NM2.pdb')

    # Load OPM structure for consistent orientation
    cmd.load(opm_file, 'mol')
    cmd.hide('everything')

    for het in REMOVE_HETATM:
        cmd.remove(f'resn {het}')
    cmd.remove('hydrogens')
    cmd.remove('resn DUM')

    # Show protein cartoon
    cmd.show('cartoon', 'mol and polymer')
    cmd.color('prot_gray', 'mol and polymer')

    # Highlight contact residues
    resi_str = '+'.join(str(r) for r in CONTACTS)
    cmd.select('contacts', f'mol and polymer and resi {resi_str}')
    cmd.color('contact_cartoon', 'contacts')
    cmd.set('cartoon_color', 'contact_cartoon', 'contacts')
    cmd.show('sticks', 'contacts')
    cmd.set('stick_radius', 0.18, 'contacts')
    cmd.color('contact_red', 'contacts and elem C')
    cmd.color('red', 'contacts and elem O')
    cmd.color('tv_blue', 'contacts and elem N')

    # Load original PDB to get lipid molecules, then align
    cmd.load(orig_file, 'orig')
    cmd.hide('everything', 'orig')
    cmd.align('orig and polymer', 'mol and polymer')

    # Remove all non-lipid HETATM from orig
    for het in REMOVE_HETATM:
        cmd.remove(f'orig and resn {het}')
    cmd.remove('orig and hydrogens')
    # Remove any other non-lipid HETATM (ligands, ions, etc.)
    cmd.remove('orig and not (polymer or resn OLA or resn OLB or resn OLC or resn CLR)')

    # Show lipids from original structure
    lipid_sel = 'orig and (resn OLA or resn OLB or resn OLC or resn CLR)'
    cmd.show('sticks', lipid_sel)
    cmd.set('stick_radius', 0.22, lipid_sel)
    cmd.color('lipid_green', f'{lipid_sel} and elem C')
    cmd.color('red', f'{lipid_sel} and elem O')
    cmd.color('tv_blue', f'{lipid_sel} and elem N')
    cmd.color('lipid_p_orange', f'{lipid_sel} and elem P')

    # Hide original protein (keep only lipids from orig)
    cmd.hide('cartoon', 'orig and polymer')

    view = standard_view()

    cmd.ray(1600, 1400)
    out = os.path.join(OUTPUT_DIR, 'fig8_panel_A_mplid.png')
    cmd.png(out, dpi=300)
    print(f"  Panel A saved: {out}")
    return view


def render_panel_b(view):
    """Panel B: OPM membrane zone labeling."""
    cmd.reinitialize()
    setup()

    opm_file = os.path.join(OPM_DIR, '5nm2_opm.pdb')
    cmd.load(opm_file, 'mol')
    cmd.hide('everything')

    for het in REMOVE_HETATM:
        cmd.remove(f'resn {het}')
    cmd.remove('hydrogens')

    # Show DUM atoms as membrane boundary
    cmd.show('nb_spheres', 'resn DUM')
    cmd.set('sphere_scale', 0.25, 'resn DUM')
    cmd.color('dum_light', 'resn DUM')
    cmd.set('sphere_transparency', 0.65, 'resn DUM')

    # Show protein cartoon
    cmd.show('cartoon', 'mol and polymer')
    cmd.color('prot_gray', 'mol and polymer')

    # Color membrane zone residues blue
    opm_membrane = get_membrane_residues()

    if opm_membrane:
        mem_sel = '+'.join(opm_membrane)
        cmd.select('in_membrane', f'mol and polymer and resi {mem_sel}')
        cmd.color('membrane_cartoon', 'in_membrane')
        cmd.set('cartoon_color', 'membrane_cartoon', 'in_membrane')

    cmd.set_view(view)

    cmd.ray(1600, 1400)
    out = os.path.join(OUTPUT_DIR, 'fig8_panel_B_opm.png')
    cmd.png(out, dpi=300)
    print(f"  Panel B saved: {out}")


def render_panel_c(view):
    """Panel C: Agreement/disagreement overlay."""
    cmd.reinitialize()
    setup()

    opm_file = os.path.join(OPM_DIR, '5nm2_opm.pdb')
    cmd.load(opm_file, 'mol')
    cmd.hide('everything')

    for het in REMOVE_HETATM:
        cmd.remove(f'resn {het}')
    cmd.remove('hydrogens')
    cmd.remove('resn DUM')

    cmd.show('cartoon', 'mol and polymer')
    cmd.color('prot_gray', 'mol and polymer')

    # Get membrane residues
    opm_membrane = get_membrane_residues()

    # Classification
    tp = CONTACTS_SET & opm_membrane          # Both: agreement
    mplid_only = CONTACTS_SET - opm_membrane  # MPLID contact, outside OPM zone
    opm_only = opm_membrane - CONTACTS_SET    # In OPM zone, no MPLID contact

    print(f"    Agreement (both): {len(tp)}")
    print(f"    MPLID-only: {len(mplid_only)}")
    print(f"    OPM-only: {len(opm_only)}")

    # Color: OPM-only in light blue
    if opm_only:
        oo_sel = '+'.join(opm_only)
        cmd.select('opm_only_sel', f'mol and polymer and resi {oo_sel}')
        cmd.color('opm_only_blue', 'opm_only_sel')
        cmd.set('cartoon_color', 'opm_only_blue', 'opm_only_sel')

    # Color: Agreement in green
    if tp:
        tp_sel = '+'.join(tp)
        cmd.select('agree_sel', f'mol and polymer and resi {tp_sel}')
        cmd.color('agree_green', 'agree_sel')
        cmd.set('cartoon_color', 'agree_green', 'agree_sel')
        cmd.show('sticks', 'agree_sel')
        cmd.set('stick_radius', 0.16, 'agree_sel')
        cmd.color('agree_green', 'agree_sel and elem C')

    # Color: MPLID-only in red (most visually striking)
    if mplid_only:
        mo_sel = '+'.join(mplid_only)
        cmd.select('mplid_only_sel', f'mol and polymer and resi {mo_sel}')
        cmd.color('mplid_only_red', 'mplid_only_sel')
        cmd.set('cartoon_color', 'mplid_only_red', 'mplid_only_sel')
        cmd.show('sticks', 'mplid_only_sel')
        cmd.set('stick_radius', 0.16, 'mplid_only_sel')
        cmd.color('mplid_only_red', 'mplid_only_sel and elem C')

    cmd.set_view(view)

    cmd.ray(1600, 1400)
    out = os.path.join(OUTPUT_DIR, 'fig8_panel_C_overlay.png')
    cmd.png(out, dpi=300)
    print(f"  Panel C saved: {out}")


def create_composite():
    """Assemble the 1x3 composite figure with annotations."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    import matplotlib.patches as mpatches

    fig, axes = plt.subplots(1, 3, figsize=(17, 7))
    fig.subplots_adjust(wspace=0.02, left=0.01, right=0.99, top=0.86, bottom=0.14)

    panels = [
        ('A', 'fig8_panel_A_mplid.png',
         'MPLID: Experimental Lipid Contacts'),
        ('B', 'fig8_panel_B_opm.png',
         'OPM: Computational Membrane Zone'),
        ('C', 'fig8_panel_C_overlay.png',
         'Agreement Overlay'),
    ]

    for idx, (label, fname, title) in enumerate(panels):
        ax = axes[idx]
        path = os.path.join(OUTPUT_DIR, fname)
        if os.path.exists(path):
            img = mpimg.imread(path)
            ax.imshow(img)
        ax.axis('off')

        ax.text(0.03, 0.97, label, transform=ax.transAxes,
                fontsize=22, fontweight='bold', va='top', ha='left',
                family='sans-serif',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                          edgecolor='none', alpha=0.85))

        ax.set_title(title, fontsize=10, fontweight='bold',
                     family='sans-serif', pad=6)

    # Legend - two rows
    legend_elements_top = [
        mpatches.Patch(facecolor=[0.95, 0.45, 0.30], edgecolor='#333',
                       label='MPLID contact residues'),
        mpatches.Patch(facecolor=[0.40, 0.75, 0.35], edgecolor='#333',
                       label='Crystallized lipid molecules'),
        mpatches.Patch(facecolor=[0.35, 0.60, 0.88], edgecolor='#333',
                       label='OPM membrane zone residues'),
    ]

    legend_elements_bottom = [
        mpatches.Patch(facecolor=[0.20, 0.70, 0.30], edgecolor='#333',
                       label='Agreement (contact + in membrane)'),
        mpatches.Patch(facecolor=[0.90, 0.25, 0.20], edgecolor='#333',
                       label='MPLID-only (contact outside membrane)'),
        mpatches.Patch(facecolor=[0.50, 0.70, 0.90], edgecolor='#333',
                       label='OPM-only (in membrane, no contact)'),
    ]

    leg1 = fig.legend(handles=legend_elements_top, loc='lower center',
                      ncol=3, fontsize=8.5, frameon=True,
                      fancybox=True, edgecolor='#cccccc',
                      bbox_to_anchor=(0.5, 0.065))

    fig.legend(handles=legend_elements_bottom, loc='lower center',
               ncol=3, fontsize=8.5, frameon=True,
               fancybox=True, edgecolor='#cccccc',
               bbox_to_anchor=(0.5, 0.005))

    fig.suptitle(
        r'MPLID Experimental Contacts vs. OPM Membrane Zone: A$_{2A}$ Adenosine Receptor (PDB: 5NM2)',
        fontsize=13, fontweight='bold', family='sans-serif', y=0.96)

    for ext in ['png', 'pdf']:
        out = os.path.join(OUTPUT_DIR, f'fig8_opm_comparison.{ext}')
        fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"  Composite: {out}")
    plt.close()


def main():
    print("MPLID Figure 8 - OPM vs. MPLID Comparison (v2)")
    print("=" * 60)

    print("\nPanel A: MPLID contacts...")
    view = render_panel_a()

    print("Panel B: OPM membrane zone...")
    render_panel_b(view)

    print("Panel C: Agreement overlay...")
    render_panel_c(view)

    print("\nComposite figure...")
    create_composite()

    print("\nDone!")


if __name__ == '__main__':
    main()

"""Per-atom and per-bond render data for the PyVista molecule view.

All distances in Bohr (matching pyscf.gto.M defaults).
"""

import numpy as np
from pyscf.data import radii

# CPK-style colors keyed by atomic number. Hex strings, no leading '#'.
# First two periods + common heavies — enough for teaching examples.
# Anything else falls back to FALLBACK_COLOR.
CPK_COLORS = {
    1:  '#f0f0f0',  # H — light grey instead of pure white so it renders against light bg
    2:  '#d9ffff',  # He
    3:  '#cc80ff',  # Li
    4:  '#c2ff00',  # Be
    5:  '#ffb5b5',  # B
    6:  '#404040',  # C — dark grey (canonical CPK)
    7:  '#3050f8',  # N
    8:  '#ff0d0d',  # O
    9:  '#90e050',  # F
    10: '#b3e3f5',  # Ne
    11: '#ab5cf2',  # Na
    12: '#8aff00',  # Mg
    13: '#bfa6a6',  # Al
    14: '#f0c8a0',  # Si
    15: '#ff8000',  # P
    16: '#ffff30',  # S
    17: '#1ff01f',  # Cl
    18: '#80d1e3',  # Ar
    19: '#8f40d4',  # K
    20: '#3dff00',  # Ca
    35: '#a62929',  # Br
    53: '#940094',  # I
}
FALLBACK_COLOR = '#ff80ff'  # bright magenta — makes "I forgot to add this element" visible

# Visual radius is a fraction of the covalent radius. Smaller than 1 because
# we draw bonds as sticks; full-covalent spheres would swallow them.
ATOM_RADIUS_SCALE = 0.4
# Atoms count as bonded if their separation is less than this multiple of
# the sum of their covalent radii.
BOND_LENGTH_FACTOR = 1.2
# Stick radius for bond cylinders (Bohr).
BOND_RADIUS = 0.08


def atom_render_data(mol):
    """Return a list of dicts {Z, symbol, pos, radius, color} for each atom in mol."""
    Zs = mol.atom_charges()
    coords = mol.atom_coords()  # Bohr
    out = []
    for i, Z in enumerate(Zs):
        out.append({
            'Z': int(Z),
            'symbol': mol.atom_symbol(i),
            'pos': coords[i],
            'radius': float(radii.COVALENT[Z]) * ATOM_RADIUS_SCALE,
            'color': CPK_COLORS.get(int(Z), FALLBACK_COLOR),
        })
    return out


def detect_bonds(mol):
    """Yield (i, j) atom-index pairs that should be drawn as bonded."""
    Zs = mol.atom_charges()
    coords = mol.atom_coords()
    n = len(Zs)
    for i in range(n):
        for j in range(i + 1, n):
            r_sum = radii.COVALENT[Zs[i]] + radii.COVALENT[Zs[j]]
            d = np.linalg.norm(coords[i] - coords[j])
            if d <= BOND_LENGTH_FACTOR * r_sum:
                yield i, j

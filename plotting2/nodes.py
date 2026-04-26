from ryven.node_env import *
import numpy as np

def get_atom_surface_points(atoms_coords, axes, size):
    theta = np.linspace(0, 2 * np.pi, 20)
    phi = np.linspace(0, np.pi, 10)
    theta, phi = np.meshgrid(theta, phi)
    r = 0.2
    x = r * np.sin(phi) * np.cos(theta)
    y = r * np.sin(phi) * np.sin(theta)
    z = r * np.cos(phi)
    for atom in atoms_coords:
        x_ = x + atom[0]
        y_ = y + atom[1]
        z_ = z + atom[2]
        axes.plot_surface(x_, y_, z_, color='red', alpha=0.8)

def _eval_mo_grid(mol, coeff_1d, bnds, grid_points=50):
    max_, min_ = max(bnds), min(bnds)
    axis = np.linspace(min_, max_, grid_points)
    spacing = (axis[1] - axis[0],) * 3
    X, Y, Z = np.meshgrid(axis, axis, axis, indexing='ij')
    coords = np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T
    ao = mol.eval_gto('GTOval', coords)
    mo = np.einsum('i,ki->k', coeff_1d, ao)
    return mo.reshape(X.shape), spacing, min_

def _safe_marching(grid, level, spacing, origin):
    from skimage import measure
    try:
        verts, faces, _, _ = measure.marching_cubes(grid, level, spacing=spacing)
    except (ValueError, RuntimeError):
        return None, None
    return verts + origin, faces

def get_isosurface(mol, mo_coeff, orbital, iso_val, bnds, beta=False):
    if orbital < 0:
        # TODO: total electron density
        return None
    if mo_coeff.ndim == 3:
        coeff_1d = mo_coeff[1 if beta else 0, :, orbital]
    else:
        coeff_1d = mo_coeff[:, orbital]
    grid, spacing, origin = _eval_mo_grid(mol, coeff_1d, bnds)
    pos_verts, pos_faces = _safe_marching(grid, iso_val, spacing, origin)
    neg_verts, neg_faces = _safe_marching(grid, -iso_val, spacing, origin)
    return (pos_verts, pos_faces, neg_verts, neg_faces)

class LinePlotNode(Node):
    """Generates a line graph"""

    title = 'LinePlot'
    #tags = ['random', 'numbers']
    init_inputs = [NodeInputType(label='1'),
                   NodeInputType(label='2'),
                   NodeInputType(label='3'),
                   NodeInputType(label='4'),
                   NodeInputType(label='5')]

    def inputs_ready(self):
        return all(self.input(i) is not None for i in range(len(self.inputs)))

    def update_event(self, inp=-1):
        #if not self.inputs_ready():
        #    return
        """Inputs can be empty just don't plot those"""

        if self.have_gui():
            self.gui.update_plot()

    def have_gui(self):
        return hasattr(self, 'gui')

class MOPlotNode(Node):
    """Generates an interactive surface"""

    title = "Molecular Orbital Plotter"
    init_inputs = [
            NodeInputType(label='Molecule'),
            NodeInputType(label='MO Coefficients')]

    beta = False

    def get_atom_surface_points(self, axes, size):
        get_atom_surface_points(self.input(0).payload.atom_coords(), axes, size)

    def inputs_ready(self):
        val =  all(hasattr(self.input(i), 'payload') for i in range(len(self.inputs)))
        return val

    def get_isosurface(self, orbital, iso_val, bnds, beta=False):
        mol = self.input(0).payload
        mo_coeff = self.input(1).payload
        return get_isosurface(mol, mo_coeff, orbital, iso_val, bnds, beta)

    def update_event(self, inp=-1):
        if not self.inputs_ready():
            return
        if self.have_gui():
            if len(self.input(1).payload.shape) == 3:
                self.gui.add_beta()
                self.beta = True
            else:
                self.gui.remove_beta()
                self.beta = False
            self.gui.update_orbitallist(self.input(0).payload.nao_nr(), self.input(0).payload.nelec)
            self.gui.update_plot()

    def have_gui(self):
        return hasattr(self, 'gui')

class AOPlotNode(Node):
    """Generates an interactive surface"""

    title = "Atomic Orbital Plotter"
    init_inputs = [
            NodeInputType(label="Molecule")
            ]

    def inputs_ready(self):
        return all(self.input(i) is not None for i in range(len(self.inputs)))

    def get_atom_surface_points(self, axes, size):
        get_atom_surface_points(self.input(0).payload.atom_coords(), axes, size)

    def get_isosurface(self, orbital, iso_val, bnds, beta=False):
        mol = self.input(0).payload
        n = mol.nao_nr()
        mo_coeff = np.asarray([np.zeros((n,n)),np.zeros((n,n))])
        mo_coeff[0,orbital,orbital] = 1
        mo_coeff[1,orbital,orbital] = 1
        return get_isosurface(mol, mo_coeff, orbital,iso_val, bnds, beta)

    def update_event(self, inp=-1):
        if not self.inputs_ready():
            return
        if self.have_gui():
            self.gui.update_orbitallist(0,0,self.input(0).payload.ao_labels())
            self.gui.update_plot()

    def have_gui(self):
        return hasattr(self, 'gui')


class PrintNode(Node):
    title = 'Print'
    init_inputs = [NodeInputType()]

    def update_event(self, inp=-1):
        print(self.input(0))

class plotting(Node):
    pass


export_nodes([
    LinePlotNode,
    MOPlotNode,
    PrintNode,
    AOPlotNode
])


@on_gui_load
def load_gui():
    from . import gui

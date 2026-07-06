from ryven.node_env import *
import numpy as np

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

def get_isosurface(mol, mo_coeff, orbital, iso_val, bnds, beta=False, grid_points=50):
    if orbital < 0:
        # TODO: total electron density
        return None
    if mo_coeff.ndim == 3:
        coeff_1d = mo_coeff[1 if beta else 0, :, orbital]
    else:
        coeff_1d = mo_coeff[:, orbital]
    grid, spacing, origin = _eval_mo_grid(mol, coeff_1d, bnds, grid_points=grid_points)
    pos_verts, pos_faces = _safe_marching(grid, iso_val, spacing, origin)
    neg_verts, neg_faces = _safe_marching(grid, -iso_val, spacing, origin)
    return (pos_verts, pos_faces, neg_verts, neg_faces)

class LinePlotNode(Node):
    """<p><b>LinePlot</b> — plots up to five scalar inputs as a line graph
    (input index on the x-axis, value on the y-axis), redrawn whenever any
    input changes.</p>
    <p>Useful for comparing a handful of values side by side &mdash; e.g.
    wire the Energy outputs of several Fock nodes in to compare guesses or
    functionals. Unconnected inputs are simply skipped.</p>
    <p><b>Inputs</b> &mdash; 1&ndash;5: any numeric scalars.</p>"""

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
    """<p><b>Molecular Orbital Plotter</b> — renders a molecular orbital
    &psi;<sub>i</sub> = &Sigma;<sub>&mu;</sub> C<sub>&mu;i</sub>
    &chi;<sub>&mu;</sub> as a 3D isosurface around the ball-and-stick
    molecule, inline in the node.</p>
    <p>Pick an orbital from the list (labeled relative to the frontier:
    Hono, Luno, ...; unrestricted coefficients get separate &alpha; and
    &beta; lists and viewers). Red and blue lobes are the +c and &minus;c
    surfaces of the wavefunction &mdash; the sign pattern that determines
    bonding vs. antibonding character. <b>Isosurface</b> sets c,
    <b>Opacity</b> the lobe transparency, <b>Grid points</b> the resolution
    of the marching-cubes grid (higher = smoother but slower, cubic cost).
    Left-drag rotates, scroll zooms.</p>
    <p><b>Inputs</b> &mdash; Molecule; MO Coefficients. Feed live SCF
    coefficients in and the picture updates every Step.</p>"""

    title = "Molecular Orbital Plotter"
    init_inputs = [
            NodeInputType(label='Molecule'),
            NodeInputType(label='MO Coefficients')]

    beta = False

    def inputs_ready(self):
        return all(hasattr(self.input(i), 'payload') for i in range(len(self.inputs)))

    def get_isosurface(self, orbital, iso_val, bnds, beta=False, grid_points=50):
        mol = self.input(0).payload
        mo_coeff = self.input(1).payload
        return get_isosurface(mol, mo_coeff, orbital, iso_val, bnds, beta,
                              grid_points=grid_points)

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
    """<p><b>Atomic Orbital Plotter</b> — renders one basis function
    &chi;<sub>&mu;</sub> (a contracted Gaussian centered on an atom) as a
    3D isosurface, inline in the node.</p>
    <p>These are the raw ingredients the SCF mixes into molecular orbitals
    &mdash; inspect them to see what "sto-3g on H" actually looks like, or
    why polarization functions matter. The list shows PySCF's AO labels
    (atom, shell, angular part). Controls are the same as the MO plotter:
    iso value, opacity, grid resolution; left-drag rotates, scroll
    zooms.</p>
    <p><b>Input</b> &mdash; Molecule (the basis comes from it).</p>"""

    title = "Atomic Orbital Plotter"
    init_inputs = [
            NodeInputType(label="Molecule")
            ]

    def inputs_ready(self):
        return all(self.input(i) is not None for i in range(len(self.inputs)))

    def get_isosurface(self, orbital, iso_val, bnds, beta=False, grid_points=50):
        # Each "MO" is a single AO basis function: C = I, so column k is e_k
        # and the iso helper picks out AO #k. Beta is meaningless for AOs.
        mol = self.input(0).payload
        n = mol.nao_nr()
        return get_isosurface(mol, np.eye(n), orbital, iso_val, bnds, beta=False,
                              grid_points=grid_points)

    def update_event(self, inp=-1):
        if not self.inputs_ready():
            return
        if self.have_gui():
            self.gui.update_orbitallist(0,0,self.input(0).payload.ao_labels())
            self.gui.update_plot()

    def have_gui(self):
        return hasattr(self, 'gui')


class PrintNode(Node):
    """<p><b>Print</b> — prints whatever arrives on its input to the
    editor's console whenever it updates. A quick probe for inspecting any
    edge of the flow.</p>
    <p><b>Input</b> &mdash; anything.</p>"""

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

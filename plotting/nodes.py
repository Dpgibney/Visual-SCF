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

def get_isosurface(mol,mo_coeff,orbital,iso_val,bnds,beta=False):
    from skimage import measure
    max_ = max(bnds)
    min_ = min(bnds)
    grid_points = 50
    X = np.linspace(min_,max_,grid_points)
    Y = np.linspace(min_,max_,grid_points)
    Z = np.linspace(min_,max_,grid_points)
    spacing = (max_-min_)/grid_points
    spacing = (spacing,spacing,spacing)
    X, Y, Z = np.meshgrid(X,Y,Z)
    coords = np.vstack([X.ravel(),Y.ravel(),Z.ravel()]).T
    ao = mol.eval_gto('GTOval', coords)

    if orbital >= 0:
        if not beta:
            if len(mo_coeff.shape) == 3:
                orbital_coeff_alpha = mo_coeff[0,:, orbital]
            else:
                orbital_coeff_alpha = mo_coeff[:, orbital]

            mo = np.einsum('i,ki->k',orbital_coeff_alpha,ao)
            mo_values_grid = mo.reshape(X.shape)

            pos_verts, pos_faces, _ ,_ = measure.marching_cubes(mo_values_grid,iso_val,spacing=spacing)
            try:
                neg_verts, neg_faces, _ ,_ = measure.marching_cubes(mo_values_grid,-iso_val,spacing=spacing)
                tmp = neg_verts[:,0] + min_
                neg_verts[:,0] = neg_verts[:,1] + min_
                neg_verts[:,1] = tmp# + min_
                neg_verts[:,2] += min_
            except:
                neg_verts, neg_faces = None, None
            tmp = pos_verts[:,0] + min_
            pos_verts[:,0] = pos_verts[:,1] + min_
            pos_verts[:,1] = tmp# + min_
            pos_verts[:,2] += min_

            alpha = (pos_verts, pos_faces, neg_verts, neg_faces)
            return alpha

        else:
            orbital_coeff_beta = mo_coeff[1,:, orbital]
            mo = np.einsum('i,ki->k',orbital_coeff_beta,ao)
            mo_values_grid = mo.reshape(X.shape)
            pos_verts, pos_faces, _ ,_ = measure.marching_cubes(mo_values_grid,iso_val,spacing=spacing)
            try:
                neg_verts, neg_faces, _ ,_ = measure.marching_cubes(mo_values_grid,-iso_val,spacing=spacing)
                tmp = neg_verts[:,0] + min_
                neg_verts[:,0] = neg_verts[:,1] + min_
                neg_verts[:,1] = tmp# + min_
                neg_verts[:,2] += min_
            except:
                neg_verts, neg_faces = None, None
            tmp = pos_verts[:,0] + min_
            pos_verts[:,0] = pos_verts[:,1] + min_
            pos_verts[:,1] = tmp# + min_
            pos_verts[:,2] += min_

            beta = (pos_verts, pos_faces, neg_verts, neg_faces)
            return beta
    else:
        #TODO impliment total electron density for this one
        #Will also need to update the labeling in the orbitallists
        pass

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


export_nodes([
    LinePlotNode,
    MOPlotNode,
    PrintNode,
    AOPlotNode
])


@on_gui_load
def load_gui():
    from . import gui

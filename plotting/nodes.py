from ryven.node_env import *
import numpy as np

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
        theta = np.linspace(0, 2 * np.pi, 20)
        phi = np.linspace(0, np.pi, 10)
        theta, phi = np.meshgrid(theta, phi)
        r = 0.2
        x = r * np.sin(phi) * np.cos(theta)
        y = r * np.sin(phi) * np.sin(theta)
        z = r * np.cos(phi)
        for atom in self.input(0).payload.atom_coords():
            x_ = x + atom[0]
            y_ = y + atom[1]
            z_ = z + atom[2]
            axes.plot_surface(x_, y_, z_, color='red', alpha=0.8)
            

    def inputs_ready(self):
        val =  all(hasattr(self.input(i), 'payload') for i in range(len(self.inputs)))
        print("inputs ready value",val)
        return val

    def get_isosurface(self,orbital,iso_val,bnds,beta=False):
        from skimage import measure
        print("Bounds: ",bnds)
        max_ = max(bnds)
        min_ = min(bnds)
        X = np.linspace(min_,max_,50)
        Y = np.linspace(min_,max_,50)
        Z = np.linspace(min_,max_,50)
        spacing = (max_-min_)/50
        spacing = (spacing,spacing,spacing)
        X, Y, Z = np.meshgrid(X,Y,Z)
        coords = np.vstack([X.ravel(),Y.ravel(),Z.ravel()]).T
        print("Coords: ",coords)
        mol = self.input(0).payload
        ao = mol.eval_gto('GTOval', coords)
        print("ao: ",ao)
        mo_coeff = self.input(1).payload

        if orbital >= 0:
            if not beta:
                if len(mo_coeff.shape) == 3:
                    orbital_coeff_alpha = mo_coeff[0,:, orbital]
                else:
                    orbital_coeff_alpha = mo_coeff[:, orbital]

                print("orbital coefficient size: ",orbital_coeff_alpha.shape)
                print("ao shape: ",ao.shape)

                mo = np.einsum('i,ki->k',orbital_coeff_alpha,ao)
                mo_values_grid = mo.reshape(X.shape)
                print("mo_values_grid: ",mo_values_grid)

                pos_verts, pos_faces, _ ,_ = measure.marching_cubes(mo_values_grid,iso_val,spacing=spacing)
                print("pos_verts: ",pos_verts)
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
            NodeInputType()]

    def inputs_ready(self):
        return all(self.input(i) is not None for i in range(len(self.inputs)))

    def get_atom_surface_points(self, axes, size):
        theta = np.linspace(0, 2 * np.pi, 20)
        phi = np.linspace(0, np.pi, 10)
        theta, phi = np.meshgrid(theta, phi)
        r = 0.2
        x = r * np.sin(phi) * np.cos(theta)
        y = r * np.sin(phi) * np.sin(theta)
        z = r * np.cos(phi)
        for atom in self.input(0).payload.atom_coords():
            x_ = x + atom[0]
            y_ = y + atom[1]
            z_ = z + atom[2]
            axes.plot_surface(x_, y_, z_, color='red', alpha=0.8)

    def get_isosurface(self,orbital,iso_val,bnds,beta=False):
        max_ = max(bnds)
        min_ = min(bnds)
        X = np.linspace(min_,max_,50)
        Y = np.linspace(min_,max_,50)
        Z = np.linspace(min_,max_,50)
        spacing = (max_-min_)/50
        spacing = (spacing,spacing,spacing)
        if orbital >= 0:
            mol = self.input(0).payload
            orbital_coeff = np.zeros(mol.nao_nr())
            orbital_coeff[orbital] = 1
            #X = np.linspace(bnds[0],bnds[1],30)
            #Y = np.linspace(bnds[2],bnds[3],30)
            #Z = np.linspace(bnds[4],bnds[5],30)
            X, Y, Z = np.meshgrid(X,Y,Z)
            coords = np.vstack([X.ravel(),Y.ravel(),Z.ravel()]).T
            from pyscf import dft
            ao = mol.eval_gto('GTOval', coords)
            print("orbital coeff shape: ",orbital_coeff.shape)
            print("AO shape:",ao.shape)
            mo = np.einsum('i,ki->k',orbital_coeff,ao)
            #for i in range(len(mo)):
            #    if mo[i] < 0:
            #        mo[i] = 0
            mo_values_grid = mo.reshape(X.shape)
            from skimage import measure
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

        return (pos_verts, pos_faces, neg_verts, neg_faces)


    def update_event(self, inp=-1):
        if not self.inputs_ready():
            return
        if self.have_gui():
            print("in ao node")
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

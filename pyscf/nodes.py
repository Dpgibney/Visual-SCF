from ryven.node_env import *
from pyscf import gto, scf, dft

class MolData(Data):
    pass

class MolNode(Node):
    """Defines the molecule and associated basis set"""

    title = 'Molecule'
    tags = ['Integrals']
    init_outputs = [NodeOutputType(label='Molecule')]

    def __init__(self,params):
        super().__init__(params)
        self.atom=""
        self.basis=""
        self.enabled = False

    def update_event(self, inp=-1):
        if not self.enabled:
            return

        mol = gto.M(atom=self.atom,basis=self.basis,verbose=5)
        
        self.set_output_val(0,
            MolData(mol)
        )
        print("Inupdate event")
        print(MolData(mol))

    def have_gui(self):
        return hasattr(self, 'gui')

class FockNode(Node):
    """Returns the fock matrix from the given 1-RDM"""

    title = 'Fock'
    tags = ['Fock']
    init_inputs = [
            NodeInputType(label='Molecule'),
            NodeInputType(label='1-RDM')]
    init_outputs = [NodeOutputType(label='Fock Matrix'),
                    NodeOutputType(label='Energy')]

    def __init__(self,params):
        super().__init__(params)
        from pyscf import scf, dft
        #self.mf = scf.RHF(mol)
        self.xc = None

    def update_mf(self,xc=None):
        #xc = self.gui.xcbox.toPlainText()
        print(xc)
        if not self.inputs_ready():
            return

        self.xc = xc

        mol = self.input(0).payload
        dm = self.input(1).payload

        if xc is not None:
            self.mf = dft.RKS(mol)
            self.mf.xc = self.xc
        else:
            self.mf = scf.RHF(mol)

        hcore = self.mf.get_hcore()
        veff = self.mf.get_veff(mol, dm)
        f = hcore + veff

        e_tot = self.mf.energy_tot(dm,hcore,veff)
        
        self.set_output_val(0,
            Data(f)
        )
        self.set_output_val(1,
            Data(e_tot)
        )

    def inputs_ready(self):
        val =  all(hasattr(self.input(i), 'payload') for i in range(len(self.inputs)))
        return val

    def update_event(self, inp=-1):
        if not self.inputs_ready():
            return

        self.update_mf(self.xc)

        #self.set_output_val(0,
        #    Data(dm=self.input(1).payload)
        #)

    def have_gui(self):
        return hasattr(self, 'gui')

class GetMOCoeffNode(Node):
    """Diagonalizes the Fock Matrix and returns the eigenvectors sorted"""

    title = 'Get MO Coefficients'
    tags = ['MO Coefficients']
    init_inputs = [
            NodeInputType(label='molecule'),
            NodeInputType(label='Fock Matrix')
            ]
    init_outputs = [NodeOutputType(label='MO Coefficients')]

    def __init__(self,params):
        super().__init__(params)
    
    def update_event(self,inp=-1):
        if not self.inputs_ready():
            return

        from scipy.linalg import eigh
        #from numpy import argsort
        ovlp = self.input(0).payload.intor('int1e_ovlp')
        Fock = self.input(1).payload
        e, c = eigh(Fock, ovlp)
        idx = e.argsort()
        e = e[idx]
        c = c[idx]
        print("eigenvalues: ",e)
        #idx = argmax(abs(c.real), axis=0)
        #c[:,c[idx,arange(len(e))].real<0] *= -1
        #c_idx 
        
        self.set_output_val(0, 
            Data(c)
        )

    def inputs_ready(self):
        val =  all(hasattr(self.input(i), 'payload') for i in range(len(self.inputs)))
        return val


class Guess1RDMNode(Node):
    """Allows the user to choose an initial guess density matrix"""

    title = 'Guess 1-RDM'
    tags = ['1-RDM']
    init_inputs = [
            NodeInputType(label='Molecule')]
    init_outputs = [NodeOutputType(label='1-RDM')]

    def __init__(self,params):
        super().__init__(params)
        from pyscf import scf

    def update_guess(self,guess='minao'):
        if not self.inputs_ready():
            return

        mol = self.input(0).payload
        mf = scf.RHF(mol)
        rdm1 = mf.get_init_guess(mol, guess)
        print(rdm1)
        self.set_output_val(0, Data(rdm1))

    def update_event(self, inp=-1):
        if not self.inputs_ready():
            return
        self.update_guess()

    def inputs_ready(self):
        return all(hasattr(self.input(i), 'payload') for i in range(len(self.inputs)))

    def have_gui(self):
        return hasattr(self, 'gui')

class Make1RDMNode(Node):
    """Makes the 1-RDM from a set of MO Coefficients and the Molecules number of electrons (alpha,beta)"""

    title = 'Make 1-RDM'
    tags = ['1-RDM']
    init_inputs = [
            NodeInputType(label='Molecules'),
            NodeInputType(label='MO Coefficients')
            ]
    init_outputs = [NodeOutputType(label='1-RDM')]

    def __init__(self, params):
        super().__init__(params)
        
    def update_event(self, inp=-1):
        from numpy import einsum, zeros
        if not self.inputs_ready():
            return

        MO_coeff = self.input(1).payload
        mol = self.input(0).payload
        if len(MO_coeff.shape) == 3:
           alpha_occs = mol.nelec[0]
           alpha_1rdm = zeros((mol.nao_nr(),mol.nao_nr()))
           for i in range(alpha_occs):
               alpha_1rdm += einsum('i,j->ij',MO_coeff[0,:,i],MO_coeff[0,:,i])
           beta_1rdm = zeros((mol.nao_nr(),mol.nao_nr()))
           beta_occs = mol.nelec[1]
           for i in range(beta_occs):
               beta_1rdm += einsum('i,j->ij',MO_coeff[1,:,i],MO_coeff[1,:,i])
           rdm1 = [alpha_1rdm,beta_1rdm]
        else:
           alpha_occs = mol.nelec[0]
           alpha_1rdm = zeros((mol.nao_nr(),mol.nao_nr()))
           for i in range(alpha_occs):
               alpha_1rdm += einsum('i,j->ij',MO_coeff[:,i],MO_coeff[:,i])
           beta_1rdm = zeros((mol.nao_nr(),mol.nao_nr()))
           beta_occs = mol.nelec[1]
           for i in range(beta_occs):
               beta_1rdm += einsum('i,j->ij',MO_coeff[:,i],MO_coeff[:,i])
           rdm1 = (alpha_1rdm+beta_1rdm)
        self.set_output_val(0,Data(rdm1))

    def inputs_ready(self):
        return all(hasattr(self.input(i), 'payload') for i in range(len(self.inputs)))

    def have_gui(self):
        return hasattr(self, 'gui')
        

class RHFNode(Node):
    title = 'RHF'
    init_inputs = [NodeInputType(label="Molecule")]
    init_outputs = [NodeOutputType(label="Mo Coefficients")]

    def update_event(self, inp=-1):
        mf = scf.RHF(self.input(0).payload)
        mf.kernel()
        print("1RDM: ",mf.make_rdm1())
        print("Fock: ",mf.get_fock())

        self.set_output_val(0,
            Data(mf.mo_coeff)
        )
        print(mf.mo_coeff)

class UHFNode(Node):
    title = 'UHF'
    init_inputs = [NodeInputType(label="Molecule")]
    init_outputs = [NodeOutputType(label="Mo Coefficients")]

    def update_event(self, inp=-1):
        mf = scf.UHF(self.input(0).payload)
        dm_alpha, dm_beta = mf.get_init_guess()
        dm_beta[:1,:1] = 0
        dm = (dm_alpha,dm_beta)
        mf.kernel(dm)

        self.set_output_val(0,
            Data(mf.mo_coeff)
        )
        print("MO coeff",mf.mo_coeff)

class RKSNode(Node):
    title = 'RKS-DFT'
    init_inputs = [NodeInputType()]

    def update_event(self, inp=-1):
        mf = dft.RKS(self.input(0).payload)
        mf.kernel()

    def have_gui(self):
        return hasattr(self, 'gui')

class UKSNode(Node):
    title = 'UKS-DFT'
    init_inputs = [NodeInputType()]

    def update_event(self, inp=-1):
        mf = dft.UKS(self.input(0).payload)
        dm_alpha, dm_beta = mf.get_init_guess()
        dm_beta[:2,:2] = 0
        dm = (dm_alpha,dm_beta)
        mf.kernel(dm)

export_nodes([
    MolNode,
    RHFNode,
    UHFNode,
    RKSNode,
    UKSNode,
    FockNode,
    Guess1RDMNode,
    GetMOCoeffNode,
    Make1RDMNode,
])


@on_gui_load
def load_gui():
    from . import gui

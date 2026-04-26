from ryven.node_env import *
from pyscf import gto, scf, dft

class MolData(Data):
    pass

class MolNode(Node):
    """Defines the molecule and associated basis set"""

    title = 'Molecule'
    tags = ['Integrals']
    init_outputs = [NodeOutputType(label='Molecule')]

    def __init__(self, params):
        super().__init__(params)
        self.atom = ""
        self.basis = ""
        self.charge = 0
        self.spin = 0  # 2S = nα - nβ; 0=singlet, 1=doublet, 2=triplet, ...
        # tracks whether the current state has been built at least once.
        # set by build(); cleared on first text/spin/charge change so we
        # know whether rebuilt() should auto-fire on project load.
        self._built = False

    def build(self):
        """Build the molecule and propagate. Called by the Build button or rebuilt()."""
        if not self.atom or not self.basis:
            return
        mol = gto.M(atom=self.atom, basis=self.basis,
                    charge=self.charge, spin=self.spin, verbose=5)
        self._built = True
        self.set_output_val(0, MolData(mol))

    def update_event(self, inp=-1):
        # called when the Build button fires self.update()
        self.build()

    def have_gui(self):
        return hasattr(self, 'gui')

    def get_state(self):
        return {
            'atom': self.atom,
            'basis': self.basis,
            'charge': self.charge,
            'spin': self.spin,
            'built': self._built,
        }

    def set_state(self, data, version):
        self.atom = data.get('atom', '')
        self.basis = data.get('basis', '')
        self.charge = data.get('charge', 0)
        self.spin = data.get('spin', 0)
        # back-compat: old projects used 'enabled' as the auto-rebuild flag.
        self._built = data.get('built', data.get('enabled', False))

    def rebuilt(self):
        # connections are now in place; if the user saved with the molecule
        # already built, rebuild it so the rest of the flow reanimates.
        if self._built:
            self.build()

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
        self.xc = ''  # empty string = HF; any non-empty = DFT XC functional name

    def update_mf(self, xc):
        # called by the widget when the XC text changes
        self.xc = xc
        self.update_event()

    def inputs_ready(self):
        return all(hasattr(self.input(i), 'payload') for i in range(len(self.inputs)))

    def update_event(self, inp=-1):
        if not self.inputs_ready():
            return

        mol = self.input(0).payload
        dm = self.input(1).payload

        if self.xc:
            self.mf = dft.RKS(mol)
            self.mf.xc = self.xc
        else:
            self.mf = scf.RHF(mol)

        hcore = self.mf.get_hcore()
        veff = self.mf.get_veff(mol, dm)
        f = hcore + veff
        e_tot = self.mf.energy_tot(dm, hcore, veff)

        self.set_output_val(0, Data(f))
        self.set_output_val(1, Data(e_tot))

    def have_gui(self):
        return hasattr(self, 'gui')

    def get_state(self):
        return {'xc': self.xc}

    def set_state(self, data, version):
        self.xc = data.get('xc', '')

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
        ovlp = self.input(0).payload.intor('int1e_ovlp')
        Fock = self.input(1).payload
        e, c = eigh(Fock, ovlp)
        idx = e.argsort()
        e = e[idx]
        c = c[:, idx]
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

    def __init__(self, params):
        super().__init__(params)
        self.guess = 'minao'

    def update_guess(self, guess):
        self.guess = guess
        self.update_event()

    def update_event(self, inp=-1):
        if not self.inputs_ready():
            return
        mol = self.input(0).payload
        rdm1 = scf.RHF(mol).get_init_guess(mol, self.guess)
        self.set_output_val(0, Data(rdm1))

    def inputs_ready(self):
        return all(hasattr(self.input(i), 'payload') for i in range(len(self.inputs)))

    def have_gui(self):
        return hasattr(self, 'gui')

    def get_state(self):
        return {'guess': self.guess}

    def set_state(self, data, version):
        self.guess = data.get('guess', 'minao')

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

        self.set_output_val(0,
            Data(mf.mo_coeff)
        )

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


class SCFStepNode(Node):
    """Single-step iterator that closes the SCF loop.

    Wire Guess1RDM → Initial 1-RDM, and Make1RDM → New 1-RDM (the feedback
    edge). The Current 1-RDM output feeds the Fock node. The Reset button
    re-emits the initial guess and zeros the step counter; Step emits the
    most recent New 1-RDM and increments the counter."""

    title = 'SCF Step'
    tags = ['1-RDM']
    init_inputs = [
        NodeInputType(label='Initial 1-RDM'),
        NodeInputType(label='New 1-RDM'),
    ]
    init_outputs = [NodeOutputType(label='Current 1-RDM')]

    def __init__(self, params):
        super().__init__(params)
        self._has_current = False
        self._step_count = 0

    def reset(self):
        if not self._port_ready(0):
            return
        self._step_count = 0
        self._has_current = True
        self.set_output_val(0, Data(self.input(0).payload))
        self._notify_gui()

    def step(self):
        if not self._port_ready(1):
            return
        self._step_count += 1
        self._has_current = True
        self.set_output_val(0, Data(self.input(1).payload))
        self._notify_gui()

    def _port_ready(self, idx):
        port = self.input(idx)
        return port is not None and hasattr(port, 'payload')

    def _notify_gui(self):
        if hasattr(self, 'gui'):
            self.gui.update_step_count(self._step_count)

    def update_event(self, inp=-1):
        # Auto-Reset on the first arrival of an Initial 1-RDM. Subsequent
        # changes on either input are silently absorbed: New 1-RDM arrivals
        # are the feedback edge (would cause infinite recursion if propagated),
        # and re-fires of Initial don't reset accidentally — the user must
        # click Reset to start over.
        if inp == 0 and not self._has_current:
            self.reset()

    def have_gui(self):
        return hasattr(self, 'gui')

    def get_state(self):
        return {'step_count': self._step_count, 'has_current': self._has_current}

    def set_state(self, data, version):
        self._step_count = data.get('step_count', 0)
        self._has_current = data.get('has_current', False)


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
    SCFStepNode,
])


@on_gui_load
def load_gui():
    from . import gui

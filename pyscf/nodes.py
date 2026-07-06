from ryven.node_env import *
from pyscf import gto, scf, dft, mcscf

class MolData(Data):
    pass

class MolNode(Node):
    """<p><b>Molecule</b> — defines the molecular system and its AO basis;
    the starting point of every flow.</p>
    <p>Enter a geometry, one <tt>Symbol x y z</tt> per line in &Aring;ngstrom
    (e.g. <tt>H 0 0 0; H 0 0 0.74</tt>), a basis set name (<tt>sto-3g</tt>,
    <tt>cc-pvdz</tt>, <tt>def2-svp</tt>, ...), the total charge, and the spin
    as 2S = n<sub>&alpha;</sub> &minus; n<sub>&beta;</sub> (0 = singlet,
    1 = doublet, ...). Edits are staged; <b>Build</b> constructs the PySCF
    molecule (integrals, electron counts) and pushes it downstream.</p>
    <p><b>Output</b> &mdash; Molecule: the built <tt>pyscf.gto.Mole</tt>
    object, passed by reference to all consumers.</p>"""

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

class GeomOptNode(Node):
    """<p><b>Optimize Geometry</b> — relaxes the structure to the nearest
    minimum of the potential energy surface: the arrangement of nuclei
    where all forces vanish.</p>
    <p>Uses PySCF's analytic nuclear gradients with the geomeTRIC optimizer
    (quasi-Newton steps in internal coordinates). Leave <b>XC
    Functional</b> blank to optimize on the Hartree&ndash;Fock surface, or
    enter a functional (<tt>pbe</tt>, <tt>b3lyp</tt>, ...) for DFT.
    Optimization runs a full SCF at every step, so it only recomputes on
    <b>Run</b>; the per-step energies appear in the console and the status
    line reports convergence.</p>
    <p>After a run the node shows the optimization trajectory: drag the
    step slider (or press Play) to watch the geometry walk downhill to the
    minimum, with an energy-vs-step curve overlaid on the 3D view.
    Left-drag rotates, scroll zooms.</p>
    <p><b>Input</b> &mdash; Molecule (charge, spin, and basis are kept).<br>
    <b>Outputs</b> &mdash; Optimized Molecule: a new molecule at the
    relaxed geometry, usable anywhere a Molecule is; Energy: the final
    total energy in Hartree.</p>"""

    title = 'Optimize Geometry'
    tags = ['Geometry']
    init_inputs = [NodeInputType(label='Molecule')]
    init_outputs = [NodeOutputType(label='Optimized Molecule'),
                    NodeOutputType(label='Energy')]

    def __init__(self, params):
        super().__init__(params)
        self.xc = ''  # empty string = HF; any non-empty = DFT XC functional name
        self.maxsteps = 100
        self._has_run = False
        # set by rebuilt() when the project was saved post-run but upstream
        # hasn't repropagated yet; consumed by update_event.
        self._rerun_on_ready = False

    def inputs_ready(self):
        return hasattr(self.input(0), 'payload')

    def _status(self, text):
        if hasattr(self, 'gui'):
            self.gui.set_status(text)

    def run(self):
        from pyscf.geomopt.geometric_solver import GeometryOptimizer
        if not self.inputs_ready():
            return

        # Optimize on a silenced copy: every optimizer step runs a full SCF,
        # and with the upstream molecule's verbose=5 the console would drown
        # in SCF logs. geomeTRIC's one-line-per-step output still shows the
        # energy trajectory.
        mol_in = self.input(0).payload
        mol = mol_in.copy()
        mol.verbose = 0

        if self.xc:
            mf = dft.RKS(mol)
            mf.xc = self.xc
        else:
            mf = scf.RHF(mol)

        # (coords in Bohr, energy) per optimizer cycle — the callback locals
        # of PySCFEngine.calc_new expose both.
        trajectory = []
        opt = GeometryOptimizer(mf)
        opt.callback = lambda envs: trajectory.append(
            (envs['coords'].copy(), envs['energy']))
        opt.max_cycle = self.maxsteps
        try:
            opt.kernel()
        except Exception as e:
            self._status(f'optimization failed: {e}')
            return

        self._has_run = True
        mol_opt = opt.mol
        # restore the input molecule's verbosity so calculations derived
        # from the optimized molecule log to the console as usual.
        mol_opt.verbose = mol_in.verbose
        e_final = trajectory[-1][1] if trajectory else float('nan')
        if opt.converged:
            self._status(f'converged in {len(trajectory)} steps, '
                         f'E = {e_final:.6f} Ha')
        else:
            self._status(f'NOT converged after {len(trajectory)} steps, '
                         f'E = {e_final:.6f} Ha')
        if hasattr(self, 'gui'):
            self.gui.show_trajectory(mol_opt, trajectory)
        self.set_output_val(0, MolData(mol_opt))
        self.set_output_val(1, Data(e_final))

    def update_event(self, inp=-1):
        # A changed input geometry invalidates the previous result, but the
        # (expensive) reoptimization is gated behind the Run button.
        # Exception: a pending rerun queued by rebuilt() fires once the
        # inputs come back.
        if self._rerun_on_ready and self.inputs_ready():
            self._rerun_on_ready = False
            self.run()
        else:
            self._status('inputs changed — press Run')

    def have_gui(self):
        return hasattr(self, 'gui')

    def get_state(self):
        return {
            'xc': self.xc,
            'maxsteps': self.maxsteps,
            'has_run': self._has_run,
        }

    def set_state(self, data, version):
        self.xc = data.get('xc', '')
        self.maxsteps = data.get('maxsteps', 100)
        self._has_run = data.get('has_run', False)

    def rebuilt(self):
        # Re-run only if the user had run this node before saving; defer to
        # update_event if upstream hasn't repropagated yet (rebuilt order
        # isn't guaranteed).
        if not self._has_run:
            return
        if self.inputs_ready():
            self.run()
        else:
            self._rerun_on_ready = True


class FockNode(Node):
    """<p><b>Fock</b> — builds the effective one-electron (Fock) matrix
    F[D] = h<sub>core</sub> + J[D] &minus; &frac12;K[D] for the density
    matrix D on its input, and evaluates the total energy of that
    density.</p>
    <p>Leave the <b>XC Functional</b> field blank for Hartree&ndash;Fock;
    enter any PySCF functional name (<tt>pbe</tt>, <tt>b3lyp</tt>,
    <tt>tpss</tt>, ...) to build the Kohn&ndash;Sham matrix
    h + J + V<sub>xc</sub> instead. Select the node to see the full
    equations in the inspector panel.</p>
    <p><b>Inputs</b> &mdash; Molecule; 1-RDM (density matrix D).<br>
    <b>Outputs</b> &mdash; Fock Matrix (AO basis); Energy: total energy
    E[D] including nuclear repulsion &mdash; watch it drop as the SCF
    iterates.</p>"""

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
    """<p><b>Get MO Coefficients</b> — solves the Roothaan&ndash;Hall
    generalized eigenvalue problem FC = SC&epsilon; in the non-orthogonal
    AO basis (via <tt>scipy.linalg.eigh</tt> with the overlap matrix S).</p>
    <p>Each eigenvector (column of C) is a molecular orbital expressed as a
    linear combination of atomic orbitals; columns are sorted by orbital
    energy &epsilon;<sub>i</sub>, lowest first.</p>
    <p><b>Inputs</b> &mdash; Molecule (supplies S); Fock Matrix.<br>
    <b>Output</b> &mdash; MO Coefficients: the n<sub>AO</sub> &times;
    n<sub>MO</sub> matrix C.</p>"""

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
    """<p><b>Guess 1-RDM</b> — produces the initial one-particle density
    matrix D&#8304; that seeds the SCF iteration.</p>
    <p>The SCF equations are nonlinear, so a starting point is required;
    a good guess converges in fewer steps, and when several SCF solutions
    exist the guess decides which one you find. Choose the scheme in the
    dropdown: <b>minao</b> (superposition of minimal-basis atomic densities,
    default), <b>atom</b> (superposition of atomic HF densities),
    <b>huckel</b> (extended H&uuml;ckel), <b>hcore</b>/<b>1e</b>
    (core Hamiltonian only &mdash; ignores electron repulsion, poor for
    molecules), <b>sap</b> (superposition of atomic potentials).</p>
    <p><b>Input</b> &mdash; Molecule.<br>
    <b>Output</b> &mdash; 1-RDM: the guess density in the AO basis. Wire it
    to the <i>Initial 1-RDM</i> input of an SCF Step node.</p>"""

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
    """<p><b>Make 1-RDM</b> — occupies the lowest MOs with the molecule's
    electrons (aufbau principle) and contracts them into a new density
    matrix D<sub>&mu;&nu;</sub> = &Sigma;<sub>i&#8712;occ</sub>
    C<sub>&mu;i</sub>C<sub>&nu;i</sub> per spin.</p>
    <p>With restricted coefficients (2-D C) both spins share the same
    orbitals and one spin-summed D is emitted; with unrestricted
    coefficients (3-D C) separate &alpha; and &beta; matrices are built
    using n<sub>&alpha;</sub> and n<sub>&beta;</sub> from the molecule.</p>
    <p><b>Inputs</b> &mdash; Molecule (electron counts); MO
    Coefficients.<br>
    <b>Output</b> &mdash; 1-RDM. Wire it back to the <i>New 1-RDM</i> input
    of the SCF Step node to close the SCF loop.</p>"""

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
    """<p><b>RHF</b> — runs a complete restricted Hartree&ndash;Fock
    calculation to convergence in one go (PySCF <tt>scf.RHF</tt>).</p>
    <p>Restricted means each spatial orbital holds an &alpha;,&beta;
    electron pair. This is the "do everything" counterpart of the
    step-by-step Guess &rarr; Fock &rarr; Diagonalize &rarr; SCF Step loop:
    use it when you want converged orbitals without watching the
    iteration.</p>
    <p><b>Input</b> &mdash; Molecule.<br>
    <b>Output</b> &mdash; MO Coefficients of the converged solution
    (energies and convergence report appear in the console).</p>"""

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
    """<p><b>UHF</b> — runs an unrestricted Hartree&ndash;Fock calculation
    to convergence (PySCF <tt>scf.UHF</tt>): &alpha; and &beta; electrons
    get independent spatial orbitals.</p>
    <p>The initial &beta; density is deliberately perturbed to break
    &alpha;/&beta; symmetry &mdash; without this, UHF collapses back onto
    the RHF solution and cannot show symmetry breaking (e.g. correct
    H&#8322; dissociation into two neutral atoms).</p>
    <p><b>Input</b> &mdash; Molecule.<br>
    <b>Output</b> &mdash; MO Coefficients as a (2, n<sub>AO</sub>,
    n<sub>MO</sub>) array: index 0 = &alpha; orbitals, 1 = &beta;.
    Downstream viewers show separate &alpha;/&beta; panels for it.</p>"""

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
    """<p><b>RKS-DFT</b> — runs a restricted Kohn&ndash;Sham DFT
    calculation to convergence (PySCF <tt>dft.RKS</tt>, default
    functional).</p>
    <p>Kohn&ndash;Sham theory replaces Hartree&ndash;Fock exchange with an
    exchange&ndash;correlation functional of the density; the orbitals are
    those of a fictitious non-interacting system reproducing the true
    density. Results go to the console; this node currently has no data
    outputs.</p>
    <p><b>Input</b> &mdash; Molecule.</p>"""

    title = 'RKS-DFT'
    init_inputs = [NodeInputType()]

    def update_event(self, inp=-1):
        mf = dft.RKS(self.input(0).payload)
        mf.kernel()

    def have_gui(self):
        return hasattr(self, 'gui')

class UKSNode(Node):
    """<p><b>UKS-DFT</b> — runs an unrestricted Kohn&ndash;Sham DFT
    calculation to convergence (PySCF <tt>dft.UKS</tt>), with independent
    &alpha; and &beta; densities each seeing its own XC potential.</p>
    <p>As in the UHF node, the &beta; initial guess is perturbed so the
    solution is allowed to break spin symmetry. Results go to the console;
    this node currently has no data outputs.</p>
    <p><b>Input</b> &mdash; Molecule.</p>"""

    title = 'UKS-DFT'
    init_inputs = [NodeInputType()]

    def update_event(self, inp=-1):
        mf = dft.UKS(self.input(0).payload)
        dm_alpha, dm_beta = mf.get_init_guess()
        dm_beta[:2,:2] = 0
        dm = (dm_alpha,dm_beta)
        mf.kernel(dm)


class CASSCFNode(Node):
    """<p><b>CASSCF</b> — Complete Active Space SCF: a full CI expansion
    over a chosen set of "active" orbitals and electrons, optimized
    simultaneously with the orbital shapes. The tool of choice where a
    single determinant fails (bond breaking, diradicals, excited
    states).</p>
    <p>Set the active space &mdash; <b>ncas</b> orbitals holding
    <b>n&alpha;</b> + <b>n&beta;</b> electrons, written CAS(n<sub>e</sub>,
    n<sub>cas</sub>) &mdash; and the number of <b>roots</b> (&gt;1 runs
    equal-weight state-averaging over ground and excited states;
    <b>Display</b> picks which root's density is emitted). The 1-RDM can be
    <b>spin-free</b> (&gamma;<sup>&alpha;</sup> + &gamma;<sup>&beta;</sup>,
    fractional occupations visible) or <b>spin-resolved</b>
    ([&gamma;<sup>&alpha;</sup>, &gamma;<sup>&beta;</sup>]). CASSCF is
    expensive, so it only recomputes on <b>Run</b>; the status line shows
    the energy or what went wrong.</p>
    <p><b>Inputs</b> &mdash; Molecule; MO Coefficients (starting orbitals
    from an upstream RHF/UHF node; UHF orbitals or n&alpha; &ne; n&beta;
    switch to the unrestricted solver).<br>
    <b>Output</b> &mdash; 1-RDM of the displayed root.</p>"""

    title = 'CASSCF'
    tags = ['1-RDM']
    init_inputs = [
        NodeInputType(label='Molecule'),
        NodeInputType(label='MO Coefficients'),
    ]
    init_outputs = [NodeOutputType(label='1-RDM')]

    def __init__(self, params):
        super().__init__(params)
        self.ncas = 2
        self.nelecas_a = 1
        self.nelecas_b = 1
        self.nroots = 1
        self.display_root = 1
        self.representation = 'spin-free'
        self._mc = None
        # number of roots the cached solver was actually run with —
        # self.nroots can drift via the spinbox without a re-run, and
        # _emit_rdm must branch on what the solver holds, not the spinbox.
        self._mc_nroots = 1
        self._has_run = False
        # set by rebuilt() when the project was saved post-run but upstream
        # hasn't repropagated yet; consumed by update_event.
        self._rerun_on_ready = False

    def inputs_ready(self):
        return (hasattr(self.input(0), 'payload')
                and hasattr(self.input(1), 'payload'))

    def _status(self, text):
        if hasattr(self, 'gui'):
            self.gui.set_status(text)

    def _validate(self, mol):
        """Return an error string if the active space can't be built for
        this molecule, else None. Cheap checks only — anything subtler is
        left to pyscf and surfaced by run()'s exception handler."""
        if self.nelecas_a > self.ncas or self.nelecas_b > self.ncas:
            return 'nα and nβ must each be ≤ ncas'
        ncore_elec = mol.nelectron - self.nelecas_a - self.nelecas_b
        if ncore_elec < 0:
            return 'active electrons exceed total electrons'
        if ncore_elec % 2:
            return 'nα+nβ must leave an even number of core electrons'
        if ncore_elec // 2 + self.ncas > mol.nao_nr():
            return 'core + active orbitals exceed the basis size'
        return None

    def run(self):
        import numpy as np
        if not self.inputs_ready():
            return
        mol = self.input(0).payload
        mo  = self.input(1).payload

        err = self._validate(mol)
        if err is not None:
            self._status(err)
            return

        # UHF starting orbitals come through as a 3D (2, nao, nmo) array;
        # spin-imbalanced active spaces also need an unrestricted solver.
        # mcscf.CASSCF can't take either — it converts a UHF reference to
        # ROHF and then rejects the 3D array — so both cases go through
        # UCASSCF, with restricted orbitals duplicated into an (α, β) pair.
        unrestricted = (
            (hasattr(mo, 'ndim') and mo.ndim == 3)
            or self.nelecas_a != self.nelecas_b
        )
        nelecas = (self.nelecas_a, self.nelecas_b)
        if unrestricted:
            mc = mcscf.UCASSCF(scf.UHF(mol), self.ncas, nelecas)
            if getattr(mo, 'ndim', 2) == 2:
                mo = np.array([mo, mo])
        else:
            mc = mcscf.CASSCF(scf.RHF(mol), self.ncas, nelecas)
        if self.nroots > 1:
            mc = mcscf.state_average_(
                mc, weights=[1.0 / self.nroots] * self.nroots
            )
        try:
            mc.kernel(mo)
        except Exception as e:
            self._status(f'CASSCF failed: {e}')
            return

        self._mc = mc
        self._mc_nroots = self.nroots
        self._has_run = True
        suffix = '' if mc.converged else '  (not converged)'
        self._status(f'E = {mc.e_tot:.6f} Ha{suffix}')
        self._emit_rdm()

    def _emit_rdm(self):
        if self._mc is None:
            return
        mc = self._mc
        spin_resolved = (self.representation == 'spin-resolved')

        if self._mc_nroots > 1:
            root = min(self.display_root, self._mc_nroots)
            # mc.fcisolver is a StateAverageFCISolver whose make_rdm1 treats
            # `ci` as a list of states. To get one root's RDM, temporarily
            # swap to the base FCI solver and to a single CI vector.
            sa_solver = mc.fcisolver
            sa_ci     = mc.ci
            mc.fcisolver = sa_solver.undo_state_average()
            mc.ci        = sa_ci[root - 1]
            try:
                if spin_resolved:
                    dm_a, dm_b = mc.make_rdm1s()
                else:
                    dm = mc.make_rdm1()
            finally:
                mc.fcisolver = sa_solver
                mc.ci        = sa_ci
        else:
            if spin_resolved:
                dm_a, dm_b = mc.make_rdm1s()
            else:
                dm = mc.make_rdm1()

        if spin_resolved:
            self.set_output_val(0, Data([dm_a, dm_b]))
        else:
            self.set_output_val(0, Data(dm))

    def update_event(self, inp=-1):
        # Upstream data changed, so the cached solver no longer matches the
        # inputs. Drop it rather than re-emit — pushing an RDM computed from
        # the old molecule/orbitals downstream would be silently wrong. The
        # Run button gates the (expensive) recompute. Exception: a pending
        # rerun queued by rebuilt() fires once the inputs come back.
        self._mc = None
        if self._rerun_on_ready and self.inputs_ready():
            self._rerun_on_ready = False
            self.run()
        else:
            self._status('inputs changed — press Run')

    def have_gui(self):
        return hasattr(self, 'gui')

    def get_state(self):
        return {
            'ncas': self.ncas,
            'nelecas_a': self.nelecas_a,
            'nelecas_b': self.nelecas_b,
            'nroots': self.nroots,
            'display_root': self.display_root,
            'representation': self.representation,
            'has_run': self._has_run,
        }

    def set_state(self, data, version):
        self.ncas           = data.get('ncas', 2)
        self.nelecas_a      = data.get('nelecas_a', 1)
        self.nelecas_b      = data.get('nelecas_b', 1)
        self.nroots         = data.get('nroots', 1)
        self.display_root   = data.get('display_root', 1)
        self.representation = data.get('representation', 'spin-free')
        self._has_run       = data.get('has_run', False)

    def rebuilt(self):
        # Re-run only if the user had run this node before saving — a fresh
        # load shouldn't pay for a CASSCF kernel the user never asked for.
        # Upstream may not have repropagated yet (rebuilt order isn't
        # guaranteed); if so, defer to the next update_event.
        if not self._has_run:
            return
        if self.inputs_ready():
            self.run()
        else:
            self._rerun_on_ready = True


class MOCoeffViewerNode(Node):
    """<p><b>MO Coefficient Viewer</b> — shows the coefficient matrix C as
    a labeled, color-coded table: each column is one molecular orbital
    written as a combination of the atomic orbitals down the rows,
    &phi;<sub>i</sub> = &Sigma;<sub>&mu;</sub> C<sub>&mu;i</sub>
    &chi;<sub>&mu;</sub>.</p>
    <p>Rows carry PySCF AO labels (atom, shell, angular part), columns are
    labeled relative to the frontier orbitals (Hono&minus;1, Hono, Luno,
    Luno+1, ...), and the top row shows occupations. Cells are shaded by
    magnitude &mdash; blue positive, red negative &mdash; so the matrix
    reads like a heat map of the AO&rarr;MO transform. Unrestricted
    coefficients get separate &alpha; and &beta; tables.</p>
    <p><b>Inputs</b> &mdash; Molecule (labels, electron counts); MO
    Coefficients.</p>"""

    title = 'MO Coefficient Viewer'
    tags = ['MO Coefficients']
    init_inputs = [
        NodeInputType(label='Molecule'),
        NodeInputType(label='MO Coefficients'),
    ]

    def inputs_ready(self):
        return all(hasattr(self.input(i), 'payload') for i in range(len(self.inputs)))

    def update_event(self, inp=-1):
        if not self.inputs_ready():
            return
        if hasattr(self, 'gui'):
            self.gui.update_view(self.input(0).payload, self.input(1).payload)

    def have_gui(self):
        return hasattr(self, 'gui')


class SCFStepNode(Node):
    """<p><b>SCF Step</b> — the hand-crank of the self-consistent field
    loop: each <b>Step</b> click advances the iteration
    D<sup>(k)</sup> &rarr; D<sup>(k+1)</sup> by one cycle, so you can watch
    energy, orbitals, and density evolve toward self-consistency.</p>
    <p><b>Wiring:</b> Guess 1-RDM &rarr; <i>Initial 1-RDM</i>; Make 1-RDM
    &rarr; <i>New 1-RDM</i> (the feedback edge); <i>Current 1-RDM</i>
    &rarr; Fock. This closes the cycle Fock &rarr; diagonalize &rarr;
    rebuild density &rarr; back to Fock.</p>
    <p><b>Buttons:</b> <b>Reset</b> re-emits the initial guess and zeros
    the step counter (also happens automatically when the guess first
    arrives); <b>Step</b> emits the most recent fed-back density and
    increments the counter. The iteration is converged when stepping no
    longer changes anything downstream.</p>
    <p><b>Output</b> &mdash; Current 1-RDM: the density of iteration k.</p>"""

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
    GeomOptNode,
    RHFNode,
    UHFNode,
    RKSNode,
    UKSNode,
    CASSCFNode,
    FockNode,
    Guess1RDMNode,
    GetMOCoeffNode,
    Make1RDMNode,
    SCFStepNode,
    MOCoeffViewerNode,
])


@on_gui_load
def load_gui():
    from . import gui

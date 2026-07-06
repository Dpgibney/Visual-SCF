"""LaTeX equation panels for node inspectors (GUI-only module).

Renders LaTeX with matplotlib's mathtext engine (no TeX install needed) into
QPixmaps shown in the inspector dock. Each equation sits on a light "card" so
it stays legible in both the dark and light editor themes.

EQUATIONS maps node class names to inspector content rows:
    ('t', '<html text>')  -- a word-wrapped rich-text caption
    ('e', r'<latex>')     -- a rendered equation (mathtext subset of LaTeX:
                             no environments, one line per entry)

NOTE: plotting2/mathinspector.py is a deliberate copy of the renderer half of
this file — Ryven nodes packages are loaded standalone by path, so
cross-package imports between pyscf/ and plotting2/ are not reliable.
"""
from io import BytesIO

from matplotlib import mathtext

from qtpy.QtCore import Qt
from qtpy.QtGui import QPixmap
from qtpy.QtWidgets import QWidget, QLabel, QVBoxLayout, QScrollArea

from ryvencore_qt import NodeInspectorWidget

EQ_COLOR = '#1b1b1b'   # near-black equation ink
CHIP_BG = '#f6f6f2'    # light card behind each equation
EQ_DPI = 230           # rendered at 2x and shown at devicePixelRatio 2

_pixmap_cache = {}


def latex_pixmap(latex, dpi=EQ_DPI, color=EQ_COLOR):
    """Render a mathtext LaTeX string to a QPixmap (cached)."""
    key = (latex, dpi, color)
    pm = _pixmap_cache.get(key)
    if pm is None:
        buf = BytesIO()
        mathtext.math_to_image(f'${latex}$', buf, dpi=dpi, format='png',
                               color=color)
        pm = QPixmap()
        pm.loadFromData(buf.getvalue(), 'PNG')
        pm.setDevicePixelRatio(2.0)
        _pixmap_cache[key] = pm
    return pm


class _EquationPanel(NodeInspectorWidget, QWidget):
    """Inspector widget listing captions and rendered equations.

    Subclasses (built by equation_inspector) set `rows`. Content lives in a
    QScrollArea because the inspector dock is in a user-resizable splitter.
    """

    rows = ()

    def __init__(self, params):
        NodeInspectorWidget.__init__(self, params)
        QWidget.__init__(self)

        content = QWidget()
        lay = QVBoxLayout()
        lay.setContentsMargins(2, 4, 6, 4)
        lay.setSpacing(6)
        for kind, body in self.rows:
            if kind == 't':
                lbl = QLabel(body)
                lbl.setWordWrap(True)
                lbl.setTextFormat(Qt.RichText)
            else:
                lbl = QLabel()
                try:
                    lbl.setPixmap(latex_pixmap(body))
                except Exception as e:  # a bad equation must not kill the inspector
                    lbl.setText(f'<i>equation failed to render: {e}</i>')
                lbl.setStyleSheet(f'background: {CHIP_BG}; border-radius: 4px;')
                lbl.setMargin(7)
            lay.addWidget(lbl)
        lay.addStretch()
        content.setLayout(lay)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.NoFrame)
        scroll.setWidget(content)

        outer = QVBoxLayout()
        outer.setContentsMargins(0, 0, 0, 0)
        outer.addWidget(scroll)
        self.setLayout(outer)


def equation_inspector(node_key):
    """Return a NodeInspectorWidget subclass showing EQUATIONS[node_key]."""
    return type(f'{node_key}EquationInspector', (_EquationPanel,),
                {'rows': tuple(EQUATIONS[node_key])})


# ---------------------------------------------------------------- content

EQUATIONS = {

    'MolNode': [
        ('t', 'This node fixes the system every downstream node works on: the '
              'clamped-nuclei electronic Hamiltonian (Born–Oppenheimer, atomic '
              'units)'),
        ('e', r'\hat{H} = -\frac{1}{2}\sum_i \nabla_i^2'
              r' - \sum_{i,A}\frac{Z_A}{r_{iA}}'
              r' + \sum_{i<j}\frac{1}{r_{ij}}'
              r' + \sum_{A<B}\frac{Z_A Z_B}{R_{AB}}'),
        ('t', 'and the finite AO basis {χ<sub>μ</sub>} in which all matrices are '
              'expressed. The basis is non-orthogonal, with overlap'),
        ('e', r'S_{\mu\nu} = \int \chi_\mu(\mathbf{r})\,'
              r'\chi_\nu(\mathbf{r})\,d\mathbf{r}'),
    ],

    'GeomOptNode': [
        ('t', 'Finds the nuclear arrangement minimizing the electronic '
              'energy — the nearest local minimum of the potential energy '
              'surface:'),
        ('e', r'\mathbf{R}^{\star} = \mathrm{arg\,min}_{\mathbf{R}}\;'
              r'E(\mathbf{R})'),
        ('t', 'At the minimum the force on every nucleus vanishes:'),
        ('e', r'\mathbf{F}_A = -\nabla_{R_A} E(\mathbf{R}) = 0'),
        ('t', 'PySCF supplies analytic gradients; geomeTRIC takes '
              'quasi-Newton steps in internal coordinates using a local '
              'quadratic model of the surface:'),
        ('e', r'E(\mathbf{R}+\Delta\mathbf{R}) \approx E(\mathbf{R})'
              r' + \mathbf{g}^{T}\Delta\mathbf{R}'
              r' + \frac{1}{2}\,\Delta\mathbf{R}^{T}\mathbf{H}\,'
              r'\Delta\mathbf{R}'),
        ('t', 'Each step is a full SCF at the trial geometry — watch the '
              'energy trajectory in the console.'),
    ],

    'Guess1RDMNode': [
        ('t', 'The SCF equations are nonlinear, so iteration needs a starting '
              'density matrix D⁽⁰⁾:'),
        ('e', r'D^{(0)} \;\longrightarrow\; D^{(k+1)} = '
              r'D\!\left[F\!\left[D^{(k)}\right]\right]'),
        ('t', 'A better guess means fewer SCF steps — and for systems with '
              'several SCF solutions, the guess decides <i>which</i> one you '
              'converge to. Options: <b>minao</b> — superposition of '
              'minimal-basis atomic densities (default); <b>atom</b> — '
              'superposition of atomic HF densities; <b>huckel</b> — extended '
              'Hückel; <b>hcore</b>/<b>1e</b> — diagonalize the core '
              'Hamiltonian, ignoring electron repulsion (poor for molecules); '
              '<b>sap</b> — superposition of atomic potentials.'),
    ],

    'FockNode': [
        ('t', 'Effective one-electron operator of the current density. With '
              'the XC field blank (Hartree–Fock):'),
        ('e', r'F_{\mu\nu}[D] = h_{\mu\nu} + J_{\mu\nu}[D]'
              r' - \frac{1}{2}K_{\mu\nu}[D]'),
        ('e', r'J_{\mu\nu} = \sum_{\lambda\sigma} D_{\lambda\sigma}'
              r'(\mu\nu|\lambda\sigma), \qquad '
              r'K_{\mu\nu} = \sum_{\lambda\sigma} D_{\lambda\sigma}'
              r'(\mu\lambda|\nu\sigma)'),
        ('t', 'With an XC functional set (Kohn–Sham DFT), exact exchange is '
              'replaced by (for hybrids: mixed with) the XC potential:'),
        ('e', r'F_{\mu\nu}[D] = h_{\mu\nu} + J_{\mu\nu}[D]'
              r' + V^{xc}_{\mu\nu}[\rho]'),
        ('t', 'The Energy output is the total energy of the <i>input</i> '
              'density (HF: E<sub>x</sub> = −¼ Tr[DK]; DFT: E<sub>xc</sub> '
              'from the functional):'),
        ('e', r'E[D] = \mathrm{Tr}[D\,h] + E_J[D] + E_{x/xc}[D] + E_{nn}'),
    ],

    'GetMOCoeffNode': [
        ('t', 'Diagonalizes the Fock matrix in the non-orthogonal AO basis — '
              'the Roothaan–Hall generalized eigenvalue problem:'),
        ('e', r'F\,C = S\,C\,\varepsilon'),
        ('t', 'Column i of C holds the AO coefficients of MO φ<sub>i</sub>, '
              'sorted by orbital energy ε<sub>i</sub>. The MOs are orthonormal '
              'with respect to the overlap metric:'),
        ('e', r'C^{T} S\, C = \mathbf{1}, \qquad '
              r'\varphi_i(\mathbf{r}) = \sum_{\mu} C_{\mu i}\,'
              r'\chi_\mu(\mathbf{r})'),
    ],

    'Make1RDMNode': [
        ('t', 'Occupies the lowest orbitals (aufbau) and contracts them into '
              'a one-particle reduced density matrix per spin:'),
        ('e', r'D^{\sigma}_{\mu\nu} = \sum_{i=1}^{N_\sigma} '
              r'C^{\sigma}_{\mu i}\, C^{\sigma}_{\nu i}, '
              r'\qquad \sigma \in \{\alpha,\beta\}'),
        ('t', 'With restricted coefficients (2-D C) both spins share the same '
              'orbitals and the spin-summed matrix D = D<sup>α</sup> + '
              'D<sup>β</sup> is emitted; with unrestricted coefficients (3-D '
              'C) the pair [D<sup>α</sup>, D<sup>β</sup>] is emitted.'),
    ],

    'SCFStepNode': [
        ('t', 'The wired graph computes one sweep of the self-consistent '
              'field map — each Step button click applies it once:'),
        ('e', r'D^{(k+1)} = D\left[\,C\left[\,F[D^{(k)}]\,\right]\right]'),
        ('t', 'Reset re-emits the initial guess (k = 0). The iteration has '
              'converged when the density stops changing — the fixed point'),
        ('e', r'D^{\star} = D\!\left[F[D^{\star}]\right]'),
        ('t', 'at which point orbitals, density, and Fock operator are '
              'mutually consistent ("self-consistent field").'),
    ],

    'RHFNode': [
        ('t', 'Restricted Hartree–Fock: every spatial orbital holds an α,β '
              'pair. PySCF iterates the Roothaan–Hall equations to '
              'convergence:'),
        ('e', r'F[D]\,C = S\,C\,\varepsilon '
              r'\;\;\Longrightarrow\;\; D^{\star}'),
        ('e', r'E_{RHF} = 2\sum_{i}^{N/2} h_{ii} '
              r'+ \sum_{ij}^{N/2}\left(2J_{ij} - K_{ij}\right) + E_{nn}'),
    ],

    'UHFNode': [
        ('t', 'Unrestricted Hartree–Fock: α and β electrons get independent '
              'spatial orbitals — the coupled Pople–Nesbet equations:'),
        ('e', r'F^{\alpha}[D^{\alpha}\!,D^{\beta}]\;C^{\alpha} = '
              r'S\,C^{\alpha}\varepsilon^{\alpha}, \qquad '
              r'F^{\beta}[D^{\alpha}\!,D^{\beta}]\;C^{\beta} = '
              r'S\,C^{\beta}\varepsilon^{\beta}'),
        ('t', 'This node deliberately perturbs the β initial guess to break '
              'α/β symmetry; otherwise UHF just reproduces the RHF solution '
              'and cannot show effects like correct H₂ dissociation.'),
    ],

    'RKSNode': [
        ('t', 'Restricted Kohn–Sham DFT: non-interacting orbitals in an '
              'effective local potential reproduce the interacting density:'),
        ('e', r'\left[-\frac{1}{2}\nabla^2 + v_{ext}(\mathbf{r})'
              r' + v_H[\rho](\mathbf{r}) + v_{xc}[\rho](\mathbf{r})\right]'
              r'\varphi_i = \varepsilon_i\,\varphi_i'),
        ('e', r'\rho(\mathbf{r}) = 2\sum_{i}^{N/2} '
              r'\left|\varphi_i(\mathbf{r})\right|^2'),
    ],

    'UKSNode': [
        ('t', 'Unrestricted Kohn–Sham: separate α and β densities, each '
              'seeing its own XC potential:'),
        ('e', r'v^{\sigma}_{xc}(\mathbf{r}) = '
              r'\frac{\delta E_{xc}[\rho^{\alpha},\rho^{\beta}]}'
              r'{\delta\rho^{\sigma}(\mathbf{r})}'),
        ('t', 'As in the UHF node, the β guess is perturbed so the solution '
              'may break spin symmetry.'),
    ],

    'CASSCFNode': [
        ('t', 'Complete Active Space SCF: full CI over the chosen active '
              'orbitals/electrons, while <i>simultaneously</i> relaxing the '
              'orbitals:'),
        ('e', r'|\Psi\rangle = \sum_{I\,\in\,CAS(n_e,\,n_{cas})} '
              r'c_I\,|\Phi_I\rangle'),
        ('e', r'E_{CASSCF} = \min_{C,\;c}\;'
              r'\langle\Psi|\hat{H}|\Psi\rangle'),
        ('t', 'With several roots, a state-averaged functional is optimized '
              '(equal weights here), giving one orbital set describing all '
              'states:'),
        ('e', r'E_{SA} = \sum_{I=1}^{n_{roots}} w_I\,E_I, '
              r'\qquad w_I = 1/n_{roots}'),
        ('t', 'The emitted 1-RDM is spin-free or spin-resolved:'),
        ('e', r'\gamma_{\mu\nu} = \gamma^{\alpha}_{\mu\nu}'
              r' + \gamma^{\beta}_{\mu\nu}'),
    ],

    'MOCoeffViewerNode': [
        ('t', 'C is the transformation from the AO basis to the MO basis — '
              'each column lists how one MO is built from AOs:'),
        ('e', r'\varphi_i(\mathbf{r}) = \sum_{\mu} C_{\mu i}\,'
              r'\chi_\mu(\mathbf{r})'),
        ('t', 'Rows are AOs (χ<sub>μ</sub>), columns are MOs labeled relative '
              'to the highest/lowest occupied natural orbital (Hono/Luno). '
              'Cell shading is proportional to |C<sub>μi</sub>| — blue '
              'positive, red negative — so the table reads like a heat map of '
              'the transform.'),
    ],
}

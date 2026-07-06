"""LaTeX equation panels for node inspectors (GUI-only module).

Renderer half is a deliberate copy of pyscf/mathinspector.py — Ryven nodes
packages are loaded standalone by path, so cross-package imports between
pyscf/ and plotting2/ are not reliable. Keep the two renderers in sync;
only the EQUATIONS content differs.
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

    'MOPlotNode': [
        ('t', 'Evaluates the selected molecular orbital on a uniform grid '
              'over the molecule:'),
        ('e', r'\psi_i(\mathbf{r}) = \sum_{\mu} C_{\mu i}\,'
              r'\chi_\mu(\mathbf{r})'),
        ('t', 'and draws the two surfaces of constant amplitude (marching '
              'cubes at the chosen iso value c):'),
        ('e', r'\{\mathbf{r} : \psi_i(\mathbf{r}) = +c\}'
              r'\quad \mathrm{and} \quad'
              r'\{\mathbf{r} : \psi_i(\mathbf{r}) = -c\}'),
        ('t', 'Red/blue lobes are the opposite signs of the wavefunction — '
              'the phase pattern that decides bonding vs. antibonding '
              'overlap. |ψ<sub>i</sub>|² is the probability density of an '
              'electron in that orbital. Higher grid resolution refines the '
              'surface at cubic cost in evaluation points.'),
    ],

    'AOPlotNode': [
        ('t', 'Displays one atomic orbital: a contracted Gaussian basis '
              'function (Cartesian form)'),
        ('e', r'\chi_\mu(\mathbf{r}) = \sum_{p} d_{\mu p}\, N_p\;'
              r'x^{l} y^{m} z^{n}\, e^{-\alpha_p r^2}'),
        ('t', 'These are the building blocks the SCF combines into molecular '
              'orbitals. Internally this node is the MO plotter evaluated '
              'with C = 1, so "orbital" k is just basis function '
              'χ<sub>k</sub>:'),
        ('e', r'\psi_k(\mathbf{r}) = \sum_\mu \delta_{\mu k}\,'
              r'\chi_\mu(\mathbf{r}) = \chi_k(\mathbf{r})'),
    ],
}

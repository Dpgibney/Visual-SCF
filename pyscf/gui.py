from qtpy.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QTextEdit,
                            QLabel, QComboBox, QSpinBox, QPushButton,
                            QTableWidget, QTableWidgetItem, QHeaderView,
                            QAbstractItemView, QLineEdit, QSlider)
from qtpy.QtGui import QColor, QPainter, QPen
from qtpy.QtCore import Qt, QTimer, QRectF, QPointF

from ryven.gui_env import *

from . import nodes
from .mathinspector import equation_inspector

# Reusing plotting2's QPainter 3D renderer. This absolute import is safe
# because the Visual-SCF root is on sys.path (Ryven's package loader puts it
# there) and 'plotting2' collides with nothing. The reverse direction — other
# packages importing the local 'pyscf' package — is NOT safe (it shadows the
# PySCF library and only resolves via pkgutil path extension), which is why
# mathinspector.py is duplicated instead.
from plotting2.qpainter3d import Painter3DCanvas

class MolInputWidget(NodeMainWidget, QWidget):
    """Atom geometry, basis set, charge, and spin (2S) for a Molecule node.
    Changes are staged on the node; clicking Build commits and propagates."""

    def __init__(self, params):
        NodeMainWidget.__init__(self, params)
        QWidget.__init__(self)

        layout = QVBoxLayout()

        layout.addWidget(QLabel('Atoms'))
        self.atomsText = QTextEdit()
        self.atomsText.setFixedSize(220, 150)
        self.atomsText.setPlaceholderText("H 0 0 0;\nH 0 0 0.74")
        layout.addWidget(self.atomsText)

        layout.addWidget(QLabel('Basis'))
        self.basisText = QTextEdit()
        self.basisText.setFixedSize(220, 40)
        self.basisText.setPlaceholderText("sto-3g, cc-pvdz, def2-svp, ...")
        layout.addWidget(self.basisText)

        charge_row = QHBoxLayout()
        charge_row.addWidget(QLabel('Charge'))
        self.chargeBox = QSpinBox()
        self.chargeBox.setRange(-10, 10)
        charge_row.addWidget(self.chargeBox)
        charge_row.addStretch()
        layout.addLayout(charge_row)

        spin_row = QHBoxLayout()
        spin_label = QLabel('Spin (2S)')
        spin_label.setToolTip("2S = nα − nβ. 0 = singlet, 1 = doublet, 2 = triplet, ...")
        spin_row.addWidget(spin_label)
        self.spinBox = QSpinBox()
        self.spinBox.setRange(0, 20)
        self.spinBox.setToolTip("2S = nα − nβ. 0 = singlet, 1 = doublet, 2 = triplet, ...")
        spin_row.addWidget(self.spinBox)
        spin_row.addStretch()
        layout.addLayout(spin_row)

        self.buildButton = QPushButton('Build')
        layout.addWidget(self.buildButton)

        # populate from node state before wiring signals so loaded projects
        # don't see spurious change events during initial population.
        self.atomsText.setPlainText(self.node.atom)
        self.basisText.setPlainText(self.node.basis)
        self.chargeBox.setValue(self.node.charge)
        self.spinBox.setValue(self.node.spin)

        self.atomsText.textChanged.connect(self._atom_changed)
        self.basisText.textChanged.connect(self._basis_changed)
        self.chargeBox.valueChanged.connect(self._charge_changed)
        self.spinBox.valueChanged.connect(self._spin_changed)
        self.buildButton.clicked.connect(self._build_clicked)

        self.setLayout(layout)

    def _atom_changed(self):
        self.node.atom = self.atomsText.toPlainText()

    def _basis_changed(self):
        self.node.basis = self.basisText.toPlainText()

    def _charge_changed(self, val):
        self.node.charge = val

    def _spin_changed(self, val):
        self.node.spin = val

    def _build_clicked(self):
        self.update_node()

    
class FockWidget(NodeMainWidget,QWidget):
    def __init__(self, params):
        NodeMainWidget.__init__(self, params)
        QWidget.__init__(self)
        self.resize(50,10)

        layout = QHBoxLayout()

        self.xclabel = QLabel('XC Functional')
        self.xcbox = QTextEdit()
        self.xcbox.setPlaceholderText("Blank for HF")
        self.xcbox.setPlainText(self.node.xc)
        self.xcbox.textChanged.connect(self.update_mf)
        layout.addWidget(self.xclabel)
        layout.addWidget(self.xcbox)

        self.setLayout(layout)

    def update_mf(self):
        self.node.update_mf(self.xcbox.toPlainText())

class GuessWidget(NodeMainWidget,QComboBox):
    def __init__(self,params):
        NodeMainWidget.__init__(self, params)
        QComboBox.__init__(self)

        self.addItems(['minao','atom','huckel','hcore','1e','sap'])
        idx = self.findText(self.node.guess)
        if idx >= 0:
            self.setCurrentIndex(idx)
        self.currentTextChanged.connect(self.update_guess)

    def update_guess(self, guess):
        self.node.update_guess(guess)


class TrajectoryCanvas(Painter3DCanvas):
    """Painter3DCanvas plus a corner overlay charting energy vs. step.

    The overlay is drawn after the base class finishes its paint pass
    (its QPainter is ended by then, so opening a second one is fine)."""

    OVERLAY_LINE = QColor('#1f77b4')
    OVERLAY_MARK = QColor('#d62728')

    def __init__(self, parent=None):
        super().__init__(parent)
        self._energies = []
        self._current_step = -1

    def set_energy_profile(self, energies, current_step):
        self._energies = list(energies)
        self._current_step = current_step
        self.update()

    def paintEvent(self, event):
        super().paintEvent(event)
        if len(self._energies) >= 2:
            self._draw_energy_overlay()

    def _draw_energy_overlay(self):
        es = self._energies
        n = len(es)
        w, h = self.width(), self.height()
        ow = min(int(w * 0.45), 230)
        oh = min(int(h * 0.34), 140)
        if ow < 70 or oh < 45:
            return
        margin = 8
        card = QRectF(w - ow - margin, margin, ow, oh)

        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing, True)

        painter.setPen(QPen(QColor(150, 150, 150), 1))
        painter.setBrush(QColor(252, 252, 250, 235))
        painter.drawRoundedRect(card, 4, 4)

        # plot area inside the card, leaving room for the title line
        px = card.x() + 8
        py = card.y() + 20
        pw = ow - 16
        ph = oh - 32
        emin, emax = min(es), max(es)
        span = max(emax - emin, 1e-12)

        def pt(k):
            x = px + pw * (k / (n - 1))
            y = py + ph * (1.0 - (es[k] - emin) / span)
            return QPointF(x, y)

        painter.setPen(QPen(self.OVERLAY_LINE, 1.6))
        for k in range(n - 1):
            painter.drawLine(pt(k), pt(k + 1))
        painter.setPen(Qt.NoPen)
        painter.setBrush(self.OVERLAY_LINE)
        for k in range(n):
            painter.drawEllipse(pt(k), 2.2, 2.2)
        if 0 <= self._current_step < n:
            painter.setBrush(self.OVERLAY_MARK)
            painter.drawEllipse(pt(self._current_step), 3.8, 3.8)

        painter.setPen(QColor(45, 45, 45))
        font = painter.font()
        font.setPointSizeF(7.5)
        painter.setFont(font)
        painter.drawText(QRectF(card.x() + 8, card.y() + 4, ow - 16, 13),
                         Qt.AlignLeft | Qt.AlignVCenter, 'E (Ha) vs. step')
        # bottom-left: a descending curve ends bottom-right, where the
        # current-step marker would collide with right-aligned text.
        painter.drawText(QRectF(card.x() + 8, card.y() + oh - 14, ow - 16, 12),
                         Qt.AlignLeft | Qt.AlignVCenter,
                         f'ΔE = {es[-1] - es[0]:+.4f}')
        painter.end()


class GeomOptWidget(NodeMainWidget, QWidget):
    """XC functional, step limit, and Run button for GeomOptNode, plus a
    trajectory viewer (3D canvas with energy overlay, step slider, Play)
    that appears after the first run. The optimization is expensive (full
    SCF per step), so it only fires on Run; convergence and the final
    energy land in the status label."""

    PLAY_INTERVAL_MS = 350

    def __init__(self, params):
        NodeMainWidget.__init__(self, params)
        QWidget.__init__(self)

        layout = QVBoxLayout()

        xc_row = QHBoxLayout()
        xc_row.addWidget(QLabel('XC Functional'))
        self.xcEdit = QLineEdit()
        self.xcEdit.setPlaceholderText('Blank for HF')
        xc_row.addWidget(self.xcEdit)
        layout.addLayout(xc_row)

        steps_row = QHBoxLayout()
        steps_row.addWidget(QLabel('Max steps'))
        self.stepsBox = QSpinBox()
        self.stepsBox.setRange(1, 500)
        steps_row.addWidget(self.stepsBox)
        steps_row.addStretch()
        layout.addLayout(steps_row)

        self.runButton = QPushButton('Run')
        layout.addWidget(self.runButton)

        self.statusLabel = QLabel('')
        self.statusLabel.setWordWrap(True)
        layout.addWidget(self.statusLabel)

        # --- trajectory viewer (hidden until the first run) --------------
        self.canvas = TrajectoryCanvas(self)
        self.canvas.setMinimumSize(320, 260)
        self.canvas.hide()
        layout.addWidget(self.canvas)

        traj_row = QHBoxLayout()
        self.playButton = QPushButton('▶')  # ▶
        self.playButton.setFixedWidth(34)
        self.trajSlider = QSlider(Qt.Horizontal)
        self.stepLabel = QLabel('')
        traj_row.addWidget(self.playButton)
        traj_row.addWidget(self.trajSlider, stretch=1)
        traj_row.addWidget(self.stepLabel)
        self._traj_widgets = (self.playButton, self.trajSlider, self.stepLabel)
        for tw in self._traj_widgets:
            tw.hide()
        layout.addLayout(traj_row)

        self._steps = []       # one Mole per optimizer cycle
        self._energies = []
        self._play_timer = QTimer(self)
        self._play_timer.setInterval(self.PLAY_INTERVAL_MS)
        self._play_timer.timeout.connect(self._play_tick)

        # populate from node state before wiring signals so loaded projects
        # don't see spurious change events during initial population.
        self.xcEdit.setText(self.node.xc)
        self.stepsBox.setValue(self.node.maxsteps)

        self.xcEdit.textChanged.connect(self._xc_changed)
        self.stepsBox.valueChanged.connect(self._steps_changed)
        self.runButton.clicked.connect(self._run_clicked)
        self.trajSlider.valueChanged.connect(self._set_step)
        self.playButton.clicked.connect(self._play_clicked)

        self.setLayout(layout)

    def _xc_changed(self, text):
        self.node.xc = text

    def _steps_changed(self, v):
        self.node.maxsteps = v

    def _run_clicked(self):
        self.node.run()

    def set_status(self, text):
        self.statusLabel.setText(text)
        # a wrapped status line grows the widget; without an explicit
        # resize + shape update the node item keeps its old geometry and
        # clips the label under the Run button.
        self.adjustSize()
        self.node.gui.update_shape()

    # ------------------------------------------------------ trajectory view

    def show_trajectory(self, mol_final, trajectory):
        """Build one Mole per optimizer cycle and reveal the viewer.
        trajectory = [(coords in Bohr, energy), ...]."""
        self._stop_play()
        base = mol_final.copy()
        base.verbose = 0   # set_geom_ would log every step's geometry
        self._steps = []
        for coords, _e in trajectory:
            m = base.copy()
            m.set_geom_(coords, unit='Bohr')
            self._steps.append(m)
        self._energies = [e for _c, e in trajectory]
        n = len(self._steps)
        if n == 0:
            return

        self.canvas.show()
        for tw in self._traj_widgets:
            tw.show()
        self.trajSlider.blockSignals(True)
        self.trajSlider.setRange(0, n - 1)
        self.trajSlider.setValue(n - 1)
        self.trajSlider.blockSignals(False)
        self._set_step(n - 1)   # land on the relaxed geometry
        self.adjustSize()
        self.node.gui.update_shape()

    def _set_step(self, k):
        if not (0 <= k < len(self._steps)):
            return
        # set_molecule refits the camera to the new geometry; keep the
        # user's view stable while stepping through the trajectory.
        first = self.canvas._mol_id is None
        cam = (self.canvas._yaw, self.canvas._pitch, self.canvas._distance)
        self.canvas.set_molecule(self._steps[k])
        if not first:
            self.canvas._yaw, self.canvas._pitch, self.canvas._distance = cam
        self.canvas.set_energy_profile(self._energies, k)
        self.stepLabel.setText(
            f'{k + 1}/{len(self._steps)}  E = {self._energies[k]:.6f}')

    # ---------------------------------------------------------- play button

    def _play_clicked(self):
        if self._play_timer.isActive():
            self._stop_play()
            return
        # restart from the beginning when already at the final step
        if self.trajSlider.value() >= len(self._steps) - 1:
            self.trajSlider.setValue(0)
        self._play_timer.start()
        self.playButton.setText('⏸')  # ⏸

    def _play_tick(self):
        k = self.trajSlider.value() + 1
        if k >= len(self._steps):
            self._stop_play()
            return
        self.trajSlider.setValue(k)

    def _stop_play(self):
        self._play_timer.stop()
        self.playButton.setText('▶')


class CASSCFWidget(NodeMainWidget, QWidget):
    """Active-space, root, and 1-RDM-representation controls for CASSCFNode.
    Structural params (ncas, nα, nβ, nroots) gate behind a Run button —
    CASSCF is too expensive to fire on every spinbox change. Display root
    and representation re-emit from the cached solver without re-running."""

    def __init__(self, params):
        NodeMainWidget.__init__(self, params)
        QWidget.__init__(self)

        layout = QVBoxLayout()

        as_row = QHBoxLayout()
        as_row.addWidget(QLabel('ncas'))
        self.ncasBox = QSpinBox(); self.ncasBox.setRange(1, 50)
        as_row.addWidget(self.ncasBox)
        as_row.addWidget(QLabel('nα'))
        self.naBox = QSpinBox(); self.naBox.setRange(0, 50)
        as_row.addWidget(self.naBox)
        as_row.addWidget(QLabel('nβ'))
        self.nbBox = QSpinBox(); self.nbBox.setRange(0, 50)
        as_row.addWidget(self.nbBox)
        as_row.addStretch()
        layout.addLayout(as_row)

        roots_row = QHBoxLayout()
        roots_row.addWidget(QLabel('Roots'))
        self.nrootsBox = QSpinBox(); self.nrootsBox.setRange(1, 20)
        roots_row.addWidget(self.nrootsBox)
        roots_row.addWidget(QLabel('Display'))
        self.displayBox = QSpinBox(); self.displayBox.setRange(1, 1)
        roots_row.addWidget(self.displayBox)
        roots_row.addStretch()
        layout.addLayout(roots_row)

        rep_row = QHBoxLayout()
        rep_row.addWidget(QLabel('1-RDM'))
        self.repBox = QComboBox()
        self.repBox.addItems(['spin-free', 'spin-resolved'])
        rep_row.addWidget(self.repBox)
        rep_row.addStretch()
        layout.addLayout(rep_row)

        self.runButton = QPushButton('Run')
        layout.addWidget(self.runButton)

        # validation errors, solver failures, and post-run energy land here
        # (via node._status) instead of dying in the Qt event loop.
        self.statusLabel = QLabel('')
        self.statusLabel.setWordWrap(True)
        layout.addWidget(self.statusLabel)

        # populate from node state before wiring signals so loaded projects
        # don't see spurious change events during initial population.
        self.ncasBox.setValue(self.node.ncas)
        self.naBox.setValue(self.node.nelecas_a)
        self.nbBox.setValue(self.node.nelecas_b)
        self.nrootsBox.setValue(self.node.nroots)
        self.displayBox.setRange(1, max(1, self.node.nroots))
        self.displayBox.setValue(self.node.display_root)
        idx = self.repBox.findText(self.node.representation)
        if idx >= 0:
            self.repBox.setCurrentIndex(idx)

        self.ncasBox.valueChanged.connect(self._ncas_changed)
        self.naBox.valueChanged.connect(self._na_changed)
        self.nbBox.valueChanged.connect(self._nb_changed)
        self.nrootsBox.valueChanged.connect(self._nroots_changed)
        self.displayBox.valueChanged.connect(self._display_changed)
        self.repBox.currentTextChanged.connect(self._rep_changed)
        self.runButton.clicked.connect(self._run_clicked)

        self.setLayout(layout)

    def _ncas_changed(self, v):
        self.node.ncas = v

    def _na_changed(self, v):
        self.node.nelecas_a = v

    def _nb_changed(self, v):
        self.node.nelecas_b = v

    def _nroots_changed(self, v):
        self.node.nroots = v
        self.displayBox.setRange(1, max(1, v))
        if self.node.display_root > v:
            self.node.display_root = v
            self.displayBox.setValue(v)

    def _display_changed(self, v):
        self.node.display_root = v
        self.node._emit_rdm()

    def _rep_changed(self, text):
        self.node.representation = text
        self.node._emit_rdm()

    def _run_clicked(self):
        self.node.run()

    def set_status(self, text):
        self.statusLabel.setText(text)
        # see GeomOptWidget.set_status — keep the node shape in sync with
        # the wrapped status text.
        self.adjustSize()
        self.node.gui.update_shape()


def _mo_label(idx, n_occ):
    """Return Hono/Luno-style label for MO at zero-based column index."""
    k = idx + 1
    if k < n_occ:
        return f"Hono-{n_occ - k}"
    if k == n_occ:
        return "Hono"
    if k == n_occ + 1:
        return "Luno"
    return f"Luno+{k - n_occ - 1}"


class _MOCoeffTable(QTableWidget):
    """One QTableWidget pre-populated with MO coefficients, AO labels as row
    headers, MO labels as column headers, and an occupation row at the top.
    Cell background is shaded blue (positive) or red (negative) by magnitude."""

    POS_BG = QColor(70, 130, 200)
    NEG_BG = QColor(220, 80, 80)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.setSelectionMode(QAbstractItemView.SingleSelection)
        self.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)

    def populate(self, coeff_2d, ao_labels, n_occ, occ_per_mo):
        nao, nmo = coeff_2d.shape
        self.clear()
        self.setRowCount(nao + 1)  # +1 for occupation row
        self.setColumnCount(nmo)
        self.setHorizontalHeaderLabels([_mo_label(j, n_occ) for j in range(nmo)])
        self.setVerticalHeaderLabels(['Occ.'] + list(ao_labels))

        for j in range(nmo):
            occ = occ_per_mo if (j + 1) <= n_occ else 0
            item = QTableWidgetItem(str(occ))
            item.setTextAlignment(Qt.AlignCenter)
            self.setItem(0, j, item)

        max_abs = max(1e-9, float(abs(coeff_2d).max()))
        for i in range(nao):
            for j in range(nmo):
                val = float(coeff_2d[i, j])
                item = QTableWidgetItem(f"{val:+.4f}")
                item.setTextAlignment(Qt.AlignCenter)
                # Alpha proportional to |val| / max(|val|), so the table reads
                # like a heat map of the AO→MO transform.
                shade = int(200 * min(abs(val) / max_abs, 1.0))
                bg = QColor(self.POS_BG if val >= 0 else self.NEG_BG)
                bg.setAlpha(shade)
                item.setBackground(bg)
                self.setItem(i + 1, j, item)

        self.resizeColumnsToContents()


class MOCoeffViewerWidget(NodeMainWidget, QWidget):
    """Two _MOCoeffTable instances (alpha + beta). Beta is hidden for closed
    shell. Header label states the AO→MO framing explicitly."""

    def __init__(self, params):
        NodeMainWidget.__init__(self, params)
        QWidget.__init__(self)

        self.title_label = QLabel('<b>MO Coefficients</b> &nbsp; — &nbsp; C: AO → MO')
        self.title_label.setTextFormat(Qt.RichText)

        self.alpha_label = QLabel('Alpha')
        self.alpha_label.hide()
        self.alpha_table = _MOCoeffTable()

        self.beta_label = QLabel('Beta')
        self.beta_label.hide()
        self.beta_table = _MOCoeffTable()
        self.beta_table.hide()

        layout = QVBoxLayout()
        layout.addWidget(self.title_label)
        layout.addWidget(self.alpha_label)
        layout.addWidget(self.alpha_table)
        layout.addWidget(self.beta_label)
        layout.addWidget(self.beta_table)
        self.setLayout(layout)
        self.setMinimumSize(420, 280)

    def update_view(self, mol, mo_coeff):
        ao_labels = mol.ao_labels()
        n_alpha, n_beta = mol.nelec
        if mo_coeff.ndim == 3:
            self.alpha_label.show()
            self.beta_label.show()
            self.beta_table.show()
            self.alpha_table.populate(mo_coeff[0], ao_labels, n_alpha, occ_per_mo=1)
            self.beta_table.populate(mo_coeff[1], ao_labels, n_beta, occ_per_mo=1)
        else:
            self.alpha_label.hide()
            self.beta_label.hide()
            self.beta_table.hide()
            # closed-shell: each occupied MO holds 2 electrons
            self.alpha_table.populate(mo_coeff, ao_labels, n_alpha, occ_per_mo=2)


class SCFStepWidget(NodeMainWidget, QWidget):
    def __init__(self, params):
        NodeMainWidget.__init__(self, params)
        QWidget.__init__(self)

        layout = QVBoxLayout()
        self.step_label = QLabel(f'Step: {self.node._step_count}')
        layout.addWidget(self.step_label)

        btn_row = QHBoxLayout()
        self.reset_btn = QPushButton('Reset')
        self.step_btn = QPushButton('Step')
        self.reset_btn.clicked.connect(self.node.reset)
        self.step_btn.clicked.connect(self.node.step)
        btn_row.addWidget(self.reset_btn)
        btn_row.addWidget(self.step_btn)
        layout.addLayout(btn_row)

        self.setLayout(layout)

    def update_step_count(self, n):
        self.step_label.setText(f'Step: {n}')


@node_gui(nodes.MolNode)
class MolNodeGui(NodeGUI):
    main_widget_class = MolInputWidget
    inspector_widget_class = equation_inspector('MolNode')
    wrap_inspector_in_default = True

@node_gui(nodes.FockNode)
class FockNodeGui(NodeGUI):
    main_widget_class = FockWidget
    inspector_widget_class = equation_inspector('FockNode')
    wrap_inspector_in_default = True

    def __init__(self, params):
        super().__init__(params)
        self.main_widget_hidden = False

    #def update_mf(self):
    #    self.main_widget().update_mf(self.xcbox.toPlainText())

@node_gui(nodes.Guess1RDMNode)
class GuessNodeGui(NodeGUI):
    main_widget_class = GuessWidget
    inspector_widget_class = equation_inspector('Guess1RDMNode')
    wrap_inspector_in_default = True

@node_gui(nodes.SCFStepNode)
class SCFStepNodeGui(NodeGUI):
    main_widget_class = SCFStepWidget
    main_widget_pos = 'below ports'
    inspector_widget_class = equation_inspector('SCFStepNode')
    wrap_inspector_in_default = True

    def update_step_count(self, n):
        self.main_widget().update_step_count(n)

@node_gui(nodes.MOCoeffViewerNode)
class MOCoeffViewerNodeGui(NodeGUI):
    main_widget_class = MOCoeffViewerWidget
    main_widget_pos = 'below ports'
    color = '#3344ff'
    inspector_widget_class = equation_inspector('MOCoeffViewerNode')
    wrap_inspector_in_default = True

    def update_view(self, mol, mo_coeff):
        self.main_widget().update_view(mol, mo_coeff)

@node_gui(nodes.CASSCFNode)
class CASSCFNodeGui(NodeGUI):
    main_widget_class = CASSCFWidget
    main_widget_pos = 'below ports'
    color = '#cc6633'
    inspector_widget_class = equation_inspector('CASSCFNode')
    wrap_inspector_in_default = True

    def set_status(self, text):
        self.main_widget().set_status(text)


# The following nodes have no main widget; they are bound to a NodeGUI
# subclass purely to attach the rendered-equation inspector panel. All
# other attributes are NodeGUI defaults, so their appearance is unchanged.

@node_gui(nodes.GeomOptNode)
class GeomOptNodeGui(NodeGUI):
    main_widget_class = GeomOptWidget
    main_widget_pos = 'below ports'
    color = '#2e8b57'
    inspector_widget_class = equation_inspector('GeomOptNode')
    wrap_inspector_in_default = True

    def set_status(self, text):
        self.main_widget().set_status(text)

    def show_trajectory(self, mol_final, trajectory):
        self.main_widget().show_trajectory(mol_final, trajectory)


@node_gui(nodes.GetMOCoeffNode)
class GetMOCoeffNodeGui(NodeGUI):
    inspector_widget_class = equation_inspector('GetMOCoeffNode')
    wrap_inspector_in_default = True

@node_gui(nodes.Make1RDMNode)
class Make1RDMNodeGui(NodeGUI):
    inspector_widget_class = equation_inspector('Make1RDMNode')
    wrap_inspector_in_default = True

@node_gui(nodes.RHFNode)
class RHFNodeGui(NodeGUI):
    inspector_widget_class = equation_inspector('RHFNode')
    wrap_inspector_in_default = True

@node_gui(nodes.UHFNode)
class UHFNodeGui(NodeGUI):
    inspector_widget_class = equation_inspector('UHFNode')
    wrap_inspector_in_default = True

@node_gui(nodes.RKSNode)
class RKSNodeGui(NodeGUI):
    inspector_widget_class = equation_inspector('RKSNode')
    wrap_inspector_in_default = True

@node_gui(nodes.UKSNode)
class UKSNodeGui(NodeGUI):
    inspector_widget_class = equation_inspector('UKSNode')
    wrap_inspector_in_default = True

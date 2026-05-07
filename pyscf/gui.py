from qtpy.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QTextEdit,
                            QLabel, QComboBox, QSpinBox, QPushButton,
                            QTableWidget, QTableWidgetItem, QHeaderView,
                            QAbstractItemView)
from qtpy.QtGui import QColor
from qtpy.QtCore import Qt

from ryven.gui_env import *

from . import nodes

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

@node_gui(nodes.FockNode)
class FockNodeGui(NodeGUI):
    main_widget_class = FockWidget

    def __init__(self, params):
        super().__init__(params)
        self.main_widget_hidden = False

    #def update_mf(self):
    #    self.main_widget().update_mf(self.xcbox.toPlainText())

@node_gui(nodes.Guess1RDMNode)
class GuessNodeGui(NodeGUI):
    main_widget_class = GuessWidget

@node_gui(nodes.SCFStepNode)
class SCFStepNodeGui(NodeGUI):
    main_widget_class = SCFStepWidget
    main_widget_pos = 'below ports'

    def update_step_count(self, n):
        self.main_widget().update_step_count(n)

@node_gui(nodes.MOCoeffViewerNode)
class MOCoeffViewerNodeGui(NodeGUI):
    main_widget_class = MOCoeffViewerWidget
    main_widget_pos = 'below ports'
    color = '#3344ff'

    def update_view(self, mol, mo_coeff):
        self.main_widget().update_view(mol, mo_coeff)

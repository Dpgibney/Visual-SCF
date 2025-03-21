from qtpy.QtWidgets import QSlider, QWidget, QVBoxLayout, QTextEdit, QCheckBox, QTextEdit, QLabel, QHBoxLayout, QComboBox
import qtpy.QtWidgets
from qtpy.QtCore import Qt


from ryven.gui_env import *

from . import nodes

class MolInputWidget(NodeMainWidget,QWidget):#,FigureCanvasQTAgg):
    """a standard QtTextEdit widget, which updates the node
    input it is attached to, every time the input"""
    
    def __init__(self, params):
        NodeMainWidget.__init__(self, params)
        QWidget.__init__(self)
        #self.resize(100,200)

        layout = QVBoxLayout()
        #self.setFixedWidth(400)
        #self.setFixedHeight(400)

        self.atomsText = QTextEdit()
        self.atomsText.setFixedWidth(200)
        self.atomsText.setFixedHeight(200)
        self.atomsText.setPlaceholderText("""Atomic coordinates go here\nH 0 0 0;\nH 1 0 0""")
        self.atomsText.textChanged.connect(self.atom_changed)
        self.basisText = QTextEdit()
        self.basisText.setFixedWidth(200)
        self.basisText.setFixedHeight(50)
        self.basisText.setPlaceholderText("""Basis set name goes here\nsto-3g, cc-pvdz, def2-svp...""")
        self.basisText.textChanged.connect(self.basis_changed)

        self.enable = QCheckBox("Enabled?")
        self.enable.toggled.connect(self.enabled_changed)

        layout.addWidget(self.atomsText)
        layout.addWidget(self.basisText)
        layout.addWidget(self.enable)

        self.setLayout(layout)

    def value_changed(self, val):
        # updates the node input this widget is attached to
        self.update_node_input(Data(val))
        self.update_node()
    
    def get_state(self) -> dict:
        # return the state of the widget
        return {'value': self.value()}
    
    def set_state(self, state: dict):
        # set the state of the widget
        self.setValue(state['value'])

    def atom_changed(self):
        self.node.atom = self.atomsText.toPlainText()
        self.update_node()

    def basis_changed(self):
        self.node.basis = self.basisText.toPlainText()
        self.update_node()

    def enabled_changed(self):
        self.node.enabled = self.sender().isChecked()
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
        self.currentTextChanged.connect(self.update_guess)

    def update_guess(self, guess):
        self.node.update_guess(guess)

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

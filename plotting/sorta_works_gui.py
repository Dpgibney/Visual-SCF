from qtpy.QtWidgets import QSlider, QWidget, QVBoxLayout
import qtpy.QtWidgets
from qtpy.QtCore import Qt


from ryven.gui_env import *

from . import nodes

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
matplotlib.use('Qt5Agg')


class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)
        #self.axes = fig.add_subplot(111)
        super().__init__(fig)
        self.setParent(parent)
        self.ax = ax


class RandSliderWidget(NodeInputWidget,QWidget):#,FigureCanvasQTAgg):
    """a standard Qt slider widget, which updates the node
    input it is attached to, every time the slider value changes"""
    
    def __init__(self, params):
        NodeInputWidget.__init__(self, params)
        #QSlider.__init__(self)
        QWidget.__init__(self)
        #QWidget.resize(600, 600)

        self.x = [0,1]
        self.y = [0,2]

        self.sc = MplCanvas(self, width=4, height=4, dpi=100)
        self.sc.ax.plot(self.x, self.y)
        layout = QVBoxLayout()

        #self.figure = plt.Figure(figsize=(255,10))
        #self.canvas = FigureCanvasQTAGG(self.figure)
        #self.setLayout(QVBoxLayout())
        #self.layout().addWidget(self.canvas)
        self.setFixedWidth(400)
        self.setFixedHeight(400)
        
        #self.setOrientation(Qt.Horizontal)
        #self.setMinimumWidth(100)
        #self.setMinimum(0)
        #self.setMaximum(100)
        #self.setValue(50)
        #self.valueChanged.connect(self.update_plot)
        self.setLayout(QVBoxLayout())
        layout.addWidget(self.sc)
        self.setLayout(layout)
        #self.layout().addWidget(self.sc)
        #self.setCentralWidget(self.sc)

    def update_plot(self):
        print("updating graph")
        self.sc.ax.cla()
        self.sc.ax.plot(self.y,self.x)
        self.sc.draw()
    
    def value_changed(self, val):
        # updates the node input this widget is attached to
        self.update_node_input(Data(val))
    
    def get_state(self) -> dict:
        # return the state of the widget
        return {'value': self.value()}
    
    def set_state(self, state: dict):
        # set the state of the widget
        self.setValue(state['value'])
        #self.updatE_plot(self)

    #def update_event(self, inp=-1):
    #    self.update_plot()
    

@node_gui(nodes.RandNode)
class RandNodeGui(NodeGUI):
    color = '#fcba03'
    
    # register the input widget class
    input_widget_classes = { 'slider': RandSliderWidget }
    
    # attach the slider widget to the first node input
    # display it _below_ the input pin
    init_input_widgets = {
        0: {'name': 'slider', 'pos': 'below'}
    }

    def update_plot(self):
       self.input_widget_classes['slider'].update_plot()



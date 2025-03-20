from qtpy.QtWidgets import QSlider, QWidget, QVBoxLayout, QHBoxLayout, QListWidget, QListWidgetItem, QLineEdit, QLabel
import qtpy.QtWidgets
from qtpy.QtCore import Qt

from ryven.gui_env import *

from . import nodes

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
matplotlib.use('Qt5Agg')
import matplotlib.style as mplstyle
mplstyle.use('fast')

import numpy as np

def get_bounding_box(atom_coords, buffer):
    x_min = atom_coords[0,0]
    x_max = atom_coords[0,0]
    y_min = atom_coords[0,1]
    y_max = atom_coords[0,1]
    z_min = atom_coords[0,2]
    z_max = atom_coords[0,2]

    for i in range(1,len(atom_coords[0])-1):
        if atom_coords[i,0] > x_max:
            x_max = atom_coords[i,0]
        if atom_coords[i,0] < x_min:
            x_min = atom_coords[i,0]
        if atom_coords[i,1] > y_max:
            y_max = atom_coords[i,1]
        if atom_coords[i,1] < y_min:
            y_min = atom_coords[i,1]
        if atom_coords[i,2] > z_max:
            z_max = atom_coords[i,2]
        if atom_coords[i,2] < z_min:
            z_min = atom_coords[i,2]

    x_min -= buffer
    x_max += buffer
    y_min -= buffer
    y_max += buffer
    z_min -= buffer
    z_max += buffer

    return (x_min, x_max, y_min, y_max, z_min, z_max)


class MplCanvas(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)
        #self.axes = fig.add_subplot(111)
        super().__init__(fig)
        self.setParent(parent)
        self.ax = ax

class MplCanvas3d(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig, ax = plt.subplots(figsize=(width, height), dpi=dpi, projection='3d')
        super().__init__(fig)
        self.setParent(parent)
        self.ax = ax

class SurfacePlotWidget(NodeMainWidget, QWidget):
    def __init__(self, params):
        NodeMainWidget.__init__(self, params)
        QWidget.__init__(self)

        self.fig_alpha = Figure()
        self.canvas_alpha = FigureCanvas(self.fig_alpha)
        self.axes_alpha = self.fig_alpha.add_subplot(111, projection='3d')

        self.fig_beta = Figure()
        self.canvas_beta = FigureCanvas(self.fig_beta)
        self.axes_beta = self.fig_beta.add_subplot(111, projection='3d')
        self.canvas_beta.hide()

        self.slider = QSlider(Qt.Horizontal)
        self.slider.setMaximum(1000)
        self.slider.setValue(100)
        self.lineedit = QLineEdit()
        self.lineedit.setText('0.1')
        self.slider.valueChanged.connect(self.update_lineedit)
        self.lineedit.textChanged.connect(self.update_slider)
        self.alpha_label = QLabel('',self)
        self.iso_label = QLabel('Isosurface Value: ',self)
        self.beta_label = QLabel('Beta Orbitals',self)
        self.beta_label.hide()
        
        self.orbitallistalpha = QListWidget()
        self.orbitallistbeta = QListWidget()
        self.orbitallistalpha.currentItemChanged.connect(self.update_plot_alpha)
        self.orbitallistbeta.currentItemChanged.connect(self.update_plot_beta)
        self.orbitallistbeta.hide()

        vlayout = QVBoxLayout()
        hlayout_top = QHBoxLayout()
        hlayout_top.addWidget(self.iso_label)
        hlayout_top.addWidget(self.slider)
        hlayout_top.addWidget(self.lineedit)

        self.vlayout_orbs_left = QVBoxLayout()
        self.vlayout_orbs_left.addWidget(self.alpha_label)
        self.vlayout_orbs_left.addWidget(self.orbitallistalpha)
        self.hlayout_left = QHBoxLayout()
        self.hlayout_left.addLayout(self.vlayout_orbs_left)
        self.hlayout_left.addWidget(self.canvas_alpha)

        self.vlayout_orbs_right = QVBoxLayout()
        self.vlayout_orbs_right.addWidget(self.beta_label)
        self.vlayout_orbs_right.addWidget(self.orbitallistbeta)
        self.hlayout_right = QHBoxLayout()
        self.hlayout_right.addWidget(self.canvas_beta)
        self.hlayout_right.addLayout(self.vlayout_orbs_right)

        self.hlayout_bottom = QHBoxLayout()
        self.hlayout_bottom.addLayout(self.hlayout_left)
        self.hlayout_bottom.addLayout(self.hlayout_right)
        vlayout.addLayout(hlayout_top)
        vlayout.addLayout(self.hlayout_bottom)

        self.setLayout(vlayout)

    def update_orbitallist(self, num_orbs, num_electrons, ao_labels=None):
        print("in update_orbitallist")
        self.orbitallistalpha.clear()
        self.orbitallistbeta.clear()
        #If ao_labels is none then it is the MO plotting
        print("ao_labels",ao_labels)
        print(ao_labels is not None)
        if ao_labels is not None:
            for label in ao_labels:
                item = QListWidgetItem(f"{label}")
                self.orbitallistalpha.addItem(item)
        else:
            for i in range(num_orbs):
                i = i+1
                if i < num_electrons[0]:
                    item = QListWidgetItem(f"Hono -{num_electrons[0]-i}")
                    self.orbitallistalpha.addItem(item)
                elif i == num_electrons[0]:
                    item = QListWidgetItem(f"Hono")
                    self.orbitallistalpha.addItem(item)
                elif i == num_electrons[0] + 1:
                    item = QListWidgetItem(f"Luno")
                    self.orbitallistalpha.addItem(item)
                else:
                    item = QListWidgetItem(f"Luno +{i-num_electrons[0]}")
                    self.orbitallistalpha.addItem(item)

            for i in range(num_orbs):
                i = i+1
                if i < num_electrons[1]:
                    item = QListWidgetItem(f"Hono -{num_electrons[1]-i}")
                    self.orbitallistbeta.addItem(item)
                elif i == num_electrons[1]:
                    item = QListWidgetItem(f"Hono")
                    self.orbitallistbeta.addItem(item)
                elif i == num_electrons[1] + 1:
                    item = QListWidgetItem(f"Luno")
                    self.orbitallistbeta.addItem(item)
                else:
                    item = QListWidgetItem(f"Luno +{i-num_electrons[1]}")
                    self.orbitallistbeta.addItem(item)
        self.orbitallistalpha.setCurrentRow(0)
        self.orbitallistbeta.setCurrentRow(0)



    def update_lineedit(self,value):
        self.lineedit.setText(str(value/1000.0))
        self.update_plot_alpha()
        self.update_plot_beta()

    def update_slider(self,text):
        try:
            value = float(text)
            value *= 1000
            if 0 <= value <= 1000:
                self.slider.setValue(value)
            self.update_plot_alpha()
            self.update_plot_beta()
        except ValueError:
            pass

    def add_beta(self):
        if self.orbitallistbeta.visibleRegion().isEmpty():
           print("add beta")
           self.orbitallistbeta.show()
           self.canvas_beta.show()
           self.beta_label.show()

    def remove_beta(self):
        if not self.orbitallistbeta.visibleRegion().isEmpty():
           print("remove beta")
           self.orbitallistbeta.hide()
           self.canvas_beta.hide()
           self.beta_label.hide()

    def update_plot_alpha(self):
        print("Updating 3d plot")
        print(self.orbitallistalpha.currentRow())
        bnds = get_bounding_box(self.node.input(0).payload.atom_coords(),1.5)
        if self.node.inputs_ready():
            alpha = self.node.get_isosurface(self.orbitallistalpha.currentRow(),float(self.lineedit.text()),bnds=bnds)
            if alpha is not None:
                self.axes_alpha.cla()
                self.axes_alpha.plot_trisurf(alpha[0][:, 0], alpha[0][:, 1], alpha[1], alpha[0][:,2],color='green', edgecolor='none',alpha=0.4)
                if alpha[2] is not None:
                    self.axes_alpha.plot_trisurf(alpha[2][:, 0], alpha[2][:, 1], alpha[3], alpha[2][:,2],color='blue', edgecolor='none',alpha=0.4)
                self.axes_alpha.set_xlim((max(bnds),min(bnds)))
                self.axes_alpha.set_ylim((max(bnds),min(bnds)))
                self.axes_alpha.set_zlim((max(bnds),min(bnds)))
                self.node.get_atom_surface_points(self.axes_alpha,1)

                self.canvas_alpha.draw()

    def update_plot_beta(self):
        print("Updating beta orbitals")
        bnds = get_bounding_box(self.node.input(0).payload.atom_coords(),1.5)
        if self.node.inputs_ready():
            beta = self.node.get_isosurface(self.orbitallistbeta.currentRow(),float(self.lineedit.text()),bnds=bnds,beta=True)
            if beta is not None:
                self.axes_beta.cla()
                self.axes_beta.plot_trisurf(beta[0][:, 0], beta[0][:, 1], beta[1], beta[0][:,2],color='green', edgecolor='none',alpha=0.4)
                if beta[2] is not None:
                    self.axes_beta.plot_trisurf(beta[2][:, 0], beta[2][:, 1], beta[3], beta[2][:,2],color='blue', edgecolor='none',alpha=0.4)
                self.node.get_atom_surface_points(self.axes_beta,1)
                self.axes_beta.set_xlim((max(bnds),min(bnds)))
                self.axes_beta.set_ylim((max(bnds),min(bnds)))
                self.axes_beta.set_zlim((max(bnds),min(bnds)))
                self.canvas_beta.draw()

    def update_plot(self):
        self.update_plot_alpha()
        if not self.orbitallistbeta.visibleRegion().isEmpty():
            self.update_plot_beta()



class LinePlotWidget(NodeMainWidget,QWidget):
    """a standard Qt slider widget, which updates the node
    input it is attached to, every time the slider value changes"""
    
    def __init__(self, params):
        NodeMainWidget.__init__(self, params)
        QWidget.__init__(self)

        self.x = [0,1]
        self.y = [0,2]

        self.sc = MplCanvas(self, width=4, height=4, dpi=100)
        self.sc.ax.plot(self.x, self.y)
        layout = QVBoxLayout()

        self.setFixedWidth(400)
        self.setFixedHeight(400)
        
        self.setLayout(QVBoxLayout())
        layout.addWidget(self.sc)
        self.setLayout(layout)

    def update_plot(self):
        print("updating graph")
        x = []
        y = []
        for i in range(len(self.node.inputs)):
            if hasattr(self.node.input(i),'payload'):
                x.append(i)
                y.append(self.node.input(i).payload)

        self.sc.ax.cla()
        self.sc.ax.plot(x,y)
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


@node_gui(nodes.LinePlotNode)
class LinePlotNodeGui(NodeGUI):
    main_widget_class = LinePlotWidget
    main_widget_pos = 'below ports'
    color = '#fcba03'

    def __init__(self, params):
        super().__init__(params)

        self.main_widget_hidden = False

    def initialized(self):
        if 'show preview' not in self.actions:  # this happens when loading a project
            self.actions['hide preview'] = {'method': self.ac_hide_mw}

    def update_plot(self):
        self.main_widget().update_plot()

    def ac_hide_mw(self):
        self.main_widget().hide()
        del self.actions['hide preview']
        self.actions['show preview'] = {'method': self.ac_show_mw}
        self.main_widget_hidden = True
        self.update_shape()

    def ac_show_mw(self):
        self.main_widget().show()
        del self.actions['show preview']
        self.actions['hide preview'] = {'method': self.ac_hide_mw}
        self.main_widget_hidden = False
        self.update_shape()

    def get_state(self):
        return {
            **super().get_state(),
            'main_widget_hidden': self.main_widget_hidden
        }

    def set_state(self, data):
        super().set_state(data)
        if data['main_widget_hidden']:
            self.ac_hide_mw()
        # shown by default
    
@node_gui(nodes.MOPlotNode)
class MOPlotNodeGui(NodeGUI):
    main_widget_class = SurfacePlotWidget
    main_widget_pos = 'below ports'
    color = '#fcba03'

    def __init__(self, params):
        super().__init__(params)

        self.main_widget_hidden = False

    def initialized(self):
        if 'show preview' not in self.actions:  # this happens when loading a project
            self.actions['hide preview'] = {'method': self.ac_hide_mw}
        self.main_widget().alpha_label.setText('Alpha Orbitals')

    def update_plot(self):
        self.main_widget().update_plot()

    def ac_hide_mw(self):
        self.main_widget().hide()
        del self.actions['hide preview']
        self.actions['show preview'] = {'method': self.ac_show_mw}
        self.main_widget_hidden = True
        self.update_shape()

    def ac_show_mw(self):
        self.main_widget().show()
        del self.actions['show preview']
        self.actions['hide preview'] = {'method': self.ac_hide_mw}
        self.main_widget_hidden = False
        self.update_shape()

    def get_state(self):
        return {
            **super().get_state(),
            'main_widget_hidden': self.main_widget_hidden
        }

    def set_state(self, data):
        super().set_state(data)
        if data['main_widget_hidden']:
            self.ac_hide_mw()
        # shown by default

    def add_beta(self):
        self.main_widget().add_beta()

    def remove_beta(self):
        self.main_widget().remove_beta()

    def update_orbitallist(self, num_orbs, num_electrons):
        self.main_widget().update_orbitallist(num_orbs, num_electrons)

@node_gui(nodes.AOPlotNode)
class AOPlotNodeGui(NodeGUI):
    main_widget_class = SurfacePlotWidget
    main_widget_pos = 'below ports'
    color = '#fcba03'

    def __init__(self, params):
        super().__init__(params)

        self.main_widget_hidden = False

    def initialized(self):
        if 'show preview' not in self.actions:  # this happens when loading a project
            self.actions['hide preview'] = {'method': self.ac_hide_mw}
        self.main_widget().alpha_label.setText('Atomic Orbitals')

    def update_plot(self):
        self.main_widget().update_plot()

    def ac_hide_mw(self):
        self.main_widget().hide()
        del self.actions['hide preview']
        self.actions['show preview'] = {'method': self.ac_show_mw}
        self.main_widget_hidden = True
        self.update_shape()

    def ac_show_mw(self):
        self.main_widget().show()
        del self.actions['show preview']
        self.actions['hide preview'] = {'method': self.ac_hide_mw}
        self.main_widget_hidden = False
        self.update_shape()

    def get_state(self):
        return {
            **super().get_state(),
            'main_widget_hidden': self.main_widget_hidden
        }

    def set_state(self, data):
        super().set_state(data)
        if data['main_widget_hidden']:
            self.ac_hide_mw()
        # shown by default

    def update_orbitallist(self, num_orbs, num_electrons, ao_labels):
        self.main_widget().update_orbitallist(num_orbs, num_electrons, ao_labels)

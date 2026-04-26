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
    mins = atom_coords.min(axis=0) - buffer
    maxs = atom_coords.max(axis=0) + buffer
    return (mins[0], maxs[0], mins[1], maxs[1], mins[2], maxs[2])


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
        # set by set_state during project load; consumed once by the next
        # update_orbitallist call (-1 means "no pending value, keep whatever's selected").
        self._pending_alpha_row = -1
        self._pending_beta_row = -1
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

    def get_state(self):
        return {
            'iso': self.lineedit.text(),
            'alpha_row': self.orbitallistalpha.currentRow(),
            'beta_row': self.orbitallistbeta.currentRow(),
            'beta_visible': self.canvas_beta.isVisibleTo(self),
        }

    def set_state(self, data):
        iso = data.get('iso')
        if iso:
            try:
                self.slider.setValue(int(round(float(iso) * 1000)))
                self.lineedit.setText(iso)
            except ValueError:
                pass
        self._pending_alpha_row = data.get('alpha_row', -1)
        self._pending_beta_row = data.get('beta_row', -1)
        if data.get('beta_visible'):
            self.add_beta()

    def update_orbitallist(self, num_orbs, num_electrons, ao_labels=None):
        # Preserve current selection across re-population (SCF iterations,
        # basis changes, etc). A pending row from set_state takes precedence
        # the first time it's available.
        prev_alpha = self._pending_alpha_row if self._pending_alpha_row >= 0 else self.orbitallistalpha.currentRow()
        prev_beta = self._pending_beta_row if self._pending_beta_row >= 0 else self.orbitallistbeta.currentRow()
        self._pending_alpha_row = -1
        self._pending_beta_row = -1
        self.orbitallistalpha.clear()
        self.orbitallistbeta.clear()
        #If ao_labels is none then it is the MO plotting
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
        self._restore_row(self.orbitallistalpha, prev_alpha)
        self._restore_row(self.orbitallistbeta, prev_beta)

    @staticmethod
    def _restore_row(list_widget, row):
        n = list_widget.count()
        if n == 0:
            return
        list_widget.setCurrentRow(row if 0 <= row < n else 0)



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
           self.orbitallistbeta.show()
           self.canvas_beta.show()
           self.beta_label.show()

    def remove_beta(self):
        if not self.orbitallistbeta.visibleRegion().isEmpty():
           self.orbitallistbeta.hide()
           self.canvas_beta.hide()
           self.beta_label.hide()

    def _draw_lobes(self, axes, surface, bnds):
        axes.cla()
        pos_verts, pos_faces, neg_verts, neg_faces = surface
        if pos_verts is not None:
            axes.plot_trisurf(pos_verts[:, 0], pos_verts[:, 1], pos_faces, pos_verts[:, 2],
                              color='green', edgecolor='none', alpha=0.4)
        if neg_verts is not None:
            axes.plot_trisurf(neg_verts[:, 0], neg_verts[:, 1], neg_faces, neg_verts[:, 2],
                              color='blue', edgecolor='none', alpha=0.4)
        axes.set_xlim((max(bnds), min(bnds)))
        axes.set_ylim((max(bnds), min(bnds)))
        axes.set_zlim((max(bnds), min(bnds)))

    def update_plot_alpha(self):
        if not self.node.inputs_ready():
            return
        bnds = get_bounding_box(self.node.input(0).payload.atom_coords(), 1.5)
        alpha = self.node.get_isosurface(self.orbitallistalpha.currentRow(),
                                         float(self.lineedit.text()), bnds=bnds)
        if alpha is None:
            return
        self._draw_lobes(self.axes_alpha, alpha, bnds)
        self.node.get_atom_surface_points(self.axes_alpha, 1)
        self.canvas_alpha.draw()

    def update_plot_beta(self):
        if not self.node.inputs_ready():
            return
        bnds = get_bounding_box(self.node.input(0).payload.atom_coords(), 1.5)
        beta = self.node.get_isosurface(self.orbitallistbeta.currentRow(),
                                        float(self.lineedit.text()), bnds=bnds, beta=True)
        if beta is None:
            return
        self._draw_lobes(self.axes_beta, beta, bnds)
        self.node.get_atom_surface_points(self.axes_beta, 1)
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
        x = []
        y = []
        for i in range(len(self.node.inputs)):
            if hasattr(self.node.input(i),'payload'):
                x.append(i)
                y.append(self.node.input(i).payload)

        self.sc.ax.cla()
        self.sc.ax.plot(x,y)
        self.sc.draw()


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

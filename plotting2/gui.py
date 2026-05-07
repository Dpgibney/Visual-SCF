from qtpy.QtWidgets import (QSlider, QWidget, QVBoxLayout, QHBoxLayout,
                            QListWidget, QListWidgetItem, QLineEdit, QLabel,
                            QSpinBox)
from qtpy.QtCore import Qt

from ryven.gui_env import *

from . import nodes
from .qpainter3d import Painter3DCanvas

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
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
        super().__init__(fig)
        self.setParent(parent)
        self.ax = ax

DEFAULT_POS_COLOR = '#d62728'  # red
DEFAULT_NEG_COLOR = '#1f77b4'  # blue
DEFAULT_OPACITY = 0.5
GRID_BUFFER = 1.5  # Bohr around the molecule for the iso grid bbox
DEFAULT_GRID_POINTS = 50
GRID_POINTS_MIN = 20
GRID_POINTS_MAX = 80
GRID_POINTS_STEP = 5


class SurfacePlotWidget(NodeMainWidget, QWidget):
    """Embedded controls (orbital list, iso/opacity sliders, grid resolution)
    plus inline 3D viewers — one for alpha orbitals, one for beta — rendered
    by Painter3DCanvas (pure QPainter, embeds inside Ryven's QGraphicsScene)."""

    def __init__(self, params):
        NodeMainWidget.__init__(self, params)
        QWidget.__init__(self)

        # --- iso slider + lineedit ---------------------------------------
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setMaximum(1000)
        self.slider.setValue(100)
        self.lineedit = QLineEdit('0.1')
        self.lineedit.setMaximumWidth(60)

        self.opacity_slider = QSlider(Qt.Horizontal)
        self.opacity_slider.setMaximum(100)
        self.opacity_slider.setValue(int(DEFAULT_OPACITY * 100))
        self.opacity_lineedit = QLineEdit(f'{DEFAULT_OPACITY:.2f}')
        self.opacity_lineedit.setMaximumWidth(60)

        self.slider.valueChanged.connect(self._iso_slider_changed)
        self.lineedit.textChanged.connect(self._iso_text_changed)
        self.opacity_slider.valueChanged.connect(self._opacity_slider_changed)
        self.opacity_lineedit.textChanged.connect(self._opacity_text_changed)

        # marching-cubes grid resolution (per-node, persisted). The QPainter
        # renderer is fast at 30³ but visibly chunky; 50³ is the sweet spot
        # for sto-3g / minimal basis MOs. Users dial down for cc-pvdz where
        # face counts climb.
        self.grid_spin = QSpinBox()
        self.grid_spin.setRange(GRID_POINTS_MIN, GRID_POINTS_MAX)
        self.grid_spin.setSingleStep(GRID_POINTS_STEP)
        self.grid_spin.setValue(DEFAULT_GRID_POINTS)
        self.grid_spin.valueChanged.connect(self.update_plot)

        self.alpha_label = QLabel('Alpha Orbitals')
        self.beta_label = QLabel('Beta Orbitals')
        self.beta_label.hide()

        # --- orbital lists -----------------------------------------------
        self.orbitallistalpha = QListWidget()
        self.orbitallistbeta = QListWidget()
        # set by set_state during project load; consumed once by the next
        # update_orbitallist call (-1 means "no pending value").
        self._pending_alpha_row = -1
        self._pending_beta_row = -1
        self._beta_visible = False
        self.orbitallistalpha.currentItemChanged.connect(self.update_plot_alpha)
        self.orbitallistbeta.currentItemChanged.connect(self.update_plot_beta)
        self.orbitallistbeta.hide()

        # --- 3D canvases -------------------------------------------------
        self.alpha_canvas = Painter3DCanvas(self)
        self.beta_canvas = Painter3DCanvas(self)
        self.beta_canvas.hide()
        # share camera so rotating/zooming one canvas mirrors the other
        self.alpha_canvas.link_camera(self.beta_canvas)

        # --- layout ------------------------------------------------------
        iso_row = QHBoxLayout()
        iso_row.addWidget(QLabel('Isosurface'))
        iso_row.addWidget(self.slider)
        iso_row.addWidget(self.lineedit)

        opacity_row = QHBoxLayout()
        opacity_row.addWidget(QLabel('Opacity'))
        opacity_row.addWidget(self.opacity_slider)
        opacity_row.addWidget(self.opacity_lineedit)

        grid_row = QHBoxLayout()
        grid_row.addWidget(QLabel('Grid points'))
        grid_row.addWidget(self.grid_spin)
        grid_row.addStretch()

        left_col = QVBoxLayout()
        left_col.addWidget(self.alpha_label)
        left_col.addWidget(self.orbitallistalpha)

        right_col = QVBoxLayout()
        right_col.addWidget(self.beta_label)
        right_col.addWidget(self.orbitallistbeta)

        lists_row = QHBoxLayout()
        lists_row.addLayout(left_col)
        lists_row.addLayout(right_col)

        canvases_row = QHBoxLayout()
        canvases_row.addWidget(self.alpha_canvas, stretch=1)
        canvases_row.addWidget(self.beta_canvas, stretch=1)

        outer = QVBoxLayout()
        outer.addLayout(iso_row)
        outer.addLayout(opacity_row)
        outer.addLayout(grid_row)
        outer.addLayout(lists_row)
        outer.addLayout(canvases_row)
        self.setLayout(outer)

        self.setMinimumWidth(420)
        self.setMinimumHeight(560)

    # ---------------------------------------------------------------- state

    def get_state(self):
        return {
            'iso': self.lineedit.text(),
            'opacity': self.opacity_lineedit.text(),
            'alpha_row': self.orbitallistalpha.currentRow(),
            'beta_row': self.orbitallistbeta.currentRow(),
            'beta_visible': self._beta_visible,
            'grid_points': self.grid_spin.value(),
        }

    def set_state(self, data):
        iso = data.get('iso')
        if iso:
            try:
                self.slider.setValue(int(round(float(iso) * 1000)))
                self.lineedit.setText(iso)
            except ValueError:
                pass
        opacity = data.get('opacity')
        if opacity:
            try:
                self.opacity_slider.setValue(int(round(float(opacity) * 100)))
                self.opacity_lineedit.setText(opacity)
            except ValueError:
                pass
        gp = data.get('grid_points')
        if isinstance(gp, int):
            self.grid_spin.setValue(max(GRID_POINTS_MIN, min(GRID_POINTS_MAX, gp)))
        self._pending_alpha_row = data.get('alpha_row', -1)
        self._pending_beta_row = data.get('beta_row', -1)
        if data.get('beta_visible'):
            self.add_beta()

    # ---------------------------------------------------------- orbital list

    def update_orbitallist(self, num_orbs, num_electrons, ao_labels=None):
        prev_alpha = self._pending_alpha_row if self._pending_alpha_row >= 0 else self.orbitallistalpha.currentRow()
        prev_beta = self._pending_beta_row if self._pending_beta_row >= 0 else self.orbitallistbeta.currentRow()
        self._pending_alpha_row = -1
        self._pending_beta_row = -1
        self.orbitallistalpha.clear()
        self.orbitallistbeta.clear()
        if ao_labels is not None:
            for label in ao_labels:
                self.orbitallistalpha.addItem(QListWidgetItem(label))
        else:
            self._populate_mo_list(self.orbitallistalpha, num_orbs, num_electrons[0])
            self._populate_mo_list(self.orbitallistbeta, num_orbs, num_electrons[1])
        self._restore_row(self.orbitallistalpha, prev_alpha)
        self._restore_row(self.orbitallistbeta, prev_beta)

    @staticmethod
    def _populate_mo_list(list_widget, num_orbs, n_occ):
        for k in range(1, num_orbs + 1):
            if k < n_occ:
                label = f"Hono -{n_occ - k}"
            elif k == n_occ:
                label = "Hono"
            elif k == n_occ + 1:
                label = "Luno"
            else:
                label = f"Luno +{k - n_occ}"
            list_widget.addItem(QListWidgetItem(label))

    @staticmethod
    def _restore_row(list_widget, row):
        n = list_widget.count()
        if n == 0:
            return
        list_widget.setCurrentRow(row if 0 <= row < n else 0)

    # -------------------------------------------------------- iso/opacity ui

    def _iso_slider_changed(self, value):
        text = f'{value / 1000.0:.3f}'
        if self.lineedit.text() != text:
            self.lineedit.setText(text)
        self.update_plot()

    def _iso_text_changed(self, text):
        try:
            value = float(text)
        except ValueError:
            return
        slider_val = int(round(value * 1000))
        if 0 <= slider_val <= 1000 and self.slider.value() != slider_val:
            self.slider.setValue(slider_val)
        self.update_plot()

    def _opacity_slider_changed(self, value):
        text = f'{value / 100.0:.2f}'
        if self.opacity_lineedit.text() != text:
            self.opacity_lineedit.setText(text)
        self.update_plot()

    def _opacity_text_changed(self, text):
        try:
            value = float(text)
        except ValueError:
            return
        slider_val = int(round(value * 100))
        if 0 <= slider_val <= 100 and self.opacity_slider.value() != slider_val:
            self.opacity_slider.setValue(slider_val)
        self.update_plot()

    def _opacity(self):
        try:
            return max(0.0, min(1.0, float(self.opacity_lineedit.text())))
        except ValueError:
            return DEFAULT_OPACITY

    # ----------------------------------------------------- beta panel toggle

    def add_beta(self):
        if not self._beta_visible:
            self._beta_visible = True
            self.orbitallistbeta.show()
            self.beta_label.show()
            self.beta_canvas.show()

    def remove_beta(self):
        if self._beta_visible:
            self._beta_visible = False
            self.orbitallistbeta.hide()
            self.beta_label.hide()
            self.beta_canvas.hide()
            self.beta_canvas.clear_lobes()

    # ------------------------------------------------------------ rendering

    def _iso_value(self):
        try:
            return float(self.lineedit.text())
        except ValueError:
            return 0.1

    def update_plot_alpha(self):
        if not self.node.inputs_ready():
            return
        mol = self.node.input(0).payload
        self.alpha_canvas.set_molecule(mol)
        bnds = get_bounding_box(mol.atom_coords(), GRID_BUFFER)
        surface = self.node.get_isosurface(self.orbitallistalpha.currentRow(),
                                           self._iso_value(), bnds=bnds,
                                           grid_points=self.grid_spin.value())
        self.alpha_canvas.set_lobes(surface, self._opacity(),
                                    pos_color=DEFAULT_POS_COLOR,
                                    neg_color=DEFAULT_NEG_COLOR)

    def update_plot_beta(self):
        if not self.node.inputs_ready() or not self._beta_visible:
            return
        mol = self.node.input(0).payload
        self.beta_canvas.set_molecule(mol)
        bnds = get_bounding_box(mol.atom_coords(), GRID_BUFFER)
        surface = self.node.get_isosurface(self.orbitallistbeta.currentRow(),
                                           self._iso_value(), bnds=bnds, beta=True,
                                           grid_points=self.grid_spin.value())
        self.beta_canvas.set_lobes(surface, self._opacity(),
                                   pos_color=DEFAULT_POS_COLOR,
                                   neg_color=DEFAULT_NEG_COLOR)

    def update_plot(self):
        self.update_plot_alpha()
        if self._beta_visible:
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

"""Embedded QPainter-based 3D viewer for atoms, bonds, and isosurfaces.

Pure Qt 2D paint (no native GL surface), so it embeds inside Ryven's
QGraphicsScene without crashing. Replaces the prior pop-out PyVista path.

Distances throughout are in Bohr (matching pyscf.gto.M defaults).
"""
from __future__ import annotations

import numpy as np
from qtpy.QtCore import Qt, QPointF
from qtpy.QtGui import (QPainter, QColor, QPen, QBrush, QPolygonF,
                        QRadialGradient)
from qtpy.QtWidgets import QWidget, QSizePolicy

from . import molecule_render as mr


# Lighting (in eye space — light moves with the camera)
_LIGHT_DIR = np.array([0.3, 0.5, 1.0])
_LIGHT_DIR = _LIGHT_DIR / np.linalg.norm(_LIGHT_DIR)
AMBIENT = 0.35
DIFFUSE = 0.65

ROT_SENSITIVITY_DEG_PER_PX = 0.5
ZOOM_FACTOR_PER_NOTCH = 0.9        # per 120 wheel units
MIN_DISTANCE = 1.0                  # Bohr
MAX_DISTANCE = 200.0
NEAR_Z = 0.05                       # eye-space z below this is "behind camera"

DEFAULT_BG = QColor(255, 255, 255)
BOND_COLOR = QColor(120, 120, 120)
ATOM_OUTLINE = QColor(40, 40, 40)


def _hex_to_rgb(h):
    h = h.lstrip('#')
    return np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)], dtype=np.float32)


def _rotation_matrix(yaw, pitch):
    cy, sy = np.cos(yaw), np.sin(yaw)
    cp, sp = np.cos(pitch), np.sin(pitch)
    Ry = np.array([[ cy, 0,  sy],
                   [  0, 1,   0],
                   [-sy, 0,  cy]], dtype=np.float64)
    Rx = np.array([[1,  0,   0],
                   [0, cp, -sp],
                   [0, sp,  cp]], dtype=np.float64)
    return Rx @ Ry


def _lighten(c, frac):
    return QColor(min(255, int(c.red()   + (255 - c.red())   * frac)),
                  min(255, int(c.green() + (255 - c.green()) * frac)),
                  min(255, int(c.blue()  + (255 - c.blue())  * frac)))


def _darken(c, frac):
    return QColor(int(c.red()   * (1 - frac)),
                  int(c.green() * (1 - frac)),
                  int(c.blue()  * (1 - frac)))


class Painter3DCanvas(QWidget):
    """Embedded QPainter-based 3D viewer.

    Painter's-algorithm depth-sorted rendering of:
      - iso-surface triangle meshes (Lambert-shaded, backface-culled),
      - bonds (depth-shaded gray sticks),
      - atoms (radial-gradient discs with CPK colors).

    Mouse left-drag rotates around the scene centroid; wheel zooms.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(220, 220)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.setAutoFillBackground(True)
        # accept focus so we can grab the mouse cleanly from inside Ryven's
        # QGraphicsProxyWidget — without this, QGraphicsScene reads click+drag
        # as "drag the node" and we never see mouseMoveEvent.
        self.setFocusPolicy(Qt.StrongFocus)

        # scene state
        self._mol_id = None
        self._atoms = []            # list of dicts from atom_render_data
        self._bonds = []            # list of (i, j) atom-index pairs
        self._atom_pos = np.zeros((0, 3))
        self._lobes = None          # (pos_v, pos_f, neg_v, neg_f) or None
        self._opacity = 0.5
        self._pos_color = _hex_to_rgb('#d62728')
        self._neg_color = _hex_to_rgb('#1f77b4')

        # camera
        self._yaw = 0.0
        self._pitch = 0.3
        self._distance = 15.0
        self._target = np.zeros(3)
        self._scene_radius = 5.0     # used to fit the view on set_molecule

        # mouse drag
        self._drag_active = False
        self._drag_start = (0, 0)
        self._drag_yaw = 0.0
        self._drag_pitch = 0.0

        # peers whose camera tracks ours (and vice versa via link_camera).
        # Used to share a camera across the alpha/beta canvases.
        self._peers = []

    # ---- public API ---------------------------------------------------

    def set_molecule(self, mol):
        """Refresh atom/bond geometry. Idempotent on the same mol object."""
        if mol is None or id(mol) == self._mol_id:
            return
        self._mol_id = id(mol)
        self._atoms = mr.atom_render_data(mol)
        self._bonds = list(mr.detect_bonds(mol))
        if self._atoms:
            self._atom_pos = np.array([a['pos'] for a in self._atoms])
            self._target = self._atom_pos.mean(axis=0)
            radii = np.linalg.norm(self._atom_pos - self._target, axis=1)
            r = float(radii.max()) if len(radii) else 0.0
            self._scene_radius = max(r + 2.0, 2.0)
            self._distance = max(self._scene_radius * 3.0, 6.0)
        else:
            self._atom_pos = np.zeros((0, 3))
            self._target = np.zeros(3)
            self._scene_radius = 5.0
            self._distance = 15.0
        self.update()

    def set_lobes(self, surface, opacity, pos_color=None, neg_color=None):
        """Set iso-surface meshes. surface = (pos_v, pos_f, neg_v, neg_f)."""
        self._lobes = surface
        self._opacity = float(opacity)
        if pos_color is not None:
            self._pos_color = _hex_to_rgb(pos_color)
        if neg_color is not None:
            self._neg_color = _hex_to_rgb(neg_color)
        self.update()

    def clear_lobes(self):
        self._lobes = None
        self.update()

    def link_camera(self, other):
        """Mirror user-driven camera changes (yaw/pitch/distance) between
        self and `other`. Bidirectional after one call."""
        if other is self or other is None:
            return
        if other not in self._peers:
            self._peers.append(other)
        if self not in other._peers:
            other._peers.append(self)

    def _broadcast_camera(self):
        for p in self._peers:
            p._yaw = self._yaw
            p._pitch = self._pitch
            p._distance = self._distance
            p.update()

    # ---- mouse --------------------------------------------------------

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            self._drag_active = True
            self._drag_start = (event.x(), event.y())
            self._drag_yaw = self._yaw
            self._drag_pitch = self._pitch
            self.setFocus(Qt.MouseFocusReason)
            event.accept()
            return
        event.ignore()

    def mouseMoveEvent(self, event):
        if self._drag_active:
            dx = event.x() - self._drag_start[0]
            dy = event.y() - self._drag_start[1]
            self._yaw = self._drag_yaw + np.deg2rad(dx * ROT_SENSITIVITY_DEG_PER_PX)
            pitch = self._drag_pitch + np.deg2rad(dy * ROT_SENSITIVITY_DEG_PER_PX)
            self._pitch = float(np.clip(pitch, -np.pi / 2 + 0.01, np.pi / 2 - 0.01))
            self.update()
            self._broadcast_camera()
            event.accept()
            return
        event.ignore()

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.LeftButton and self._drag_active:
            self._drag_active = False
            event.accept()
            return
        event.ignore()

    def wheelEvent(self, event):
        notches = event.angleDelta().y() / 120.0
        self._distance *= ZOOM_FACTOR_PER_NOTCH ** notches
        self._distance = float(np.clip(self._distance, MIN_DISTANCE, MAX_DISTANCE))
        self.update()
        self._broadcast_camera()
        event.accept()

    # ---- rendering ----------------------------------------------------

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing, True)
        painter.fillRect(self.rect(), DEFAULT_BG)

        w, h = self.width(), self.height()
        if w < 4 or h < 4:
            painter.end()
            return

        R = _rotation_matrix(self._yaw, self._pitch)
        cam_offset = np.array([0.0, 0.0, self._distance])
        focal = 0.5 * min(w, h) * (self._distance / max(self._scene_radius, 1.0))

        # iso lobes first (back), then bonds, then atoms (front).
        # Atoms+bonds are intentionally drawn over the lobes — for a
        # translucent orbital this keeps nuclei visible, which is what
        # the teaching context wants. Strict per-primitive depth sorting
        # would give physically-correct occlusion at much higher cost.
        if self._lobes is not None:
            pos_v, pos_f, neg_v, neg_f = self._lobes
            self._draw_tris(painter, pos_v, pos_f, self._pos_color,
                            R, cam_offset, focal, w, h)
            self._draw_tris(painter, neg_v, neg_f, self._neg_color,
                            R, cam_offset, focal, w, h)

        self._draw_bonds(painter, R, cam_offset, focal, w, h)
        self._draw_atoms(painter, R, cam_offset, focal, w, h)

        painter.end()

    # ---- per-primitive ------------------------------------------------

    def _draw_tris(self, painter, verts, faces, base_rgb,
                   R, cam_offset, focal, w, h):
        if verts is None or faces is None or len(faces) == 0:
            return
        eye = (verts - self._target) @ R.T + cam_offset       # (V, 3)
        tri = eye[faces]                                       # (M, 3, 3)

        v0 = tri[:, 0]
        e1 = tri[:, 1] - v0
        e2 = tri[:, 2] - v0
        n = np.cross(e1, e2)
        n_len = np.linalg.norm(n, axis=1)
        n_unit = n / np.clip(n_len, 1e-12, None)[:, None]

        # backface cull: face visible if its normal points roughly toward
        # the camera (negative of view ray from origin to centroid).
        centroid = tri.mean(axis=1)
        view_len = np.linalg.norm(centroid, axis=1)
        view_unit = centroid / np.clip(view_len, 1e-12, None)[:, None]
        facing = -np.einsum('ij,ij->i', n_unit, view_unit)

        in_front = (tri[:, :, 2] > NEAR_Z).all(axis=1)
        keep = (facing > 0) & in_front & (n_len > 1e-12)
        if not keep.any():
            return

        tri = tri[keep]
        n_unit = n_unit[keep]

        depth = tri[:, :, 2].mean(axis=1)
        order = np.argsort(-depth)         # back-to-front
        tri = tri[order]
        n_unit = n_unit[order]

        # Lambert shading
        lambert = np.clip(n_unit @ _LIGHT_DIR, 0.0, 1.0)
        intensity = AMBIENT + DIFFUSE * lambert            # (M,)

        # project to screen
        z = tri[:, :, 2]
        sx = w * 0.5 + focal * tri[:, :, 0] / z
        sy = h * 0.5 - focal * tri[:, :, 1] / z

        rgb = np.clip(base_rgb[None, :] * intensity[:, None], 0, 255).astype(np.int32)
        alpha = int(np.clip(self._opacity * 255, 0, 255))

        painter.setPen(Qt.NoPen)
        for k in range(tri.shape[0]):
            painter.setBrush(QColor(int(rgb[k, 0]), int(rgb[k, 1]),
                                    int(rgb[k, 2]), alpha))
            poly = QPolygonF([QPointF(sx[k, 0], sy[k, 0]),
                              QPointF(sx[k, 1], sy[k, 1]),
                              QPointF(sx[k, 2], sy[k, 2])])
            painter.drawPolygon(poly)

    def _draw_bonds(self, painter, R, cam_offset, focal, w, h):
        if not self._bonds or len(self._atom_pos) == 0:
            return
        eye = (self._atom_pos - self._target) @ R.T + cam_offset
        pen = QPen(BOND_COLOR)
        pen.setCapStyle(Qt.RoundCap)
        for i, j in self._bonds:
            za, zb = eye[i, 2], eye[j, 2]
            if za <= NEAR_Z or zb <= NEAR_Z:
                continue
            sx_a = w * 0.5 + focal * eye[i, 0] / za
            sy_a = h * 0.5 - focal * eye[i, 1] / za
            sx_b = w * 0.5 + focal * eye[j, 0] / zb
            sy_b = h * 0.5 - focal * eye[j, 1] / zb
            z_mid = 0.5 * (za + zb)
            tk = max(mr.BOND_RADIUS * focal / z_mid, 1.5)
            pen.setWidthF(tk)
            painter.setPen(pen)
            painter.drawLine(QPointF(sx_a, sy_a), QPointF(sx_b, sy_b))

    def _draw_atoms(self, painter, R, cam_offset, focal, w, h):
        if len(self._atom_pos) == 0:
            return
        eye = (self._atom_pos - self._target) @ R.T + cam_offset
        order = np.argsort(-eye[:, 2])     # far ones first
        pen = QPen(ATOM_OUTLINE)
        pen.setWidthF(0.5)
        painter.setPen(pen)
        for k in order:
            z = eye[k, 2]
            if z <= NEAR_Z:
                continue
            sx = w * 0.5 + focal * eye[k, 0] / z
            sy = h * 0.5 - focal * eye[k, 1] / z
            rs = self._atoms[k]['radius'] * focal / z
            if rs < 0.5:
                continue
            base = QColor(self._atoms[k]['color'])
            grad = QRadialGradient(QPointF(sx - rs * 0.3, sy - rs * 0.3),
                                   rs * 1.4)
            grad.setColorAt(0.0, _lighten(base, 0.45))
            grad.setColorAt(0.7, base)
            grad.setColorAt(1.0, _darken(base, 0.5))
            painter.setBrush(QBrush(grad))
            painter.drawEllipse(QPointF(sx, sy), rs, rs)

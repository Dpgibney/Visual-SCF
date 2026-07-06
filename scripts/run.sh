#!/usr/bin/env bash
#
# Launch Ryven with the Visual-SCF nodes packages registered.
#
# QT_QPA_PLATFORM=xcb used to be load-bearing because the orbital plotter
# embedded PyVista/VTK, which made raw X11 calls and crashed on Wayland.
# The plotter is now a pure-QPainter widget (plotting2/qpainter3d.py) with
# no native GL surface, so xcb is no longer required for stability. We
# keep the export to match prior behavior; remove it if you've confirmed
# everything renders correctly under your native Qt platform.
#
# Pass any extra args through to ryven (e.g. a saved project file).
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"

QT_QPA_PLATFORM=xcb \
  "$REPO/.venv/bin/ryven" \
  -q "pyside6" \
  -n "$REPO/pyscf" \
  -n "$REPO/plotting2" \
  -n "$REPO/matrices" \
  "$@"

# Visual-SCF

PySCF nodes for the [Ryven](https://github.com/leon-thomm/Ryven) visual
programming editor, with an embedded QPainter-based 3D viewer for atoms and
orbital iso-surfaces. The goal is a hands-on way to teach electronic
structure theory: build a molecule, run an SCF iteration step by step, and
watch molecular orbitals appear as 3D isosurfaces.

## Setup

This repo expects an editable Ryven checkout alongside it so the editor itself can
be patched as Visual-SCF evolves. One-time setup:

```bash
# clone Ryven (and ryvencore separately — it lives in its own repo)
git clone https://github.com/Dpgibney/Ryven.git ~/Projects/Ryven
git clone --branch v0.5.0 https://github.com/leon-thomm/ryvencore.git ~/Projects/ryvencore

# Visual-SCF venv (Python 3.12)
cd ~/Projects/Visual-SCF
python3.12 -m venv .venv
.venv/bin/pip install --upgrade pip wheel setuptools

# editable Ryven installs
.venv/bin/pip install -e ~/Projects/ryvencore
.venv/bin/pip install -e ~/Projects/Ryven/ryvencore-qt
.venv/bin/pip install -e ~/Projects/Ryven/ryven-editor

# runtime deps
.venv/bin/pip install -r requirements.txt
```

## Running

Use the launcher script:

```bash
./scripts/run.sh                       # opens the start dialog
./scripts/run.sh path/to/project.json  # opens a saved project
```

The script still exports `QT_QPA_PLATFORM=xcb` for historical reasons. It
used to be load-bearing because the old orbital plotter embedded VTK, which
makes raw X11 calls and crashes on Wayland. The plotter is now a pure-QPainter
widget — no native GL surface — so the xcb override is no longer required.
You can drop the export from `scripts/run.sh` if you'd like to run under your
native Qt platform.

## Packages

The repo registers three nodes packages with Ryven:

- `pyscf/` — `Molecule`, `Optimize Geometry` (relaxes the structure via
  analytic gradients + geomeTRIC, on the HF or DFT surface), `Guess 1-RDM`,
  `Fock`, `Get MO Coefficients`, `Make 1-RDM`, `SCF Step` (the iterator
  that closes the SCF loop), plus high-level `RHF` / `UHF` / `RKS-DFT` /
  `UKS-DFT` "do everything" nodes and a `CASSCF` node with a selectable
  active space, optional state averaging, and spin-free or spin-resolved
  1-RDM output per root.
- `plotting2/` — `Molecular Orbital Plotter`, `Atomic Orbital Plotter`,
  `LinePlot`, `Print`. The orbital plotters render the molecule and
  iso-surface lobes inline inside the node (left-drag to rotate, scroll
  to zoom); grid resolution is a per-node setting.
- `matrices/` — `Show Matrix`, `MO Coefficient Viewer`.

## Example

Symmetry-broken vs symmetric H₂:

![image](https://github.com/user-attachments/assets/ef062c0c-cdfa-4edc-853a-02bbe6dc23e7)

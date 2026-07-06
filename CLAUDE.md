# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Visual-SCF is a set of three Ryven nodes packages (`pyscf/`, `plotting2/`, `matrices/`) that turn the [Ryven](https://github.com/leon-thomm/Ryven) visual programming editor into a teaching tool for electronic structure theory: build a molecule, run the SCF loop one step at a time, watch MOs appear as 3D isosurfaces. Visual-SCF is *not* an application — it's plugin code loaded by an editable Ryven checkout.

## Running

```bash
./scripts/run.sh                       # opens the start dialog
./scripts/run.sh path/to/project.json  # opens a saved project
```

`scripts/run.sh` invokes `.venv/bin/ryven` with `-n <pkg>` for each of the three packages and exports `QT_QPA_PLATFORM=xcb`. The xcb export is now belt-and-suspenders only — it was load-bearing under the old PyVista plotter (VTK made raw X11 calls that segfaulted on Wayland), but the orbital viewer is now a pure-QPainter widget (`plotting2/qpainter3d.py`) with no native GL surface. Removing the export should be safe; keep it if you want to match the previous behavior exactly.

There are no tests, no linter, and no build step. The launcher is the dev loop.

## Environment setup

Python 3.12 venv at `.venv/`. Ryven is installed editable from sibling clones at `~/Projects/Ryven` and `~/Projects/ryvencore` so the editor itself can be patched alongside Visual-SCF. See README.md for the one-time install incantation.

Pin notes that matter when touching `requirements.txt`:
- `PySide6<6.7` — Ryven 3.5.0's UI breaks on newer PySide.
- `numpy<2` — the shiboken6 wheel bundled with PySide6 6.6 is built against numpy 1.x.

## Architecture

Each of the three top-level directories is a Ryven nodes package, registered separately at launch via `-n` flags. The split is by domain, not layer:

- **`pyscf/`** — quantum chemistry nodes. `MolNode` (Molecule), `GeomOptNode` (Optimize Geometry, via analytic gradients + geomeTRIC), `Guess1RDMNode`, `FockNode`, `GetMOCoeffNode`, `Make1RDMNode`, `SCFStepNode` (the iterator that closes the SCF feedback loop), plus convenience "do everything" nodes `RHFNode` / `UHFNode` / `RKSNode` / `UKSNode`. `MOCoeffViewerNode` also lives here.
- **`plotting2/`** — visualization nodes. `MOPlotNode`, `AOPlotNode`, `LinePlotNode`, `PrintNode`. The orbital plotters render the molecule and iso-surface lobes inline inside the node via `qpainter3d.Painter3DCanvas` (pure QPainter, no GL). `molecule_render.py` is the per-atom/per-bond CPK render data shared by both.
- **`matrices/`** — table viewers. `ShowMOCoeff` and friends, built on `MatrixNodeBase` in `matrices/matrices.py`.

Each package follows the same shape:
- `nodes.py` — pure-logic `Node` subclasses; ends with `export_nodes([...])` and an `@on_gui_load` hook that does `from . import gui`. Nothing in `nodes.py` may import Qt — it must be importable in headless contexts.
- `gui.py` — the Qt widgets (`NodeMainWidget` subclasses) and `@node_gui(NodeClass)` bindings. This is where PySide6 / matplotlib / `qpainter3d` live.
- `mathinspector.py` (pyscf and plotting2 only) — rendered-LaTeX equation panels for the inspector dock. `EQUATIONS` maps node class names to rows of `('t', html)` captions and `('e', latex)` equations; `equation_inspector('<NodeClass>')` builds the `NodeInspectorWidget` subclass that `gui.py` attaches via `inspector_widget_class` + `wrap_inspector_in_default = True`. Equations render through matplotlib **mathtext** (no TeX install; no `\tfrac`, no environments — one line per entry). The renderer half of the two copies is intentionally duplicated (packages can't reliably cross-import); keep them in sync.

Cross-package imports are one-way: other packages' `gui.py` may `from plotting2.qpainter3d import ...` (the Visual-SCF root is on `sys.path` and `plotting2` collides with nothing — `GeomOptWidget`'s trajectory canvas does this), but never import the local `pyscf` package from outside it — that name is shadowed by the PySCF library and only resolves through pkgutil path extension.

Node docstrings are user-facing: Ryven shows `__doc__` as the add-node-list tooltip and in the inspector's description area, both rendered as Qt rich text. Write them as HTML starting with `<p>` (Qt's rich-text sniffing keys on the leading tag), using the Qt subset only (`<p> <b> <i> <tt> <br> <sub> <sup>`, entities). Every exported node must have one.

### The Ryven node contract used here

Conventions repeated across nodes — match them when adding new ones:

- `inputs_ready()` checks `hasattr(self.input(i), 'payload')` (not `is not None`) before dereferencing payloads. Nodes guard `update_event` with this; payloads of `None`-valued inputs would otherwise crash.
- `update_event(self, inp=-1)` is the single entry point for upstream changes. The `inp` arg is the input port index that fired (or `-1`); `SCFStepNode` uses it to distinguish initial-1RDM arrival from feedback.
- GUI updates go through `if hasattr(self, 'gui'): self.gui.<method>(...)` — never assume a GUI exists. `have_gui()` exists as a helper but the explicit `hasattr` is also used inline.
- Node ↔ GUI communication is method calls, not signals. The widget calls `self.node.<method>` for user actions; the node calls `self.gui.<method>` for view updates.
- Payloads are wrapped in `Data(...)` (or `MolData` / `MatrixData`) before `set_output_val`; downstream reads `.payload`.
- Persistence pairs: every stateful node implements `get_state` / `set_state`; `set_state` must default-load missing keys (back-compat with older `testsave.json` projects) and `rebuilt()` is the post-load hook that runs *after* edges are reconnected — use it to re-fire computation, not `set_state`.

### The SCF feedback loop

`SCFStepNode` is the only node that intentionally breaks Ryven's normal forward-only data flow. It has two inputs (`Initial 1-RDM`, `New 1-RDM`) and one output (`Current 1-RDM`). The user wires `Guess1RDM → Initial`, `Make1RDM → New`, `Current → Fock`, closing the loop. `update_event` deliberately *does not* propagate when `New 1-RDM` changes (that would recurse infinitely); progress only happens when the user clicks Step or Reset on the widget. If you change `update_event`, preserve this — the auto-reset on the first `Initial` arrival is the only automatic behavior allowed.

### Isosurface plotting

`plotting2/nodes.py::get_isosurface` and `_eval_mo_grid` are shared by `MOPlotNode` and `AOPlotNode`; the AO plotter calls them with `np.eye(nao)` so column k picks out AO #k. Both accept a `grid_points` kwarg threaded through from the widget's "Grid points" QSpinBox — that's the user-facing perf knob.

The 3D rendering lives in `plotting2/qpainter3d.py::Painter3DCanvas`. It's a `QWidget` with a custom `paintEvent`: a per-frame painter's-algorithm pipeline that vectorizes projection, backface culling, and Lambert shading via numpy, then iterates `QPainter.drawPolygon` for triangles, `drawLine` for bonds, and `QRadialGradient`-filled `drawEllipse` for atoms. Projection is **orthographic** by default (wheel = uniform zoom via the `_distance` state); the class attribute `projection = 'perspective'` restores the pinhole camera, whose distance-coupled focal length makes the wheel a dolly-zoom instead. There is no GL surface — the canvas embeds inside Ryven's `QGraphicsScene` like any other QWidget. Atoms and bonds are drawn *over* the lobes (not depth-sorted against them) — this is intentional so nuclei stay visible through translucent orbitals.

Distances throughout the rendering code are in **Bohr** (matching `pyscf.gto.M` defaults), not Ångström.

## Saved projects

`testsave.json` is a working sample project used during development; it doubles as a back-compat smoke test for `set_state`. If you change a node's persisted schema, load `testsave.json` after to confirm old projects still rebuild — the `built`/`enabled` rename in `MolNode.set_state` is the pattern to follow when fields move.

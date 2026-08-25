# ESP Visualization — VMD

**A lean, independent rendering path for molecular electrostatic potentials,
built on VMD.**

| | |
|---|---|
| Input | Turbomole `pointval` grids (`td.xyz`, `tp.xyz`) + a structure file |
| Output | Gaussian cube files + a ready-to-run VMD scene (`esp.tcl`) |
| Software | VMD 1.9.4 + Python 3 + NumPy, all free |
| Sister project | [`esp_visualization`](https://github.com/FlxLfr/esp_visualization) — the PyMOL pipeline |

> **This project visualises. It does not calculate.**
> No V_S,min, no V_S,max, no σ-hole value, no surface statistics. Those numbers
> are produced by the PyMOL pipeline and are not reimplemented here — a second
> implementation of the same quantity is a duplicate that drifts silently.
> See [section 3](#3-scope-what-is-deliberately-missing).

---

## Contents

1. [Why a second pipeline](#1-why-a-second-pipeline)
2. [What changes when the viewer changes](#2-what-changes-when-the-viewer-changes)
3. [Scope: what is deliberately missing](#3-scope-what-is-deliberately-missing)
4. [Usage](#4-usage)
5. [Installation](#5-installation)
6. [Roadmap](#6-roadmap)
7. [Repository layout](#7-repository-layout)
8. [Notes and troubleshooting](#8-notes-and-troubleshooting)

---

## 1. Why a second pipeline

The PyMOL workflow is finished and does what it was built for: it turns
Turbomole `pointval` output into a reproducible, publication-quality picture of
the ESP on the ρ = 0.001 a.u. density isosurface, and it reports the numbers that
belong to that picture.

There are two reasons to do the picture again in VMD.

**The practical one.** The large majority of published ESP figures in this field
are made with VMD. A workflow is only adopted if it speaks the tool its intended
users already have open, and a figure is only comparable to the literature if it
can be produced the way the literature produced it. A pipeline that exists only
for PyMOL asks every reader to switch software before they can use it.

**The methodological one.** Two independent viewers reading the *same* cube files
are a control on each other. Isosurface extraction, interpolation of the
potential onto that surface, and the mapping of values to colour are all
implementation details, and every one of them can change what a figure shows. If
the same data at the same isovalue on the same colour scale produces the same
picture in both, that agreement is evidence that the figure shows a property of
the quantum-chemical data rather than an artefact of one renderer. If the two
disagree, the disagreement is itself worth reporting — and a single-viewer study
would never have surfaced it.

That makes this repository more than a port. The comparison is part of the
result.

---

## 2. What changes when the viewer changes

| | PyMOL | VMD |
|---|---|---|
| Scene language | Python / `.pml` | Tcl / `.tcl` |
| ESP on surface | `ramp_new` + `set surface_color` | second volumetric dataset via `mol addfile`, `Isosurface` rep coloured by `Volume` |
| Colour scale | ramp with explicit min/max | `color scale method RWB` + `mol scaleminmax` |
| Structure file | loaded separately | not needed — the cube carries the atoms |
| Batch mode | `pymol -cq script.py --` | `vmd -dispdev text -e script.tcl` |
| Renderer | built-in ray tracer | Tachyon / TachyonLOptiXInternal |

The essential mechanic is the second row: **VMD has no ramp object.** Both grids
are loaded into the *same* molecule ID — the density first, the potential as an
additional volumetric dataset — and the isosurface drawn from volume 0 is then
coloured by the values of volume 1. Everything else follows from that.

The fourth row is a quiet win. In PyMOL the structure file is a second source of
coordinates, and if its units disagree with the grid, the molecule floats next to
its own surface. A cube file carries its atom block in Bohr alongside the grid,
so in VMD there is only one source and nothing to misalign.

---

## 3. Scope: what is deliberately missing

`xyzToCubeToVMDVis.py` converts grids and writes a scene. It computes no physical
quantity, and that is a decision rather than an omission: V_S,min, V_S,max and the
σ-hole potential are already implemented, parameter-studied and documented in the
PyMOL project. A second implementation here would be a duplicate that nobody
notices drifting until two chapters of the same thesis quote different numbers
for the same molecule.

The most visible consequence is the **colour scale**. It is not derived from the
data; it is a parameter (`--esp-range`, default 0.035 a.u.). That is the point:
only a scale fixed from outside and applied identically in both pipelines makes
their images comparable at all. Take the value from the PyMOL run's
`<molecule>_settings.txt`, line "Farbskala".

What is carried over unchanged, so the two image sets can be laid side by side:

* **The surface.** ρ = 0.001 a.u. electron-density isosurface (Politzer/Murray).
* **The colour convention.** Red negative, blue positive, white at zero.
* **The three views.** π face, in-plane edge, σ-hole along the C–X axis.
* **`reference/` is committed and not edited; `sandbox/` is ignored.**

---

## 4. Usage

```bash
python scripts/xyzToCubeToVMDVis.py --struct <structure> td.xyz tp.xyz [options]
vmd -e esp.tcl
```

Concretely, from inside a molecule folder:

```bash
python ../../scripts/xyzToCubeToVMDVis.py \
    --struct brombenzol_aro_opti.mol td.xyz tp.xyz --esp-range 0.035
```

This writes `td.cube`, `tp.cube` and `esp.tcl` next to the input grids. Expect a
minute or two per grid — the pointval files are ~1.25 GB each and the cubes come
out around 200 MB.

**Options that change the cube files**

| | |
|---|---|
| `--struct`, `-s` | structure file for the cube's atom block: `.xyz`, `.mol`, `.sdf` (required) |
| `--struct-unit` | `angstrom` (default) or `bohr`; ignored for `.mol`/`.sdf`, which are Å by definition |
| `--stride N` | keep only every N-th grid point — `2` gives an 8× smaller file |
| `--outdir`, `-o` | write elsewhere than next to the input |

**Options that change only `esp.tcl`**

| | |
|---|---|
| `--tcl-only` | rewrite `esp.tcl` from the existing cubes, touching no grid |
| `--esp-range` | half-width of the colour scale in a.u. (default 0.035) |
| `--iso` | isovalue of the density surface (default 0.001) |
| `--opacity` | surface opacity 0…1 (default 0.50 — see section 8 on how it relates to PyMOL's) |
| `--no-colorbar` | do not show the colour bar on start |
| `--scale` | zoom of the start view (default 0.12; useful range 0.08–0.18) |
| `--color-scale` | VMD colour scale, default `RWB`; `BWR` reverses it |
| `--no-vmd` | cubes only, no scene |

**Never re-run the conversion just to change the scene.** Converting the grids
takes minutes; `--tcl-only` reads the two cube headers, works out which is which
and rewrites the scene in well under a second. `--struct` is not needed then:

```bash
python ../../scripts/xyzToCubeToVMDVis.py td.cube tp.cube --tcl-only --opacity 0.5
```

**Inside VMD**, `esp.tcl` defines four commands:

| | |
|---|---|
| `esp_view pi \| edge \| sigma` | the three standard views |
| `esp_iso <value>` | change the isovalue live |
| `esp_range <half-width>` | change the colour scale live |
| `esp_opacity <0…1>` | change the surface opacity live |
| `esp_colorbar [gap] [height] [textsize]` | redraw / move the colour bar |
| `esp_colorbar_off` | remove the colour bar |
| `esp_snapshot <name>` | ray-traced image via Tachyon |

Tune with the live commands first, then bake the value you settled on into the
file with `--tcl-only` so the scene reproduces itself.

---

## 5. Installation

### VMD

Download from the [TCBG test release page](https://www.ks.uiuc.edu/Research/vmd/alpha/)
(free registration) and take **Version 1.9.4, "Windows 64-bit, CUDA, OptiX,
OSPray"**.

Two notes on that choice:

* The page labels the Windows build "Windows 10". That is the tested minimum, not
  an upper limit — it installs and runs on Windows 11 natively.
* Version 2.0.0 exists as a monthly-updated alpha with a *"fully overhauled plugin
  system"* still to come. Its new features do not touch what this pipeline needs,
  and a viewer rebuilt every month is a poor foundation for a workflow that has
  to be reproducible. 1.9.x is what the published literature used. To try 2.0,
  install it alongside — the two do not conflict.

The Windows installer does **not** put VMD on the `PATH`, so `vmd -e esp.tcl`
will not work in a fresh shell until you add it — see section 8.

### Python

Only NumPy is required.

```bash
conda env create -f environment.yml
conda activate esp-vmd
```

The `esp` environment from the PyMOL project works just as well — it already has
NumPy, and VMD is not a conda package either way.

---

## 6. Roadmap

- [x] Install VMD 1.9.4, confirm `Extensions → Tk Console` opens
- [x] Load a cube pair by hand, get one ESP-coloured isosurface on screen
- [x] Confirm `mol addfile` + colour-by-`Volume` is the right mechanic
- [x] `xyzToCubeToVMDVis.py` — pointval → cube → `esp.tcl`
- [ ] Convert the halobenzene set and check all three views
- [ ] `scripts/render_esp.tcl` — the standard image set, headless
- [ ] Fill `reference/` with a known-good example rendered in VMD
- [ ] **Cross-check PyMOL vs. VMD** on identical cubes at an identical colour scale

---

## 7. Repository layout

```
VMD_esp_visualization/
├── README.md                     this document
├── environment.yml               conda environment (Python side only)
├── .gitignore
├── scripts/
│   ├── xyzToCubeToVMDVis.py      pointval -> cube -> esp.tcl
│   ├── constants.py              unit conversions
│   └── render_esp.tcl            standard image set, headless      (planned)
├── reference/                    known-good example — output, not input
├── results/                      the delivered image sets
├── tools/                        helper scripts
├── docs/                         background and the PyMOL/VMD comparison
└── sandbox/                      your own data and experiments, not tracked
```

**Large files are deliberately not tracked.** `.gitignore` excludes `*.cube`,
`*.dx` and the raw `td.xyz`/`tp.xyz` grids — a full-resolution cube is ~200 MB
and GitHub rejects anything above 100 MB. Regenerate them with
`xyzToCubeToVMDVis.py`. `esp.tcl` is generated too and is ignored inside molecule
folders; only the scripts under `scripts/` are tracked.

**`reference/` and `sandbox/` do different jobs.** `reference/` holds a
known-good example and is committed; `sandbox/` is where your own molecules and
the large raw data live, and git ignores it entirely. If a run goes wrong,
`reference/` tells you whether the problem is your installation or your data.

The name `xyzToCubeToVMDVis.py` is deliberately not `xyzToCube.py`: the PyMOL
project has a script by that name with a different feature set, and two files
with one name in two repositories is how the wrong one gets edited.

---

## 8. Notes and troubleshooting

**`vmd` is not recognised as a command (Windows).**
The installer registers a Start Menu entry but does not extend the `PATH`. Find
the executable — it is `vmd.exe` in `C:\Program Files\University of Illinois\VMD`
or the `Program Files (x86)` equivalent — and add that folder once, in PowerShell:

```powershell
$vmd = "C:\Program Files (x86)\University of Illinois\VMD"   # adjust if needed
$old = [Environment]::GetEnvironmentVariable("Path", "User")
[Environment]::SetEnvironmentVariable("Path", "$old;$vmd", "User")
```

Open a new shell afterwards. Until then, start VMD from the Start Menu and load
the scene from the Tk Console instead:

```tcl
cd {C:/Users/felix/.../sandbox/brombenzol}
source esp.tcl
```

Forward slashes and braces, not backslashes — Tcl reads `\` as an escape
character.

**"Unable to load molecule" when opening an `.xyz` in VMD.**
The Turbomole-derived structure files here have no XYZ header — no atom count,
no comment line, just coordinates. PyMOL guesses; VMD's reader does not and
aborts. You do not need to fix the file: load the **cube** instead. It carries
the atom block, so `mol new td.cube type cube` gives you geometry and grid at
once.

**The molecule sits small and off-centre in the window.**
`display resetview` fits the view to everything drawn, and the isosurface carries
the extent of the whole grid box (30 × 30 × 30 Bohr), inside which the molecule
is rarely centred. `esp.tcl` handles this in `esp_center`, which centres on the
atoms instead. If you lose the view, call `esp_view pi`.

**The skeleton does not show through at all, at any opacity.**
This one is not about opacity. **VMD draws representations in index order and
writes depth for transparent surfaces too.** With the isosurface on rep 0 it is
drawn first, and the skeleton drawn afterwards fails the depth test everywhere
the surface lies in front of it — which is everywhere. You get a visibly
translucent surface with nothing behind it, and lowering the opacity does not
help, because the geometry was discarded before blending. The tell-tale symptom:
the skeleton appears only where the near clipping plane cuts the surface open
when you zoom in.

The fix is the order. `esp.tcl` adds the opaque Licorice rep **first** (rep 0)
and the transparent isosurface **second** (rep 1), which is why the procs address
the surface as `$REP_SURF`. If you build a scene by hand in the GUI, the same
rule applies: transparent reps last.

**The skeleton shows through but is too faint.**
Two effects stack here, and both are worth knowing when comparing to PyMOL.

*You look through two layers.* The isosurface is closed: the line of sight
enters at the front and leaves at the back. At opacity *a* only (1−*a*)² of what
lies behind survives — at *a* = 0.6 that is 16 %, and the skeleton disappears
although the number sounds like "half transparent". PyMOL's `transparency 0.15`
is nominally *a* = 0.85, yet its images show the sticks; the two programs weight
transparency differently. **Identical parameters do not give identical images**,
which is exactly the kind of renderer detail this project exists to document.
Around 0.3 in VMD gives the impression PyMOL gives at 0.15.

*A strong specular highlight fills the gap back in.* A milky sheen sits on top of
whatever the opacity lets through, so `esp.tcl` keeps specular low (0.10).

**The colour bar sits off-screen, or is the wrong size.**
VMD has no legend of its own, so `esp.tcl` draws one from graphics primitives in
a separate, empty molecule and detaches it from the mouse with `mol fix` — else
it would tilt away with the molecule.

Its placement is **relative to the molecule, in ångström**: below the molecule's
bounding sphere, sharing the molecule's centring and scale. That matters. Placing
it in absolute screen coordinates requires knowing VMD's visible world extent,
which depends on window size, aspect ratio and zoom — guess it and the bar lands
just outside the frame with no error message to tell you. Anchored to the
molecule, it cannot miss. The start view is zoomed out slightly (`--scale 0.10`)
to leave room for it.

Its three arguments are the gap to the molecule (Å), the bar height (Å) and the
label size:

```tcl
esp_colorbar 1.5          ;# further below the molecule
esp_colorbar 0.5 1.0      ;# thicker bar
esp_colorbar 0.5 0.7 1.2  ;# larger labels
```

`esp_range` and `esp_view` redraw the bar automatically, so the labels always
match the scale actually in use and the bar keeps up with a change of view.

**Zooming in cuts a hole in the surface.**
That is the near clipping plane, not transparency — VMD's default (0.5) slices a
disc out of the surface as you approach, and you find yourself looking into a
cross-section. `esp.tcl` sets `display nearclip set 0.01`. If you have an old
scene open, run that line in the Tk Console.

**The surface looks grainy, like a bad print.**
VMD's default render mode fakes transparency by dropping pixels ("screen door").
`esp.tcl` requests `display rendermode GLSL`, which does real alpha blending; if
your driver cannot, the request fails silently and you keep the dotted look on
screen. The rendered images are unaffected — Tachyon does proper transparency
regardless.

**VMD shows two frames after loading.**
The second cube brings its own copy of the coordinates, which VMD appends as an
extra frame. `esp.tcl` deletes it.

**The blue σ-hole looks pale.**
It should. On a ±0.035 a.u. scale a σ-hole of +0.016 a.u. reaches under half
saturation. Rescaling to make it vivid would break comparability with the PyMOL
images — that is what `esp_range` is for when you want to look, and why the
delivered images keep the fixed scale.

**The three standard views are rotated wrong.**
`esp_view` assumes a planar molecule in the xy-plane with the C–X axis along y,
which is how the Turbomole optimisation leaves the halobenzenes here. For other
geometries, rotate by hand and use `esp_snapshot`.

# ESP Visualization — VMD

**A lean, independent rendering path for molecular electrostatic potentials,
built on VMD.**

| | |
|---|---|
| Input | Turbomole `pointval` grids (`td.xyz`, `tp.xyz`) + a structure file |
| Output | cube files, a ready-to-run VMD scene (`esp.tcl`), the image set in `images/` |
| Software | VMD 1.9.4 + Python 3 + NumPy, all free |
| Sister project | [`Pymol_esp_visualization`](https://github.com/FlxLfr/Pymol_esp_visualization) — the PyMOL pipeline |

> **This project visualises.** It computes V_S,min and V_S,max, because without
> them there is no colour scale and therefore no meaningful picture. It does not
> compute the σ-hole potential and does not interpolate onto the isosurface —
> those live in the PyMOL pipeline. See
> [section 3](#3-scope-what-is-deliberately-missing).

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

The line runs between *what the picture needs* and *what the study reports*.

**Computed here:** V_S,min and V_S,max on the ρ = 0.001 a.u. shell, and the
symmetric colour range rounded up from them. Without those there is no colour
scale, and a picture with an arbitrary scale is not a result. `--esp-range auto`
is the default; a number overrides it.

The values come from **grid points near the shell**, not from an interpolated
surface, and no σ-hole is located. That is the deliberate cut: locating the
σ-hole maximum needs trilinear interpolation onto the isosurface, and that code
is already written, parameter-studied and documented in the PyMOL project. A
second implementation here would be a duplicate nobody notices drifting until
two chapters of the same thesis quote different numbers for the same molecule.
**Quote the PyMOL numbers.** The ones here size a colour bar.

The two agree where they overlap. On the committed reference dataset
(4-bromoacetophenone) this pipeline finds 173 shell points,
V_S,min = −0.06379 a.u., V_S,max = +0.04693 a.u., range ±0.0650 — identical to
every digit in the PyMOL project's `reference/summary.csv`.

The computed values are stamped into the first line of `tp.cube`, so
`--tcl-only` recovers them without re-reading a 200 MB grid.

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
    --struct brombenzol_aro_opti.mol td.xyz tp.xyz
```

This writes `td.cube`, `tp.cube` and `esp.tcl` next to the input grids. Expect a
minute or two per grid — the pointval files are ~1.25 GB each and the cubes come
out around 200 MB.

**Options that change the cube files**

| | |
|---|---|
| `--struct`, `-s` | structure file for the cube's atom block: `.xyz`, `.mol`, `.sdf` (required unless `--tcl-only`) |
| `--struct-unit` | `angstrom` (default) or `bohr`; ignored for `.mol`/`.sdf`, which are Å by definition |
| `--stride N` | keep only every N-th grid point — `2` gives an 8× smaller file |
| `--outdir`, `-o` | write elsewhere than next to the input |

**Options that change only `esp.tcl`**

| | |
|---|---|
| `--tcl-only` | rewrite `esp.tcl` from the existing cubes, touching no grid |
| `--esp-range` | half-width of the colour scale in a.u., or `auto` (default): from V_S,min/V_S,max on the shell |
| `--iso` | isovalue of the density surface (default 0.001) |
| `--opacity` | surface opacity 0…1 (default 0.50 — see section 8 on how it relates to PyMOL's) |
| `--scale` | zoom, a number or `auto` (default): from the molecule's size and the window height |
| `--fill` | fraction of the image height the molecule fills at `--scale auto` (default 0.85) |
| `--color-scale` | VMD colour scale, default `RWB`; `BWR` reverses it |
| `--no-vmd` | cubes only, no scene |

**Never re-run the conversion just to change the scene.** Converting the grids
takes minutes; `--tcl-only` reads the two cube headers, works out which is which
and rewrites the scene in well under a second. `--struct` is not needed then:

```bash
python ../../scripts/xyzToCubeToVMDVis.py td.cube tp.cube --tcl-only --opacity 0.5
```

### Everything in one go

```bash
python scripts/run_allVMD.py --root sandbox
```

Walks the tree, and for every molecule folder converts the grids, writes
`esp.tcl`, renders the three views, builds the colour bar and drops a
`summary.csv` at the root.

Without `--root` it runs on `reference/` instead — that is the smoke test, see
[section 5](#5-installation).

**It analyses first and renders once.** V_S,min and V_S,max come out of the cube
files, so no rendering is needed to find them: the script collects the values of
all molecules first, then renders a single pass.

**No second pass — a recommendation instead.** Images are only comparable on one
shared scale, and the PyMOL pipeline renders twice to get there (`--two-pass`).
In VMD an image costs a minute or two, so the first set would be waste. Every
run therefore ends with the smallest scale that covers all molecules, and the
ready-made command line to render with it:

```
Jedes Molekuel hat seine eigene Skala - 0.0350, 0.0500, 0.0700 a.u.
Nebeneinanderlegen darf man die Bilder so NICHT.
Kleinster Wert, der alle abdeckt: +/- 0.0700 a.u.
Zum Rendern eines vergleichbaren Satzes:

    python run_allVMD.py --root ../sandbox --esp-range 0.0700
```

You read it, decide what the figure should show, and start the run that counts —
the smallest scale that clips nothing is not always the most informative one. A
single molecule with a very negative oxygen flattens the contrast of every other
molecule in the set, and that is a judgement call, not something the script
should make for you. Set a fixed value that is too small and the run says so:
which molecules it clips, and what they would need.

| | |
|---|---|
| `--esp-range auto` | each molecule on its own scale (default) |
| `--esp-range 0.035` | fixed value — this is what you pass on the second run, and what the PyMOL `summary.csv` value goes into for a direct comparison |
| `--esp-range common` | one scale for all, taken straight from this run — the recommendation applied without asking |
| `--only NAME…` | restrict to these folders, wildcards allowed |
| `--stride`, `--struct-unit`, `--force-convert` | conversion, as in the single-molecule script |
| `--iso`, `--opacity`, `--scale`, `--fill`, `--color-scale` | scene |
| `--headless` | no VMD window (see below) |
| `--res`, `--ao`, `--shadows`, `--keep-tga`, `--dpi`, `--vmd` | rendering |
| `--no-render` | convert and write `esp.tcl` only |
| `--images-dir`, `--summary` | output locations |

`summary.csv` carries the numbers *and* what each image was rendered with —
resolution, renderer, ambient occlusion, shadows, and per view which pass
produced it. When Tachyon falls back, that must not disappear from the record.

### The image set

From a molecule folder, after `esp.tcl` exists:

```bash
python ../../scripts/render_espVMD.py
```

Renders the three views with Tachyon (ambient occlusion and shadows on),
converts VMD's TGA output to PNG, and writes the colour bar — into `images/`:

```
images/<molecule>_pi.png  _edge.png  _sigma.png  _colorbar.png  _settings.txt
```

The same five files the PyMOL pipeline produces, deliberately: only then can the
two sets be laid side by side. Options: `--res 2400x1920`, `--ao` (slower,
deeper), `--no-vmd` (only convert and build the bar), `--vmd <path>` if VMD is
not on the `PATH`.

**`--headless` keeps the window shut.** VMD normally opens its OpenGL window for
every pass, which over nine molecules is a dozen windows stealing focus.
`--headless` starts it with `-dispdev text`: no window, no OpenGL. Tachyon still
renders — it works from the scene graph, not from the framebuffer — and the
image size then comes from `-size` on the command line rather than
`display resize`, so it is no longer capped by your screen. On a 1061-pixel-tall
screen, `--res 1600x1280` finally gives 1280.

The one thing that genuinely needs the window is the snapshot fallback, which
copies the OpenGL buffer. So `--headless` does not remove it, it only defers it:
the window appears solely when Tachyon has failed on a view and the alternative
would be losing that image. On the current data that happens for roughly a third
of the views.

**The colour bar comes from matplotlib, not from VMD** — and not from PyMOL
either. Neither program has a legend object, so both pipelines draw it the same
way, as a separate PNG. The scene itself deliberately draws none: doing it in VMD
means putting graphics primitives into a dummy molecule and detaching it from the
mouse, which is fiddly and never looks quite right.

**Ambient occlusion is the single biggest step** from preview to publication
image. It barely shows in the OpenGL window and strongly in the ray-traced
render, laying soft shadows into the recesses between the CH bulges. Compare a
VMD OpenGL screenshot against a PyMOL `ray` image and you are comparing a
preview against a render — not the two programs.

**Inside VMD**, `esp.tcl` defines these commands:

| | |
|---|---|
| `esp_view pi \| edge \| sigma` | the three standard views |
| `esp_iso <value>` | change the isovalue live |
| `esp_range <half-width>` | change the colour scale live |
| `esp_opacity <0…1>` | change the surface opacity live |
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

NumPy for the conversion; matplotlib and Pillow for the image set.

```bash
conda env create -f environment.yml
conda activate esp-vmd
```

The `esp` environment from the PyMOL project works just as well — it already has
NumPy, and VMD is not a conda package either way.

### Verify & smoke test

Check the Python side first:

```bash
python -c "import numpy, matplotlib, PIL; print('ok')"
```

Then run the whole chain on the bundled reference molecule — `run_allVMD.py`
with no arguments does exactly that:

```bash
cd scripts
python run_allVMD.py
```

It converts `reference/brombenzol/td.xyz` and `tp.xyz` to cubes, writes
`esp_check.tcl`, renders the three views into `images_check/` and leaves a
`reference/summary_check.csv`. The `_check` names are deliberate: the smoke test
must never overwrite the committed reference files, otherwise you can no longer
tell whether the reference is still the reference.

Compare your `images_check/` with the committed `images/`, and your numbers with
`reference/summary.csv`. On the decimated reference grid expect

| | |
|---|---|
| V_S,min | −0.01863 a.u. |
| V_S,max | +0.03070 a.u. |
| colour scale | ±0.0350 a.u. |
| shell points | 313 |

**These are the same numbers the PyMOL pipeline reports**, because both
reference sets are built from the same bromobenzene data with the same
parameters (`--stride 5`, `--margin 2.5`) and both derive V_S from grid points
near the ρ = 0.001 shell. That is the point: it is a cross-check of the two
pipelines, not just of this one. If your numbers match, the installation is
sound and any remaining difference between the two image sets is the viewer,
not the data.

The reference grid is decimated to 0.60 Bohr — five times coarser than the
delivered data — and is far too coarse for a V_S value you would quote. It
answers "does the installation run and do the documented numbers come out", not
"how big is the σ-hole". Rebuild it, or make one for another molecule, with:

```bash
python tools/make_reference.py sandbox/brombenzol --name brombenzol
```

It crops to the bounding box of the ρ > iso/2 region plus a margin, so the whole
isosurface stays inside the grid, and then keeps every n-th point (`--stride`).
It ships as raw `pointval` files rather than cubes on purpose: unit conversion
and index reordering are the steps most likely to break, and ready-made cubes
would skip exactly those.

---

## 6. Roadmap

- [x] Install VMD 1.9.4, confirm `Extensions → Tk Console` opens
- [x] Load a cube pair by hand, get one ESP-coloured isosurface on screen
- [x] Confirm `mol addfile` + colour-by-`Volume` is the right mechanic
- [x] `xyzToCubeToVMDVis.py` — pointval → cube → `esp.tcl`
- [x] Colour scale from V_S,min / V_S,max, cross-checked against the PyMOL reference
- [ ] Convert the halobenzene set and check all three views
- [x] `render_esp.py` / `render_esp.tcl` — the standard image set into `images/`
- [x] `run_allVMD.py` — batch driver, analyse-then-render
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
│   ├── xyzToCubeToVMDVis.py      pointval -> cube -> esp.tcl   (535 lines)
│   ├── esp_template.tcl          the VMD scene, with @@placeholders@@
│   ├── render_espVMD.py          image set -> images/  (drives VMD)
│   ├── render_esp.tcl            the VMD half: 3 views via Tachyon
│   ├── run_allVMD.py             batch: convert + scene + images + summary.csv
│   └── constants.py              unit conversions
├── reference/                    known-good example — output, not input
│   ├── summary.csv
│   └── brombenzol/
│       ├── brombenzol_aro_opti.mol
│       ├── td.xyz                raw pointval grids, decimated to 0.60 Bohr
│       ├── tp.xyz                (32×37×24, 1.7 MB each)
│       └── images/               reference images
├── results/                      the delivered image sets
├── tools/
│   └── make_reference.py         decimated reference set from a full folder
├── docs/                         background and the PyMOL/VMD comparison
└── sandbox/                      your own data and experiments, not tracked
```

**The scene lives in its own file.** `esp_template.tcl` is real Tcl with
`@@placeholder@@` markers that the Python script fills in. Keeping 240 lines of
Tcl inside a Python string cost syntax highlighting, made the converter twice as
long as the job warrants, and hid mistakes in the escaping.

**Large files are deliberately not tracked.** `.gitignore` excludes `*.cube`,
`*.dx` and the raw `td.xyz`/`tp.xyz` grids — a full-resolution cube is ~200 MB
and GitHub rejects anything above 100 MB. Regenerate them with
`xyzToCubeToVMDVis.py`. `esp.tcl` is generated too and is ignored inside molecule
folders; only the scripts under `scripts/` are tracked.

**`reference/` and `sandbox/` do different jobs.** `reference/` holds a
known-good example and is committed; `sandbox/` is where your own molecules and
the large raw data live, and git ignores it entirely. If a run goes wrong,
`reference/` tells you whether the problem is your installation or your data.
Its decimated grids are the one exception to the rule above and are tracked on
purpose — that is what makes a fresh clone self-testing. The smoke test's own
output (`images_check/`, `esp_check.tcl`, `summary_check.csv`) is ignored again.

The name `xyzToCubeToVMDVis.py` is deliberately not `xyzToCube.py`: the PyMOL
project has a script by that name with a different feature set, and two files
with one name in two repositories is how the wrong one gets edited.

---

## 8. Notes and troubleshooting

**`vmd` is not recognised as a command (Windows).**
The installer registers a Start Menu entry but does not extend the `PATH`. Find
the executable and add its folder once, in PowerShell:

```powershell
$sh = New-Object -ComObject WScript.Shell     # ask the Start Menu shortcut
Get-ChildItem "$env:ProgramData\Microsoft\Windows\Start Menu\Programs" -Recurse `
  -Filter "*VMD*.lnk" -ErrorAction SilentlyContinue |
  ForEach-Object { $sh.CreateShortcut($_.FullName).TargetPath }

$vmd = "C:\Program Files\VMD"               # the folder the line above reported
$old = [Environment]::GetEnvironmentVariable("Path", "User")
if ($old -notlike "*$vmd*") {
    [Environment]::SetEnvironmentVariable("Path", "$old;$vmd", "User")
}
```

Take `$old` rather than appending to `$env:Path`, which would copy every system
entry into your user PATH as well. `render_esp.py` also looks at `VMDDIR`, which
the installer does set, so it usually finds VMD even without this.

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

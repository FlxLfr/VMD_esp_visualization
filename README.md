# ESP Visualization — VMD

**A second, independent rendering path for molecular electrostatic potentials,
built on VMD.**

| | |
|---|---|
| Input | the same Gaussian cube files the PyMOL pipeline uses (`td.cube`, `tp.cube`) |
| Output | the same standard image set (`*_pi`, `*_sigma`, `*_edge`, `*_colorbar`) plus `*_settings.txt` |
| Software | VMD 1.9.4 + Python 3 + NumPy, all free |
| Sister project | [`esp_visualization`](https://github.com/FlxLfr/esp_visualization) — the PyMOL pipeline |

> **Status: scaffold.** This repository currently holds the folder structure,
> the ignore rules and the viewer-independent conversion scripts copied from the
> PyMOL pipeline. The VMD side — `esp.tcl`, `render_esp.tcl`, the batch driver —
> is not written yet. Section 4 is the roadmap.

---

## Contents

1. [Why a second pipeline](#1-why-a-second-pipeline)
2. [What changes when the viewer changes](#2-what-changes-when-the-viewer-changes)
3. [What stays fixed](#3-what-stays-fixed)
4. [Roadmap](#4-roadmap)
5. [Installation](#5-installation)
6. [Repository layout](#6-repository-layout)
7. [Setting up the git remote](#7-setting-up-the-git-remote)

---

## 1. Why a second pipeline

The PyMOL workflow is finished and does what it was built for: it turns
Turbomole `pointval` output into a reproducible, publication-quality picture of
the ESP on the ρ = 0.001 a.u. density isosurface, without a single manual step.

There are two reasons to do the same thing again in VMD.

**The practical one.** The large majority of published ESP figures in this field
are made with VMD. A workflow is only adopted if it speaks the tool its intended
users already have open, and a result is only checkable against the literature
if it can be produced the way the literature produced it. A pipeline that exists
only for PyMOL asks every reader to switch software before they can use it.

**The methodological one, and the more interesting of the two.** Two independent
viewers reading the *same* cube files are a control on each other. Isosurface
extraction, interpolation of the potential onto that surface, and the mapping of
values to colour are all implementation details, and every one of them can shift
a number or a picture. If the σ-hole maximum and the visual impression agree
between PyMOL and VMD, that agreement is evidence that what we are looking at is
a property of the quantum-chemical data rather than an artefact of one renderer.
If they disagree, the disagreement is itself a finding worth reporting — and one
that a single-viewer study would never have surfaced.

That makes this repository more than a port. The comparison is part of the
result.

---

## 2. What changes when the viewer changes

| | PyMOL | VMD |
|---|---|---|
| Scene language | Python / `.pml` | Tcl / `.tcl` |
| ESP on surface | `ramp_new` + `set surface_color` | second volumetric dataset via `mol addfile`, `Isosurface` rep coloured by `Volume` |
| Colour scale | ramp with explicit min/max | `mol scaleminmax`, colour scale `BWR` |
| Batch mode | `pymol -cq script.py --` | `vmd -dispdev text -e script.tcl` |
| Renderer | built-in ray tracer | Tachyon / TachyonLOptiXInternal |
| Cube units | needed care (`--struct-unit`) | read natively, Bohr assumed by the format |

The essential mechanic to get right is the first one in the table: VMD does not
have a ramp object. Both grids are loaded into the *same* molecule ID — the
density first, the potential as an additional volumetric dataset — and the
isosurface drawn from grid 0 is then coloured by the values of grid 1. Everything
else in the pipeline follows from that.

---

## 3. What stays fixed

These are carried over unchanged from the PyMOL pipeline, deliberately, so that
the two image sets can be laid next to each other:

* **The surface.** ρ = 0.001 a.u. electron-density isosurface. Rationale and
  literature in `docs/ESP_Visualization_Background.docx` of the sister repo.
* **The image set.** `*_pi.png`, `*_sigma.png`, `*_edge.png`, `*_colorbar.png` —
  same views, same names.
* **The colour scale.** Same numeric min/max across a molecule set, so colours
  mean the same thing in every picture.
* **`*_settings.txt` next to every image.** An image never travels without the
  parameters that produced it.
* **`reference/` is committed and not edited; `sandbox/` is ignored.** Same split,
  same reasoning as upstream.

---

## 4. Roadmap

- [ ] Install VMD 1.9.4 and confirm `Extensions → Tk Console` opens
- [ ] Load a cube pair by hand, get one ESP-coloured isosurface on screen
- [ ] Verify that `mol addfile` + colour-by-`Volume` reproduces the PyMOL picture
- [ ] `scripts/esp.tcl` — generated per-molecule interactive scene (the `esp.pml` analogue)
- [ ] `scripts/render_esp.tcl` — the standard image set, headless
- [ ] Teach `xyzToCube.py` to emit `esp.tcl` instead of `esp.pml`
- [ ] `scripts/run_all.py` — batch driver + `summary.csv`
- [ ] Fill `reference/` with the 4-bromoacetophenone example, rendered in VMD
- [ ] **Cross-check PyMOL vs. VMD** on identical cubes: σ-hole value, colour scale, visual agreement

---

## 5. Installation

### VMD

Download from the [TCBG test release page](https://www.ks.uiuc.edu/Research/vmd/alpha/)
(free registration required) and take **Version 1.9.4, "Windows 64-bit, CUDA,
OptiX, OSPray"**.

Two notes on that choice:

* The page labels the Windows build "Windows 10". That is the tested minimum,
  not an upper limit — it installs and runs on Windows 11 natively, no
  compatibility mode.
* Version 2.0.0 exists as a monthly-updated alpha with a *"fully overhauled
  plugin system"* still to come. Its new features (UI redesign, faster surfaces,
  live ray tracing) do not touch what this pipeline needs, and a viewer that is
  rebuilt every month is a poor foundation for a workflow that has to be
  reproducible. 1.9.x is what the published literature used. If you want to try
  2.0, install it alongside — the two do not conflict.

Keep the default install path so that `vmd` is resolvable from the command line;
the batch driver depends on it.

### Python part

```bash
conda env create -f environment.yml
conda activate esp-vmd
```

Only NumPy is required — the heavy lifting happens inside VMD.

---

## 6. Repository layout

```
VMD_esp_visualization/
├── README.md                     this document
├── environment.yml               conda environment (Python side only)
├── .gitignore
├── scripts/
│   ├── xyzToCube.py              Turbomole pointval -> Gaussian cube   [copied]
│   ├── constants.py              unit conversions                      [copied]
│   ├── ansi.py                   console colours                       [copied]
│   ├── UPSTREAM.txt              where the copies came from, and the no-drift rule
│   ├── esp.tcl                   interactive scene                     (planned)
│   ├── render_esp.tcl            standard image set, headless          (planned)
│   └── run_all.py                batch driver + summary.csv            (planned)
├── reference/                    known-good example — output, not input
├── results/                      the delivered image sets
├── tools/                        helper scripts (own environment)
├── docs/                         background, method notes, the PyMOL/VMD comparison
└── sandbox/                      your own data and experiments, not tracked
```

**Three scripts are copies, not forks.** `xyzToCube.py`, `constants.py` and
`ansi.py` are byte-identical to commit `ed3b0d0` of the PyMOL repository. Grid
conversion is viewer-independent, and copying keeps this repository able to clone
and run on its own. The rule that keeps the copies honest is in
`scripts/UPSTREAM.txt`: fix them upstream, then re-copy — do not edit them here.
The single planned exception is the scene-file writer inside `xyzToCube.py`,
which has to emit Tcl instead of `.pml`; at that point the file becomes a real
fork and `UPSTREAM.txt` says so.

**Large files are deliberately not tracked.** `.gitignore` excludes `*.cube`,
`*.dx` and the raw `td.xyz`/`tp.xyz` grids — a full-resolution cube is ~200 MB
and GitHub rejects anything above 100 MB. Regenerate them with `xyzToCube.py`.
The decimated reference grids are the exception and are committed on purpose.

---

## 7. Setting up the git remote

```bash
git init
git add .
git commit -m "Scaffold: folder structure, ignore rules, shared conversion scripts"
git branch -M main
git remote add origin https://github.com/FlxLfr/VMD_esp_visualization.git
git push -u origin main
```

Create the (empty, no README, no .gitignore) repository on GitHub first —
initialising it there and then pushing produces an unrelated-histories conflict.

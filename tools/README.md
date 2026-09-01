# `tools/` — the reference dataset and the σ-hole

Two scripts that sit **outside** the rendering workflow. Neither is called by
`run_allVMD.py`, and nothing in `scripts/` imports them.

| Script | Purpose | Documented in |
|---|---|---|
| `make_reference.py` | builds the committed reference dataset under `reference/` | this file |
| `SigmaHoleCalc.py` | σ-hole per halogen, by ray marching on finished cube files | §4 of this file |

---

## 1  What the reference dataset is for

`reference/brombenzol/` holds a deliberately small bromobenzene dataset that is
committed to the repository. Running the batch driver without arguments

```bash
cd scripts
python run_allVMD.py
```

processes it and writes to separate `_check` names (`images_check/`,
`esp_check.tcl`, `summary_check_*.csv`), so the committed reference files stay
untouched.

The self test answers exactly one question: **does the installation run, and do
the documented numbers come out?** It does not answer "is this an accurate
V_S value" — see §2.

The expected result on this dataset:

```
grid 32x37x24, spacing 0.6000 Bohr, 12 atoms, isovalue 0.001
V_S,min = -0.01863   V_S,max = +0.03070 a.u.  (313 points)  -> +/- 0.0350
```

If those numbers come out, unit conversion, index reordering, the shell band and
the colour-scale rule all work. If they do not, the problem is the installation
or a change to the scripts — not your own data.

---

## 2  Why the grid is 0.60 Bohr, and what that costs

0.60 Bohr is 0.317 Å — five times coarser than the production grids at 0.12 Bohr
(0.063 Å). That is a deliberate choice, and it is measurable:

| | reference, 0.60 Bohr | production, 0.12 Bohr | deviation |
|---|---|---|---|
| shell points | 313 | 38 010 | 0.8 % of the points |
| V_S,min | −0.01863 | −0.01882 | 1.0 % |
| V_S,max | +0.03070 | +0.03154 | 2.7 % |
| σ-hole | +0.01528 (9.59 kcal/(mol·e)) | +0.01629 (10.22) | **6.2 %** |

The V_S values hold to within one to three per cent; the σ-hole falls short by
six. That is the smooth, one-sided degradation the grid study in the background
document describes: the interpolated density smooths the isosurface outwards,
and the potential is less positive there. Sharp features (the σ-hole) suffer,
broad flat ones (V_S,min) barely notice.

**The reason is git, not physics.** The dataset has to be committed so that a
fresh clone can test itself. At 0.60 Bohr with cropping that is 32 × 37 × 24
points, about 1.7 MB per `pointval` file. The same molecule at full resolution
is 15.8 million points and roughly 200 MB per cube — and GitHub rejects anything
above 100 MB. A middle course at 0.30 Bohr would be some 14 MB per file:
possible, but heavy for a repository whose self test is meant to finish in
seconds.

**And the test asks for determinism, not accuracy.** Only the reproducibility of
a documented number matters, so a coarse grid is a feature: it runs in seconds
and the file stays small.

`render_espVMD.py`'s sister script in the PyMOL project says so at run time
whenever the spacing goes past 0.30 Bohr:

```
! grid spacing 0.60 Bohr - too coarse for a trustworthy sigma-hole value;
  expect it to be a few per cent low. Compute finer (a smaller --stride).
```

> **Do not quote a V_S or σ-hole value from the reference dataset.** It exists to
> prove the chain works. Numbers that go into a table come from the full-resolution
> grids under `sandbox/` or `results/`.

One consequence worth knowing: the self test built into `SigmaHoleCalc.py` is
pinned to the coarse value +0.01528 a.u. For its purpose — are the two
implementations still the same code? — that is enough, because both compute the
same six-per-cent-low number. It does not exercise the fine-grid regime.

---

## 3  `make_reference.py`

Builds such a dataset from a full molecule folder.

```bash
python tools/make_reference.py sandbox/brombenzol --name brombenzol
```

It writes a decimated `td.xyz` / `tp.xyz` pair in Turbomole `pointval` format to
`reference/<name>/`, and copies the structure file along — without it
`run_allVMD.py` will not recognise the folder.

### Options

| Option | Default | Effect |
|---|---|---|
| `source` | *required* | molecule folder with `td`/`tp` and a structure file |
| `--name` | the folder name | name under `reference/` |
| `--outdir` | the repository's `reference/` | write somewhere else |
| `--stride N` | `5` | keep every N-th grid point per axis |
| `--margin` | `2.5` | margin around the isosurface, in Bohr |
| `--iso` | `0.001` | the isovalue that must stay fully contained |

### Two steps, both necessary

**Cropping.** The full grid is a 30 Bohr box and the molecule with its
ρ = 0.001 surface fills only the middle of it. The crop is derived from the
**density**: the bounding box of all points with ρ > iso/2, plus `--margin`.
That guarantees the isosurface is fully contained and the images are not cut
off at the edge.

**Decimating.** Then only every `--stride`-th point per axis. Together the two
bring 1.25 GB down to a few MB.

### Why `pointval` and not cube

The self test is meant to exercise the whole chain, `xyzToCubeToVMDVis.py`
included — and unit conversion and index reordering are the two steps most
likely to break. Shipping ready-made cubes would skip exactly those. It also
prefers reading an existing `.cube` if one is present in the source folder,
because that takes seconds where parsing the raw `pointval` files takes minutes;
the output is identical either way.

### If you rebuild the reference

Changing the dataset changes the documented numbers, and three places have to
follow:

- the expected values in §1 above and in the repository `README.md`,
- `REFERENCE_SIGMA` in `SigmaHoleCalc.py`, whose self test is pinned to them,
- the reference numbers quoted in `docs/` and in the comparison document.

The same reference molecule exists in the sister PyMOL project. Keeping the two
identical is what makes the cross-pipeline comparison meaningful, so rebuild
both or neither.

---

## 4  `SigmaHoleCalc.py`

Locates the σ-hole of every halogen from a finished pair of cube files.

```bash
cd tools
python SigmaHoleCalc.py --folder ../results/brombenzol
```

It reads `td.cube` and `tp.cube`. The atom block sits in the cube header, so no
structure file is needed and the alignment question of the converter document
does not arise. Run without `--folder` and it works on `reference/brombenzol`
as a self test (§1).

### Why rays instead of grid points

The σ-hole is a peak *on* the C–X axis. Whether a grid point happens to sit
there *and* inside the thin ρ = iso shell at the same time is luck. So 400 rays
are cast from the halogen into a cone around the axis; on each one the radius at
which ρ crosses the isovalue is located and V read off there, both by trilinear
interpolation. Both numbers are printed, so the difference stays visible:

```
  Cl21
    sigma hole  = +0.01709 a.u. =   +44.9 kJ/(mol*e) =  +10.7 kcal/(mol*e)   [interpolated, 3.8 degrees off the C-Cl axis]
    grid points = +0.02675 a.u. =   +70.2 kJ/(mol*e) =  +16.8 kcal/(mol*e)   [point-based, 117 points in the cap]
    belt        = -0.01836 a.u. =   -48.2 kJ/(mol*e) =  -11.5 kcal/(mol*e)   [508 points]
```

The bracket is the same one `render_esp.py` prints in the PyMOL project, so the
two consoles can be read side by side. The ray count is appended only when rays
were lost — "400 of 400" every time is noise, a shortfall is not.

**Read the angle.** It is the quality control. `0.0 degrees` means the maximum
sits on the axis, the normal case. A value at the rim of the cone (≈36.9° at the
default `--cone`) means there is no maximum inside the cone at all; fluorine
reports that reliably, because it has no σ-hole.

Molecules with several halogens are measured one by one and reported strongest
first. A molecule without a halogen says so and stops.

### Options

| Option | Default | Effect |
|---|---|---|
| `--folder` | `reference/brombenzol` (self test) | molecule folder with `td.cube` and `tp.cube` |
| `--iso` | `0.001` | **isovalue the σ-hole is read on — changes the measured value.** See below. |
| `--rays` | `400` | rays per halogen |
| `--cone` | `0.80` | cosine of the cone half-angle, i.e. 36.9° |
| `--step` | `0.02` | step along a ray in Bohr |
| `--csv [PATH]` | off | also write the result as CSV. Bare `--csv` writes `sigma_holes_<molecule>.csv` next to the cube files. |

### `--iso`: reproducing values computed on another surface

ρ = 0.001 a.u. is the Politzer/Murray convention and the default here, but it is
a convention, not a constant of nature. Part of the literature reports σ-holes on
ρ = 0.002 a.u., and those numbers are not comparable with 0.001 values — the
isovalue is the single parameter that moves the result most. Bromobenzene on the
same 251³ grid:

| `--iso` | σ-hole | |
|---|---|---|
| 0.0005 | +0.0098 a.u. | +6.2 kcal/(mol·e) |
| **0.0010** | **+0.01629** | **+10.2**  ← default, Politzer/Murray |
| 0.0020 | +0.02652 | +16.6 |
| 0.0040 | +0.0434 | +27.2 |

A factor of 4.4 across that range. A larger isovalue puts the surface closer to
the nuclei, where less of the positive nuclear contribution has been screened by
the electrons.

```bash
python SigmaHoleCalc.py --folder ../sandbox/brombenzol --iso 0.002
```

`--iso` exists to reproduce someone else's choice, not to improve on 0.001. Every
value the script prints carries its isovalue in the header line, and the CSV
records it in the `iso_au` column, so a number can always be traced back to the
surface it was measured on.

> The same applies in the PyMOL pipeline, where `--iso` additionally moves the
> surface that is drawn. Here nothing is rendered, so `--iso` affects the
> measurement only.

### On the duplication

The same method exists in the PyMOL project (`render_esp.py`). Two
implementations of one numerical procedure have to agree forever, so the self
test pins them together: run without arguments, `SigmaHoleCalc.py` compares its
result on `reference/brombenzol` against the value the PyMOL pipeline measures on
the same cubes and fails loudly if they have drifted apart. Both give
+0.01528 a.u. there, and they agree to eight decimals on every molecule in
`results/`.

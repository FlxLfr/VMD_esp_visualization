#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SigmaHoleCalc.py - the sigma hole from a pair of cube files
===========================================================

    cd tools
    python SigmaHoleCalc.py --folder ../results/brombenzol

Deliberately outside the rendering workflow. ``run_allVMD.py`` does not call
it, and nothing in the image pipeline depends on it: the scene needs
V_S,min and V_S,max for the colour scale, and those come from grid points near
the shell (``shell_range`` in the converter). Locating the sigma hole is an
analysis of its own, it is slower, and it answers a question the picture does
not ask. Keeping it separate means a broken sigma hole can never cost you an
image set.

It reads the FINISHED cubes, not the pointval files. td.cube carries the atom
block in its header, so no structure file is needed and the question of
whether structure and grid belong together (section 8 of the converter
document) does not arise here - a cube is internally consistent by
construction. A full 251^3 pair is read in about fifteen seconds.


Why the rays, and not the grid points
-------------------------------------
The sigma hole is a peak ON the C-X axis. Whether a grid point happens to sit
there AND inside the thin rho = iso shell at the same time is luck. For
bromobenzene on the 126^3 grid the best point lay 1.14 Bohr off the axis and
the value came out 28 % too low, although 144 points lay inside the cap.

So instead: rays from the halogen into a cone around the axis, on each ray the
radius at which rho crosses the isovalue, and V read off there - both by
trilinear interpolation. The result no longer depends on where the grid points
happen to fall. Both numbers are reported, so the difference stays visible.

Note that trilinear interpolation is the right tool HERE and the wrong one at a
nucleus. On the isosurface the density is smooth and locally almost linear.
At a nucleus it has a cusp, and linear interpolation can only cut a cusp off -
measured at 46 % too low, which is why ``check_alignment`` in the converter
takes the largest value in a 3x3x3 neighbourhood instead. Same grid, two
places, two different correct answers.


On the duplication
------------------
The method below also exists in the sister project (render_esp.py). That is a
liability, not a feature: two implementations of one numerical procedure have
to agree forever, and the kind of bug fixed here in August - the ray that ran
into a different part of the molecule - is exactly what gets repaired in one
copy and silently kept in the other. A cross-repository import was not an
option, because both repositories have to stay clonable on their own.

The answer is the self test. Run without ``--folder`` this script works on
reference/brombenzol and compares its result against the value the PyMOL
pipeline measures on the same cubes (REFERENCE_SIGMA below). If the two
implementations drift apart, the run says so instead of quietly reporting a
different number.

And to be clear about what this is not: section 1 of Details.docx argues that
two independent viewers are a control on each other. This script is a PORT,
not an independent derivation, so it does not carry that weight. It makes the
sigma hole available where until now only rendering happened.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "scripts"))

import xyzToCubeToVMDVis as conv                                # noqa: E402
from constants import (BOHR_PER_ANGSTROM, HARTREE_TO_KJ,        # noqa: E402
                       HARTREE_TO_KCAL)

# van der Waals radii in Angstrom, halogens only - the distance limit below is
# derived from them and nothing else needs them.
VDW_ANGSTROM = {9: 1.47, 17: 1.75, 35: 1.85, 53: 1.98}

# The value the PyMOL pipeline measures on reference/brombenzol, and the
# tolerance the self test allows. Filled in by measurement, not by hand - see
# the module docstring.
REFERENCE_SIGMA = 0.01527802
REFERENCE_TOL = 5.0e-5


# ----------------------------------------------------------------------------
# The shell, as grid points
# ----------------------------------------------------------------------------

def shell_points(density, esp, origin, voxel, iso=0.001):
    """Coordinates and ESP values of the grid points on the rho = iso shell.

    The same band the converter uses for the colour scale, and the same
    widening when it turns out too thin - so the point-based comparison value
    below is measured on exactly the points the rest of the pipeline sees.

    Only the shell points are materialised. A full coordinate array of a 251^3
    grid would be several hundred MB.
    """
    mask = np.abs(density - iso) < iso * conv.SHELL_TOL
    if mask.sum() < conv.SHELL_MIN_POINTS:
        mask = np.abs(density - iso) < iso * conv.SHELL_TOL_WIDE
    idx = np.argwhere(mask)                        # (N, 3) grid indices
    return origin + idx @ voxel, esp[mask]         # Bohr, a.u.


# ----------------------------------------------------------------------------
# Which halogens, and where do they point
# ----------------------------------------------------------------------------

def halogen_axes(atoms):
    """Every halogen with its C-X axis.

    A molecule can carry more than one, and each is evaluated separately.

    Returns dicts with
      index    0-based atom index
      symbol   element symbol
      label    symbol plus 1-based number, e.g. "Cl21", as in the output
      pos      coordinates of the halogen, Bohr
      axis     normalised C->X direction; it points at the sigma hole
      r_limit  the largest distance at which the halogen's OWN surface can lie

    Halogens with no carbon in reach are skipped.

    Why r_limit exists: in a folded molecule the cone around the C-X axis does
    not point into empty space but at another part of the molecule, and that
    part's surface gets measured as well. Triazolam is such a case - the cone
    around Cl21 hits the methyl group on the triazole ring, which produced a
    "sigma hole" of +18.8 instead of +10.7 kcal/(mol*e). A rho = 0.001 surface
    sits at roughly 1.1 to 1.2 van der Waals radii; the factor 1.6 leaves room
    and excludes everything beyond.
    """
    coords = np.array([[a[1], a[2], a[3]] for a in atoms], dtype=float)
    znums = np.array([a[0] for a in atoms], dtype=int)
    carbons = [i for i, z in enumerate(znums) if z == 6]
    out = []
    if not carbons:
        return out
    for hi, z in enumerate(znums):
        if int(z) not in conv.HALOGENS:
            continue
        d = np.linalg.norm(coords[carbons] - coords[hi], axis=1)
        ci = carbons[int(np.argmin(d))]
        axis = coords[hi] - coords[ci]
        n = np.linalg.norm(axis)
        if n < 1e-6:
            continue
        out.append({"index": hi,
                    "symbol": conv.HALOGENS[int(z)],
                    "label": f"{conv.HALOGENS[int(z)]}{hi + 1}",
                    "pos": coords[hi],
                    "axis": axis / n,
                    "r_limit": 1.6 * VDW_ANGSTROM[int(z)] * BOHR_PER_ANGSTROM})
    return out


# ----------------------------------------------------------------------------
# Interpolation and ray directions
# ----------------------------------------------------------------------------

def trilinear(vol, origin, delta, pts):
    """Trilinear interpolation on an axis-aligned, regular grid."""
    f = (pts - origin) / delta
    i0 = np.floor(f).astype(int)
    n = np.array(vol.shape)
    i0 = np.clip(i0, 0, n - 2)
    frac = np.clip(f - i0, 0.0, 1.0)

    out = np.zeros(len(pts))
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                w = ((frac[:, 0] if dx else 1 - frac[:, 0]) *
                     (frac[:, 1] if dy else 1 - frac[:, 1]) *
                     (frac[:, 2] if dz else 1 - frac[:, 2]))
                out += w * vol[i0[:, 0] + dx, i0[:, 1] + dy, i0[:, 2] + dz]
    return out


def cone_directions(axis, cone_cos, n=400):
    """Directions in a cap around ``axis``; the first one IS the axis.

    The rest is a Fibonacci spiral over the spherical cap - even coverage,
    without the crowding near the axis that spherical coordinates produce.

    Why the axis is prepended instead of falling out of the spiral: the + 0.5
    in k is the midpoint rule for an EQUAL-AREA distribution, every ray sitting
    in the middle of a ring of the same size. That is correct when averaging
    over the cap. But this is a search for a MAXIMUM, and for the sigma hole
    that maximum sits exactly on the axis. With the offset alone the innermost
    ray was 1.281 degrees away, the axis itself was never evaluated, and every
    axially symmetric molecule stubbornly reported "1.3 degrees" - a lower
    bound of the sampling grid, not a measurement. The error in the value was
    small (4-bromoacetophenone: 0.008 kcal/(mol*e)); the reported angle was
    misleading, and the angle is the quality control.
    """
    axis = np.asarray(axis, dtype=float)

    m = max(1, n - 1)                              # n-1 spiral directions
    k = np.arange(m) + 0.5
    cosv = 1.0 - (1.0 - cone_cos) * k / m          # cone_cos .. 1
    phi = np.pi * (1 + 5 ** 0.5) * k               # golden angle
    sinv = np.sqrt(np.maximum(0.0, 1 - cosv ** 2))

    tmp = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(tmp, axis)) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    e1 = np.cross(axis, tmp)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)

    spiral = (cosv[:, None] * axis
              + (sinv * np.cos(phi))[:, None] * e1
              + (sinv * np.sin(phi))[:, None] * e2)
    return np.vstack([axis[None, :], spiral])


# ----------------------------------------------------------------------------
# The measurement
# ----------------------------------------------------------------------------

def sigma_hole_rays(density, esp, origin, voxel, hal, iso=0.001,
                    cone_cos=0.80, n_rays=400, dr=0.02, r_max=14.0):
    """The sigma hole of ONE halogen, by ray marching.

    Returns (value in a.u., angle to the axis in degrees, number of rays that
    actually found the surface) or (None, None, 0).
    """
    diag = np.diag(voxel)
    if not np.allclose(voxel, np.diag(diag)):
        return None, None, 0                  # skewed grid - no cheap lookup
    delta = diag

    dirs = cone_directions(hal["axis"], cone_cos, n_rays)
    radii = np.arange(1.0, min(r_max, hal["r_limit"]), dr)

    best_v, best_cos, hits = None, None, 0
    for d in dirs:
        pts = hal["pos"] + radii[:, None] * d[None, :]
        rho = trilinear(density, origin, delta, pts)
        # The INNERMOST crossing, going outwards. The ray starts deep inside
        # the halogen's own density, so the first drop below the isovalue is
        # its own surface. Taking the outermost crossing - the original
        # version - lets the ray dive through and measure whatever lies behind
        # it in a folded molecule.
        if rho[0] < iso:
            continue                          # ray already starts outside
        below = np.nonzero(rho < iso)[0]
        if below.size == 0:
            continue                          # surface not hit within r_limit
        j = below[0] - 1
        if j < 0 or j + 1 >= len(radii):
            continue
        r0, r1 = radii[j], radii[j + 1]
        y0, y1 = rho[j], rho[j + 1]
        rs = r0 + (iso - y0) * (r1 - r0) / (y1 - y0) if y1 != y0 else r0
        v = float(trilinear(esp, origin, delta,
                            (hal["pos"] + rs * d)[None, :])[0])
        hits += 1
        if best_v is None or v > best_v:
            best_v, best_cos = v, float(np.dot(d, hal["axis"]))

    if best_v is None:
        return None, None, 0
    return best_v, float(np.degrees(np.arccos(min(1.0, best_cos)))), hits


def cap_and_belt(pos, vals, hal, cone_cos=0.80, belt_cos=0.35,
                 belt_factor=1.5):
    """The point-based cap maximum and the belt minimum for one halogen.

    The cap value is the comparison number for the ray result. The belt is the
    negative ring perpendicular to the C-X axis - the other half of the
    anisotropy that makes a halogen bond directional, and it has no ray-based
    counterpart because it is not a peak on a known axis.

    The distance limit is the same one as in halogen_axes: without it, points
    on a completely different part of a folded molecule count towards the cap.
    """
    out = {}
    rel = pos - hal["pos"]
    r = np.linalg.norm(rel, axis=1)
    r[r == 0] = 1e-9
    cos = (rel @ hal["axis"]) / r

    cap = (cos > cone_cos) & (r < hal["r_limit"])
    if cap.sum() >= 5:
        out["cap_max"] = float(vals[cap].max())
        out["cap_points"] = int(cap.sum())
        r_cap = float(r[cap].mean())
        belt = (np.abs(cos) < belt_cos) & (r < belt_factor * r_cap)
        if belt.sum() >= 5:
            out["belt_min"] = float(vals[belt].min())
            out["belt_points"] = int(belt.sum())
    return out


def analyse(folder, iso=0.001, cone_cos=0.80, n_rays=400, dr=0.02,
            verbose=True):
    """Every halogen of one molecule folder, strongest sigma hole first."""
    td = os.path.join(folder, "td.cube")
    tp = os.path.join(folder, "tp.cube")
    for path in (td, tp):
        if not os.path.exists(path):
            raise SystemExit(
                f"{folder}: {os.path.basename(path)} is missing. This script "
                f"works on the finished cubes - run run_allVMD.py on the "
                f"folder first.")

    if verbose:
        print(f"reading {td} ...", flush=True)
    density = conv.read_cube(td)
    atoms, origin, voxel = conv.read_cube_frame(td)
    if verbose:
        print(f"reading {tp} ...", flush=True)
    esp = conv.read_cube(tp)
    if density.shape != esp.shape:
        raise SystemExit("Density and ESP cube are on different grids.")

    spacing = float(np.max(np.abs(np.diag(voxel))))
    if verbose:
        print(f"grid {'x'.join(map(str, density.shape))}, {len(atoms)} atoms, "
              f"spacing {spacing:.4f} Bohr, isovalue {iso:g}\n", flush=True)

    pos, vals = shell_points(density, esp, origin, voxel, iso=iso)

    rows = []
    for hal in halogen_axes(atoms):
        v, angle, hits = sigma_hole_rays(density, esp, origin, voxel, hal,
                                         iso=iso, cone_cos=cone_cos,
                                         n_rays=n_rays, dr=dr)
        row = {"label": hal["label"], "symbol": hal["symbol"],
               "index": hal["index"], "sigma_max": v, "sigma_angle": angle,
               "rays_hit": hits, "n_rays": n_rays, "spacing_bohr": spacing,
               "grid": "x".join(map(str, density.shape)), "iso": iso}
        row.update(cap_and_belt(pos, vals, hal, cone_cos=cone_cos))
        rows.append(row)

    # Strongest first: with several halogens that is the one the sigma view of
    # the scene is oriented on, so the two outputs agree on what "the" sigma
    # hole of the molecule is.
    rows.sort(key=lambda e: -(e["sigma_max"] if e["sigma_max"] is not None
                              else -math.inf))
    return rows


# ----------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------

def _three_units(v):
    return (f"{v:+.5f} a.u. = {v * HARTREE_TO_KJ:+7.1f} kJ/(mol*e) "
            f"= {v * HARTREE_TO_KCAL:+6.1f} kcal/(mol*e)")


def report(rows, folder):
    name = os.path.basename(os.path.normpath(folder))
    print("-" * 70)
    if not rows:
        print(f"{name}: no halogen with a bonded carbon - nothing to measure.")
        print("-" * 70)
        return
    for r in rows:
        print(f"  {r['label']}")
        if r["sigma_max"] is None:
            print("    sigma hole  : not evaluable - no ray reached the "
                  "isosurface within the distance limit")
        else:
            # Same bracket as render_esp.py in the sister project, so the two
            # consoles can be read side by side: how the value was obtained,
            # then the angle.
            #
            # The angle is the quality control. 0 degrees means the maximum
            # sits on the axis, which is the normal case. A value at the rim of
            # the cone (about 36.9 degrees at the default opening) means there
            # is no maximum inside the cone at all - fluorine reports that
            # reliably, because it has no sigma hole.
            #
            # The ray count is appended only when rays were LOST. "400 of 400"
            # every time is noise; a shortfall means part of the cone found no
            # surface within the distance limit, and that is worth seeing.
            note = (f"interpolated, {r['sigma_angle']:.1f} degrees off the "
                    f"C-{r['symbol']} axis")
            if r["rays_hit"] < r["n_rays"]:
                note += f", only {r['rays_hit']} of {r['n_rays']} rays"
            print(f"    sigma hole  = {_three_units(r['sigma_max'])}"
                  f"   [{note}]")
        if "cap_max" in r:
            print(f"    grid points = {_three_units(r['cap_max'])}"
                  f"   [point-based, {r['cap_points']} points in the cap]")
        if "belt_min" in r:
            print(f"    belt        = {_three_units(r['belt_min'])}"
                  f"   [{r['belt_points']} points]")
        print()
    print("-" * 70)
    if any(r["sigma_max"] is not None and "cap_max" in r for r in rows):
        print("The 'grid points' line is the point-based value on the same "
              "shell.\nIt depends on where the grid points happen to fall; "
              "the ray value does not.")
        print("-" * 70)


FIELDS = ["molecule", "halogen", "atom_index", "grid", "spacing_bohr",
          "iso_au", "sigma_hole_au", "sigma_hole_kJ", "sigma_hole_kcal",
          "sigma_angle_deg", "sigma_method", "rays_hit", "n_rays",
          "pointbased_au", "cap_points", "belt_min_au", "belt_points"]


def write_csv(path, rows, folder):
    name = os.path.basename(os.path.normpath(folder))
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for r in rows:
            s = r["sigma_max"]
            w.writerow({
                "molecule": name,
                "halogen": r["label"],
                "atom_index": r["index"] + 1,
                "grid": r["grid"],
                "spacing_bohr": f"{r['spacing_bohr']:.4f}",
                "iso_au": f"{r['iso']:g}",
                "sigma_hole_au": "" if s is None else f"{s:.5f}",
                "sigma_hole_kJ": "" if s is None else f"{s * HARTREE_TO_KJ:.1f}",
                "sigma_hole_kcal": ("" if s is None
                                    else f"{s * HARTREE_TO_KCAL:.2f}"),
                "sigma_angle_deg": ("" if r["sigma_angle"] is None
                                    else f"{r['sigma_angle']:.2f}"),
                # Constant by construction: this script has no point-based
                # fallback for the sigma hole - if no ray finds the surface the
                # value stays empty. The column exists so the table can be laid
                # beside the PyMOL summary, which does carry both methods.
                "sigma_method": "" if s is None else "interpolated",
                "rays_hit": r["rays_hit"],
                "n_rays": r["n_rays"],
                "pointbased_au": ("" if "cap_max" not in r
                                  else f"{r['cap_max']:.5f}"),
                "cap_points": r.get("cap_points", ""),
                "belt_min_au": ("" if "belt_min" not in r
                                else f"{r['belt_min']:.5f}"),
                "belt_points": r.get("belt_points", ""),
            })
    print(f"CSV: {path}")


# ----------------------------------------------------------------------------
# Main program
# ----------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        description="Determine the sigma hole from td.cube and tp.cube. Runs "
                    "outside the rendering workflow; run_allVMD.py does not "
                    "call it.")
    p.add_argument("--folder", default=None,
                   help="molecule folder with td.cube and tp.cube. Without it "
                        "the self test runs on reference/brombenzol and "
                        "compares against the value of the PyMOL pipeline.")
    p.add_argument("--iso", type=float, default=0.001,
                   help="isovalue of the surface the sigma hole is read on "
                        "(default: 0.001). Changes the measured value - see "
                        "section 6.2 of the converter document.")
    p.add_argument("--rays", type=int, default=400,
                   help="rays per halogen (default: 400)")
    p.add_argument("--cone", type=float, default=0.80,
                   help="cosine of the cone half angle (default: 0.80, "
                        "i.e. 36.9 degrees)")
    p.add_argument("--step", type=float, default=0.02,
                   help="step along a ray in Bohr (default: 0.02)")
    # nargs="?" so that a bare --csv works: the rest of the project names its
    # output files by itself (summary_<time>_<date>.csv, stride_sweep_<name>.csv),
    # and being made to invent a file name for a one-line table is friction.
    p.add_argument("--csv", nargs="?", const=True, default=None,
                   metavar="PATH",
                   help="also write the result as CSV. Without a path: "
                        "sigma_holes_<molecule>.csv next to the cube files.")
    args = p.parse_args(argv)

    self_test = args.folder is None
    folder = os.path.abspath(
        args.folder or os.path.join(_HERE, "..", "reference", "brombenzol"))

    print("=" * 70)
    print("SigmaHoleCalc.py - sigma hole from the cube files")
    print("=" * 70)
    if self_test:
        print("Self test on the reference dataset. The value is compared "
              "against the\none the PyMOL pipeline measures on the same "
              "cubes; a difference means\nthe two implementations have "
              "drifted apart.")
    print(f"folder: {folder}\n")

    rows = analyse(folder, iso=args.iso, cone_cos=args.cone,
                   n_rays=args.rays, dr=args.step)
    report(rows, folder)

    if args.csv:
        out = args.csv
        if out is True:
            name = os.path.basename(os.path.normpath(folder))
            out = os.path.join(folder, f"sigma_holes_{name}.csv")
        write_csv(out, rows, folder)

    if not self_test:
        return 0

    # --- the cross-check -------------------------------------------------
    measured = rows[0]["sigma_max"] if rows else None
    if measured is None:
        print("! self test: no sigma hole measured on the reference dataset.")
        return 1
    diff = abs(measured - REFERENCE_SIGMA)
    print(f"self test: {measured:+.5f} a.u. against {REFERENCE_SIGMA:+.5f} "
          f"from the PyMOL pipeline, difference {diff:.2e}")
    if diff > REFERENCE_TOL:
        print(f"! The two implementations disagree by more than "
              f"{REFERENCE_TOL:g} a.u.")
        print("! One of them has been changed without the other. Compare with "
              "sigma_hole_interpolated in the PyMOL project's render_esp.py.")
        return 1
    print("ok - both implementations agree.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

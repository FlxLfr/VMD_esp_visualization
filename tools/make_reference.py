#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_reference.py - a small reference dataset from a full molecule folder

    python tools/make_reference.py sandbox/brombenzol --name brombenzol

Writes a decimated td.xyz / tp.xyz in Turbomole pointval format plus the
structure file to reference/<name>/. That is the dataset for the self test:
`python run_allVMD.py` without arguments runs on it.

Why pointval and not cube
-------------------------
The self test is meant to check the whole chain, xyzToCubeToVMDVis.py
included - unit conversion and index reordering are the places most likely to
break. If it found ready-made cubes, it would skip exactly that step.

What is decimated
-----------------
Two steps, both necessary:

  1. **Cropping.** The full grid is a 30 Bohr box, and the molecule with its
     rho = 0.001 surface fills only the middle of it. The crop is determined
     from the DENSITY: the bounding box of all points with rho > iso/2, plus a
     margin. That guarantees the isosurface is fully contained and the images
     do not look cut off.
  2. **Decimating.** Only every n-th point per axis.

Together the two bring 1.25 GB down to a few MB. The dataset is deliberately
too coarse for a citable V_S value - it answers the question "does the
installation run and do the documented numbers come out", not "how large is
the sigma hole".

The same parameters as in the sister project
--------------------------------------------
The PyMOL pipeline builds its reference set with the same settings
(brombenzol, --stride 5, --margin 2.5). Only that way do both repositories
hold numbers that may be laid side by side: if the grid resolution differs,
V_S,min and V_S,max differ by one to three per cent, and the comparison then
measures the decimation instead of the viewers.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "scripts"))
sys.path.insert(0, _HERE)

GRIDS = ("td", "tp")
QUANTITY = {"td": "density", "tp": "electrostatic potential"}


def read_cube(path):
    """Gaussian cube -> (values[i1,i2,i3], origin, delta) in Bohr."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        fh.readline()
        fh.readline()
        parts = fh.readline().split()
        natoms, origin = abs(int(parts[0])), [float(v) for v in parts[1:4]]
        dims, delta = [], []
        for i in range(3):
            p = fh.readline().split()
            dims.append(int(p[0]))
            delta.append(float(p[1 + i]))     # diagonal element
        for _ in range(natoms):
            fh.readline()
        # np.fromstring parses in C and builds no intermediate object in
        # Python. A full grid has 15.8 million values; via .split() that would
        # be a list of 15.8 million str objects, a good gigabyte for the
        # intermediate step alone.
        values = np.fromstring(fh.read(), dtype=np.float32, sep=" ")
    if values.size != dims[0] * dims[1] * dims[2]:
        raise ValueError(f"{path}: {values.size} values, expected "
                         f"{dims[0] * dims[1] * dims[2]}")
    return values.reshape(dims), np.array(origin), np.array(delta)


def load(folder, tag, verbose=True):
    """Prefer the cube (seconds), otherwise the pointval file (minutes)."""
    cube = os.path.join(folder, f"{tag}.cube")
    if os.path.exists(cube):
        if verbose:
            print(f"    reading {tag}.cube")
        return read_cube(cube)
    raw = os.path.join(folder, f"{tag}.xyz")
    if not os.path.exists(raw):
        raise SystemExit(f"{folder}: neither {tag}.cube nor {tag}.xyz")
    if verbose:
        print(f"    reading {tag}.xyz (pointval, this takes a while)")
    import xyzToCubeToVMDVis as conv
    info, data = conv.read_values(raw, verbose=verbose)
    origin = np.array([info["grid"][i][0] for i in range(3)])
    delta = np.array([info["grid"][i][1] for i in range(3)])
    return data, origin, delta


def write_pointval(path, data, origin, delta, quantity, title="101"):
    """Write Turbomole pointval: x varies fastest."""
    n1, n2, n3 = data.shape
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        # #origin stays at zero, the offset sits in the #grid lines.
        # Turbomole writes it the same way, and xyzToCubeToVMDVis.py adds both
        # of them: were the offset in both places, the grid would sit twice as
        # far from the origin as the atoms - the classic "molecule floats next
        # to its surface" error.
        fh.write(f"#origin      {0.0:14.6f}{0.0:14.6f}{0.0:14.6f}\n")
        for i in range(3):
            v = [0.0, 0.0, 0.0]
            v[i] = 1.0
            fh.write(f"#vector{i+1}     {v[0]:14.6f}{v[1]:14.6f}{v[2]:14.6f}\n")
        for i, n in enumerate((n1, n2, n3)):
            fh.write(f"#grid{i+1}  start  {origin[i]:.6f}  delta  "
                     f"{delta[i]:.6f}  points  {n}\n")
        fh.write(f"#title for this grid {title}\n")
        fh.write(f"#{quantity}\n")
        fh.write("#plotdata\n")
        fh.write("# cartesian coordinates x,y,z and f(x,y,z)\n")

        xs = origin[0] + delta[0] * np.arange(n1)
        ys = origin[1] + delta[1] * np.arange(n2)
        zs = origin[2] + delta[2] * np.arange(n3)
        buf = []
        for k in range(n3):                      # z outermost ...
            for j in range(n2):                  # ... x innermost
                col = data[:, j, k]
                for i in range(n1):
                    buf.append(f"   {xs[i]:14.8f} {ys[j]:14.8f} {zs[k]:14.8f} "
                               f"{col[i]:14.8f}\n")
                if len(buf) >= 65536:
                    fh.write("".join(buf))
                    buf.clear()
        fh.write("".join(buf))


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Build a decimated reference dataset for the self test.")
    p.add_argument("source", help="molecule folder with td/tp and a structure file")
    p.add_argument("--name", help="name under reference/ (default: folder name)")
    p.add_argument("--outdir", default=None,
                   help="default: the repository's reference/")
    p.add_argument("--stride", type=int, default=5,
                   help="every n-th grid point per axis (default 5)")
    p.add_argument("--margin", type=float, default=2.5,
                   help="margin around the isosurface in Bohr (default 2.5)")
    p.add_argument("--iso", type=float, default=0.001,
                   help="isovalue that must be fully contained")
    args = p.parse_args(argv)

    src = os.path.normpath(args.source)
    name = args.name or os.path.basename(src)
    outdir = args.outdir or os.path.normpath(
        os.path.join(_HERE, "..", "reference"))
    dest = os.path.join(outdir, name)
    os.makedirs(dest, exist_ok=True)

    print("=" * 66)
    print(f"make_reference.py - {src}  ->  {dest}")
    print("=" * 66)

    dens, origin, delta = load(src, "td")
    esp, o2, d2 = load(src, "tp")
    if dens.shape != esp.shape or not np.allclose(origin, o2):
        raise SystemExit("td and tp are not on the same grid.")
    print(f"[1] grid {dens.shape[0]}x{dens.shape[1]}x{dens.shape[2]}, "
          f"delta {delta[0]:.4f} Bohr, origin {origin[0]:.2f}")

    # Crop from the density: everything the isosurface needs.
    mask = dens > args.iso / 2.0
    if not mask.any():
        raise SystemExit(f"No points with rho > {args.iso/2:g} found.")
    lo, hi = [], []
    for ax in range(3):
        idx = np.nonzero(mask.any(axis=tuple(a for a in range(3) if a != ax)))[0]
        pad = int(np.ceil(args.margin / delta[ax]))
        lo.append(max(0, idx[0] - pad))
        hi.append(min(dens.shape[ax] - 1, idx[-1] + pad))
    sl = tuple(slice(lo[a], hi[a] + 1, args.stride) for a in range(3))
    sub_d, sub_e = dens[sl], esp[sl]
    new_origin = origin + delta * np.array(lo)
    new_delta = delta * args.stride
    print(f"[2] crop {sub_d.shape[0]}x{sub_d.shape[1]}x{sub_d.shape[2]}, "
          f"delta {new_delta[0]:.4f} Bohr "
          f"({sub_d.size:,} instead of {dens.size:,} points)")
    print(f"    rho {sub_d.min():.2e} .. {sub_d.max():.4g} | "
          f"ESP {sub_e.min():+.4g} .. {sub_e.max():+.4g}")

    for tag, arr in (("td", sub_d), ("tp", sub_e)):
        out = os.path.join(dest, f"{tag}.xyz")
        write_pointval(out, arr, new_origin, new_delta, QUANTITY[tag])
        print(f"[3] -> {out}  ({os.path.getsize(out)/1024**2:.1f} MB)")

    # Take the structure file along - without it run_allVMD.py does not find
    # the folder.
    for fn in sorted(os.listdir(src)):
        stem, ext = os.path.splitext(fn)
        if ext.lower() in (".mol", ".sdf", ".xyz") and stem not in GRIDS:
            shutil.copy(os.path.join(src, fn), os.path.join(dest, fn))
            print(f"[4] -> {os.path.join(dest, fn)}")
            break
    else:
        print("    ! no structure file found - add one by hand",
              file=sys.stderr)

    print("\nSelf test:  python run_allVMD.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())

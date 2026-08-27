#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
xyzToCubeToVMDVis.py - Turbomole pointval -> Gaussian cube -> VMD scene

    python xyzToCubeToVMDVis.py --struct brombenzol_aro_opti.mol td.xyz tp.xyz
    vmd -e esp.tcl

Writes td.cube, tp.cube and an esp.tcl that loads both cubes, builds the
isosurface of the electron density and colours it by the electrostatic
potential. The scene itself sits next to this file in esp_template.tcl.

Scope: for the colour scale, V_S,min and V_S,max are determined on the
rho = iso shell - without the scale there is no meaningful figure. The grid
points in the band rho = iso +/- 12 % are used for that, with the same
constants as xyzToCube.py in the sister project Pymol_esp_visualization. It is
the same calculation, and the numbers agree accordingly - which is exactly
what the self test checks.

Deliberately NOT included is the sigma hole. Finding it requires trilinear
interpolation along the C-X axis, and that code lives over there, parameter
studied and documented; a second, independently maintained copy of it would be
precisely the kind of duplicate that drifts apart unnoticed. So for sigma hole
values look into the PyMOL pipeline. V_S,min and V_S,max are the same here.

Two pitfalls in the data format that the script catches:
  1. Axis order - Turbomole varies x fastest, cube varies z fastest.
  2. Units - the grid is in Bohr, the structure file usually in Angstrom.

The structure file only provides the atom block of the cube header (the
pointval file contains no atoms). After that VMD does not need it any more.

Only numpy is required.
"""

from __future__ import annotations

import argparse
import datetime
import math
import os
import re
import sys
import time

import numpy as np

from constants import BOHR_PER_ANGSTROM, HARTREE_TO_KJ

ELEMENTS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co "
    "Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb "
    "Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re "
    "Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es "
    "Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()
SYMBOL_TO_Z = {s.upper(): i + 1 for i, s in enumerate(ELEMENTS)}

CHUNK_BYTES = 1 << 24          # 16 MB read buffer
TEMPLATE = "esp_template.tcl"  # next to this script

# Shell thickness relative to the isovalue. Thin enough that the values really
# come from the isosurface; if that leaves it too sparsely populated, widen it
# once.
SHELL_TOL, SHELL_TOL_WIDE, SHELL_MIN_POINTS = 0.12, 0.30, 50


# ----------------------------------------------------------------------------
# Structure file -> atoms in Bohr
# ----------------------------------------------------------------------------

def _z(sym, path):
    key = sym.strip().capitalize().upper()
    if key in SYMBOL_TO_Z:
        return SYMBOL_TO_Z[key]
    if re.fullmatch(r"\d+", sym.strip()):
        return int(sym)
    raise ValueError(f"Unknown element '{sym}' in {path}.")


def _read_xyz(lines, path):
    """``symbol x y z``, with or without the two header lines.

    The Turbomole structure files here have none - which is also why they
    cannot be loaded into VMD directly ("Unable to load molecule"). Taking the
    detour via the cube makes that irrelevant, the atoms end up in the cube
    header.
    """
    start = 2 if lines and re.fullmatch(r"\d+", lines[0].strip()) else 0
    atoms = []
    for line in lines[start:]:
        p = line.split()
        if len(p) < 4:
            continue
        try:
            atoms.append((_z(p[0], path), float(p[1]), float(p[2]), float(p[3])))
        except ValueError:
            continue
    return atoms


def _read_molfile(lines, path):
    """MDL molfile V2000/V3000. Coordinates come BEFORE the symbol, in Angstrom."""
    if any("V3000" in ln for ln in lines[:8]):
        atoms, inside = [], False
        for line in lines:
            s = line.strip()
            if s.startswith("M  V30 BEGIN ATOM"):
                inside = True
            elif s.startswith("M  V30 END ATOM"):
                break
            elif inside and s.startswith("M  V30"):
                p = s.split()
                if len(p) >= 7:
                    atoms.append((_z(p[3], path),
                                  float(p[4]), float(p[5]), float(p[6])))
        return atoms

    if len(lines) < 5:
        raise ValueError(f"{path}: too short for a molfile.")
    try:
        natoms = int(lines[3][0:3])          # counts line
    except ValueError:
        raise ValueError(f"{path}: counts line not readable: {lines[3].strip()!r}")

    atoms = []
    for line in lines[4:4 + natoms]:
        p = line.split()
        if len(p) < 4:
            raise ValueError(f"{path}: atom line incomplete: {line.strip()!r}")
        atoms.append((_z(p[3], path), float(p[0]), float(p[1]), float(p[2])))
    return atoms


def read_structure(path, unit="angstrom"):
    """.xyz / .mol / .sdf -> a list of (Z, x, y, z) in Bohr."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    if not lines:
        raise ValueError(f"{path} is empty.")

    if os.path.splitext(path)[1].lower() in (".mol", ".sdf", ".sd"):
        for i, ln in enumerate(lines):
            if ln.startswith("$$$$"):        # SD file: the first record only
                lines = lines[:i]
                break
        atoms, unit = _read_molfile(lines, path), "angstrom"
    else:
        atoms = _read_xyz(lines, path)

    if not atoms:
        raise ValueError(f"No atoms found in {path}.")
    if unit == "angstrom":
        b = BOHR_PER_ANGSTROM
        atoms = [(z, x * b, y * b, zz * b) for (z, x, y, zz) in atoms]
    elif unit != "bohr":
        raise ValueError("unit must be 'angstrom' or 'bohr'")
    return atoms


# ----------------------------------------------------------------------------
# Turbomole grid file
# ----------------------------------------------------------------------------

def parse_header(fh):
    """Reads the ``#`` header. Returns (info, the first data line already read)."""
    info = {"origin": np.zeros(3), "vectors": np.eye(3),
            "grid": [None] * 3, "title": "", "quantity": ""}
    while True:
        line = fh.readline()
        if not line:
            raise ValueError("File ends inside the header - no data found.")
        if not line.startswith("#"):
            break
        body = line[1:].strip()
        low = body.lower()
        if low.startswith("origin"):
            info["origin"] = np.array([float(v) for v in body.split()[1:4]])
        elif low.startswith("vector"):
            info["vectors"][int(body[6]) - 1] = [float(v) for v in body.split()[1:4]]
        elif low.startswith("grid"):
            m = re.search(r"start\s+(\S+)\s+delta\s+(\S+)\s+points\s+(\d+)",
                          body, re.I)
            if not m:
                raise ValueError(f"Grid line not readable: {line!r}")
            info["grid"][int(body[4]) - 1] = (float(m.group(1)),
                                              float(m.group(2)), int(m.group(3)))
        elif low.startswith("title"):
            info["title"] = body
        elif ("potential" in low or "density" in low) \
                and not low.startswith("cartesian"):
            info["quantity"] = body

    if any(g is None for g in info["grid"]):
        raise ValueError("Header incomplete: #grid1/#grid2/#grid3 are missing.")
    return info, line


def read_values(path, verbose=True):
    """The 4th column of a pointval file as a float32 grid [i1, i2, i3].

    Block by block instead of line by line - for the 1.25 GB files that is
    about an order of magnitude faster.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        info, remainder = parse_header(fh)
        n1, n2, n3 = (info["grid"][i][2] for i in range(3))
        total = n1 * n2 * n3
        values = np.empty(total, dtype=np.float32)
        filled, t0 = 0, time.time()

        while True:
            chunk = fh.read(CHUNK_BYTES)
            if not chunk:
                break
            chunk = remainder + chunk
            cut = chunk.rfind("\n")
            if cut == -1:
                remainder = chunk
                continue
            remainder, block = chunk[cut + 1:], chunk[:cut]

            tokens = block.split()
            if not tokens:
                continue
            if len(tokens) % 4:
                raise ValueError(f"{path}: expected 4 columns per line.")
            arr = np.asarray(tokens, dtype=np.float32).reshape(-1, 4)[:, 3]
            if filled + arr.size > total:
                raise ValueError(f"{path}: more data points than in the header.")
            values[filled:filled + arr.size] = arr
            filled += arr.size
            if verbose:
                sys.stdout.write(f"\r    reading {os.path.basename(path)}: "
                                 f"{100.0 * filled / total:5.1f} %")
                sys.stdout.flush()

        tokens = remainder.split()
        if tokens:
            if len(tokens) % 4:
                raise ValueError(f"{path}: last line incomplete.")
            arr = np.asarray(tokens, dtype=np.float32).reshape(-1, 4)[:, 3]
            values[filled:filled + arr.size] = arr
            filled += arr.size

    if verbose:
        sys.stdout.write(f"\r    reading {os.path.basename(path)}: 100.0 %  "
                         f"({filled:,} points in {time.time() - t0:.1f} s)\n")
    if filled != total:
        raise ValueError(f"{path}: read {filled} values, expected {total}.")

    # Turbomole: x fastest -> [i3, i2, i1]. Cube: z fastest.
    return info, np.ascontiguousarray(
        np.transpose(values.reshape(n3, n2, n1), (2, 1, 0)))


# ----------------------------------------------------------------------------
# Colour scale from the isosurface
# ----------------------------------------------------------------------------

def shell_range(density, esp, iso=0.001, step=0.005):
    """V_S,min / V_S,max on the rho = iso shell and the symmetric scale.

    Do NOT take the range from the whole grid - there the nuclear
    singularities dominate with several hundred a.u. Only grid points near the
    isosurface count, and the result is rounded up symmetrically to a multiple
    of ``step``.

    Returns (vmin, vmax, half_width, number_of_points) or None.
    """
    if density is None or esp is None or density.shape != esp.shape:
        return None
    mask = np.abs(density - iso) < iso * SHELL_TOL
    if mask.sum() < SHELL_MIN_POINTS:
        mask = np.abs(density - iso) < iso * SHELL_TOL_WIDE
    if not mask.any():
        return None
    shell = esp[mask]
    vmin, vmax = float(shell.min()), float(shell.max())
    amp = math.ceil(max(abs(vmin), abs(vmax)) / step) * step
    return vmin, vmax, amp, int(mask.sum())


STATS_RE = re.compile(r"V_S,min=(\S+)\s+V_S,max=(\S+)\s+range=(\S+)"
                      r"(?:\s+iso=\S+\s+n=(\d+))?")


def stats_comment(stats, iso):
    """Statistics for the first cube line, so that --tcl-only finds them again.

    Eight decimals, not six. The values are reported to five, and six decimals
    put the stamp one digit away from that: a value like -0.018825 sits exactly
    on the rounding boundary, so reading it back and rounding a second time can
    land on -0.01883 where the full float gives -0.01882. Two extra digits move
    the stamp far enough from the boundary that the second rounding cannot
    disagree with the first. Cubes stamped with six decimals still parse.
    """
    vmin, vmax, amp, npts = stats
    return (f" | V_S,min={vmin:+.8f} V_S,max={vmax:+.8f} range={amp:.4f}"
            f" iso={iso:g} n={npts}")


def read_cube(path):
    """Gaussian cube -> float32 grid [i1, i2, i3].

    Only as much parser as the shell statistics need. It is needed for cubes
    that do not yet carry the statistics in their header (older runs) - the
    ~15 s of reading are still better than minutes of reconversion.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        fh.readline()
        fh.readline()
        natoms = abs(int(fh.readline().split()[0]))
        dims = [int(fh.readline().split()[0]) for _ in range(3)]
        for _ in range(natoms):
            fh.readline()
        values = np.asarray(fh.read().split(), dtype=np.float32)
    if values.size != dims[0] * dims[1] * dims[2]:
        raise ValueError(f"{path}: {values.size} values, expected "
                         f"{dims[0] * dims[1] * dims[2]}.")
    return values.reshape(dims)


def read_cube_atoms(path):
    """The atom block of a cube -> [(Z, x, y, z), ...], lengths in Bohr.

    The cube carries the atoms itself - which is why the orientation can be
    determined with --tcl-only as well, without finding the structure file
    again.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        fh.readline()
        fh.readline()
        natoms = abs(int(fh.readline().split()[0]))
        for _ in range(3):
            fh.readline()
        atoms = []
        for _ in range(natoms):
            p = fh.readline().split()
            # column 1 atomic number, column 2 nuclear charge, then x y z
            atoms.append((int(p[0]), float(p[2]), float(p[3]), float(p[4])))
    return atoms


def read_stats(cube_path):
    """Statistics from the first line of a cube written earlier."""
    with open(cube_path, "r", encoding="utf-8", errors="replace") as fh:
        m = STATS_RE.search(fh.readline())
    if not m:
        return None
    return (float(m.group(1)), float(m.group(2)), float(m.group(3)),
            int(m.group(4)) if m.group(4) else 0)


# ----------------------------------------------------------------------------
# Orientation of the three views
#
# This calculation is taken over from the PyMOL pipeline (render_esp.py,
# molecular_frame/view_matrix) and deliberately yields the same axes: only
# then do both image sets show the same thing and the comparison measures the
# viewer.
#
# The orientation used to sit in the Tcl template as a fixed rotation - "pi,
# then rotate x by -90". That assumes the molecule lies planar in the xy plane
# with the C-X axis pointing along -y. For the halobenzenes, as Turbomole
# writes them out, that happened to hold. For the substituted pyridines it did
# not: there the sigma view looked past the hole. It is hardly noticeable,
# because in every direction there is a round, coloured surface to be seen -
# and one takes the far side for the sigma hole.
# ----------------------------------------------------------------------------

HALOGENS = {9: "F", 17: "Cl", 35: "Br", 53: "I"}


def molecular_frame(atoms, halogen_index=None):
    """A reproducible molecular frame derived from the geometry.

    Returns (normal, axis, sigma_axis, label)
      normal      surface normal: the principal axis of smallest extent over
                  the heavy atoms (hydrogen wobbles and contributes nothing)
      axis        the C->halogen axis PROJECTED INTO THE PLANE; together with
                  normal it spans a right-handed frame and orients pi and edge
      sigma_axis  the TRUE, unprojected C->halogen axis
      label       where the axis comes from, for the message in the script

    Why two axes: for planar molecules they are identical and the projection
    is rounding cosmetics. But if the C-X bond points out of the best-fit
    plane, it really does turn the axis away - and the sigma view looked past
    the hole by exactly that angle.

    Without a halogen there is no sigma hole. The axis then takes the largest
    inertial extent, that is the long axis; there the sigma view is the look
    from the narrow side and says nothing about a hole.
    """
    coords = np.array([[a[1], a[2], a[3]] for a in atoms], dtype=float)
    znums = np.array([a[0] for a in atoms])

    heavy = coords[znums > 1]
    if len(heavy) < 3:
        heavy = coords
    _, _, vt = np.linalg.svd(heavy - heavy.mean(axis=0), full_matrices=False)
    normal = vt[2]                                  # smallest extent
    long_axis = vt[0]                               # largest extent

    axis, label = None, "longest inertial axis (no halogen)"
    hal_idx = [i for i, z in enumerate(znums) if z in HALOGENS]
    if hal_idx:
        hi = halogen_index if halogen_index in hal_idx else hal_idx[0]
        carbons = [i for i, z in enumerate(znums) if z == 6]
        if carbons:
            d = np.linalg.norm(coords[carbons] - coords[hi], axis=1)
            ci = carbons[int(np.argmin(d))]
            axis = coords[hi] - coords[ci]          # C -> X, points at the hole
            label = (f"C{ci + 1}->{HALOGENS[int(znums[hi])]}{hi + 1}")

    if axis is None:
        axis = long_axis.copy()

    axis = axis / np.linalg.norm(axis)
    normal = normal / np.linalg.norm(normal)
    sigma_axis = axis.copy()

    axis = axis - normal * float(np.dot(axis, normal))
    if np.linalg.norm(axis) < 1e-6:                 # C-X is perpendicular
        axis = long_axis
    axis = axis / np.linalg.norm(axis)

    return normal, axis, sigma_axis, label


def view_matrix(forward, up):
    """Rotation matrix for VMD's ``rotate_matrix``.

    ``forward`` points from the molecule to the camera, ``up`` points up in
    the image.

    VMD computes image coordinates as R * (world - centre). For the camera to
    look along forward, R has to map exactly that vector onto +z - so the ROWS
    of R are the camera basis. (PyMOL's set_view wants the same matrix
    transposed; that is why there is a .T over there and none here.)
    """
    z = np.asarray(forward, dtype=float)
    z = z / np.linalg.norm(z)
    up = np.asarray(up, dtype=float)
    up = up - z * float(np.dot(up, z))
    if np.linalg.norm(up) < 1e-8:                   # up parallel to forward
        alt = np.array([0.0, 0.0, 1.0])
        if abs(float(np.dot(alt, z))) > 0.9:
            alt = np.array([1.0, 0.0, 0.0])
        up = alt - z * float(np.dot(alt, z))
    up = up / np.linalg.norm(up)
    right = np.cross(up, z)
    right = right / np.linalg.norm(right)
    return np.array([right, up, z])


def tcl_matrix(m):
    """3x3 rotation -> 4x4 matrix as a Tcl list, the way molinfo expects it."""
    rows = []
    for i in range(3):
        rows.append("{%s 0.000000}" % " ".join(f"{v:.6f}" for v in m[i]))
    rows.append("{0.000000 0.000000 0.000000 1.000000}")
    return "{%s}" % " ".join(rows)


# Fallback without atoms: exactly the axes the earlier fixed rotations gave
# (pi = looking from +z, edge = rotate y by 90, sigma = rotate x by -90). That
# way nothing changes if the geometry is missing for once.
_FALLBACK = {
    "pi":    ((0.0, 0.0, 1.0), (0.0, 1.0, 0.0)),
    "edge":  ((-1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
    "sigma": ((0.0, -1.0, 0.0), (0.0, 0.0, 1.0)),
}


def view_matrices(atoms):
    """The three views as Tcl matrices plus a note on the sigma axis."""
    if not atoms:
        mats = {k: view_matrix(f, u) for k, (f, u) in _FALLBACK.items()}
        return {k: tcl_matrix(m) for k, m in mats.items()}, \
               "world axes (no geometry found)", None

    normal, axis, sigma_axis, label = molecular_frame(atoms)
    mats = {
        # looking perpendicular onto the plane, C-X axis pointing down
        "pi":    view_matrix(forward=normal, up=-axis),
        # looking in the plane, perpendicular to the C-X axis; C-X horizontal
        "edge":  view_matrix(forward=np.cross(normal, axis), up=normal),
        # looking in from outside along the TRUE C-X axis onto the sigma hole
        "sigma": view_matrix(forward=sigma_axis, up=normal),
    }
    # Tilt of the C-X axis against the best-fit plane: if it is large, sigma
    # and pi/edge diverge, and that should be on record.
    tilt = abs(90.0 - math.degrees(math.acos(
        min(1.0, abs(float(np.dot(sigma_axis, normal)))))))
    return {k: tcl_matrix(m) for k, m in mats.items()}, label, tilt


# ----------------------------------------------------------------------------
# Writing the cube
# ----------------------------------------------------------------------------

def write_cube(path, info, data, atoms, stride=1, comment=""):
    """Gaussian cube, all lengths in Bohr.

    It is written to <name>.part first and renamed at the end. A full cube is
    200 MB and takes minutes; if the run is aborted in that time - Ctrl-C, a
    closed window, a full disk - half a file with a valid header would
    otherwise be left behind. The next run takes that for finished, skips the
    conversion and renders an isosurface with its rear half missing. The
    rename is the moment the file comes into being - before that it does not
    exist under its name.
    """
    if stride > 1:
        data = data[::stride, ::stride, ::stride]
    n = data.shape
    vecs = info["vectors"]
    starts = [info["grid"][i][0] for i in range(3)]
    origin = info["origin"] + sum(starts[i] * vecs[i] for i in range(3))
    voxel = [info["grid"][i][1] * stride * vecs[i] for i in range(3)]

    tmp = path + ".part"
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(f"{comment or 'written by xyzToCubeToVMDVis.py'}\n")
        fh.write(f"{info.get('quantity') or 'volumetric data'} | "
                 f"{info.get('title', '')} | units: Bohr\n")
        fh.write(f"{len(atoms):5d} {origin[0]:12.6f} {origin[1]:12.6f} "
                 f"{origin[2]:12.6f}\n")
        for i in range(3):
            fh.write(f"{n[i]:5d} {voxel[i][0]:12.6f} {voxel[i][1]:12.6f} "
                     f"{voxel[i][2]:12.6f}\n")
        for (znum, x, y, z) in atoms:
            fh.write(f"{znum:5d} {float(znum):12.6f} {x:12.6f} {y:12.6f} "
                     f"{z:12.6f}\n")

        # Values: z fastest, 6 per line. A precompiled format pattern per z
        # row is considerably faster than a loop over all individual values.
        nz = n[2]
        full, rest = divmod(nz, 6)
        row_fmt = ("%13.5E" * 6 + "\n") * full + \
                  ("%13.5E" * rest + "\n" if rest else "")
        flat = data.reshape(n[0] * n[1], nz)
        buf = []
        for idx in range(flat.shape[0]):
            buf.append(row_fmt % tuple(flat[idx].tolist()))
            if len(buf) >= 4096:
                fh.write("".join(buf))
                buf.clear()
        fh.write("".join(buf))
    # Only here does <name>.cube exist. os.replace also replaces an existing
    # file and is atomic within one file system - there is no moment in which
    # the old cube is gone and the new one is not yet there.
    os.replace(tmp, path)
    return n


# ----------------------------------------------------------------------------
# VMD scene
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Colour ramps
#
# The same anchor colours as in the PyMOL pipeline (render_esp.py, RAMP_HEX),
# so that both image sets really show the same ramp and not just the same
# direction. VMD's built-in scales know only three anchor colours; the five of
# the rainbow are therefore written into the 1024 entries of the colour table
# via esp_ramp in esp_template.tcl.
#
# Red-white-blue stays VMD's own RWB for now (pure red and blue instead of
# PyMOL's slightly darker #d40000/#0030d4). Switching would change every
# standard set already rendered, if only slightly.
# ----------------------------------------------------------------------------

# The VMD scale underneath the ramp. Hard-wired, not selectable: RED is
# negative, BLUE positive - the convention of Politzer/Murray, of the
# literature and of the PyMOL pipeline. There used to be --color-scale for
# this. It was taken out, because a reversed ramp does not give an image that
# LOOKS different, it gives one that appears to show the opposite statement -
# and because the colour bar comes from matplotlib and would not have followed
# the reversal: the images would have been mirrored, the legend beside them
# not.
#
# With --rainbow, RGB is only the underlay; esp_ramp afterwards overwrites the
# colour table with the five anchor colours from RAMP_HEX.
COLOR_SCALE = {"redblue": "RWB", "rainbow": "RGB"}

RAMP_HEX = {
    "rainbow": ["#d40000", "#f0e000", "#00a000", "#00c8d4", "#0030d4"],
    "redblue": ["#d40000", "#ffffff", "#0030d4"],
}


def ramp_stops(name):
    """Anchor colours as a Tcl list {{r g b} {r g b} ...}, values 0..1."""
    hexes = RAMP_HEX.get(name)
    if not hexes:
        return ""
    out = []
    for h in hexes:
        h = h.lstrip("#")
        out.append("{%s}" % " ".join(
            f"{int(h[i:i+2], 16) / 255.0:.4f}" for i in (0, 2, 4)))
    return " ".join(out)


def write_vmd_script(path, rho_cube, esp_cube, esp_range, stats, iso=0.001,
                     opacity=0.50, scale="auto", fill=0.85,
                     sources="", atoms=None, rainbow=False):
    """Fills in esp_template.tcl.

    Via @@PLACEHOLDER@@ rather than str.format or string.Template: Tcl uses
    braces and dollar signs as syntax, and both mechanisms would choke on
    exactly that.

    rho_cube is loaded first and is therefore Volume 0, esp_cube Volume 1.
    """
    tpl_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), TEMPLATE)
    with open(tpl_path, "r", encoding="utf-8") as fh:
        text = fh.read()

    note = ""
    if stats:
        note = (f"\\nV_S,min = {stats[0]:+.5f}   V_S,max = {stats[1]:+.5f} a.u. "
                f"(grid points on the shell)")

    # Orientation from the geometry. Without atoms passed in, read them from
    # the cube - it carries them, with --tcl-only as well.
    if atoms is None:
        try:
            atoms = read_cube_atoms(
                os.path.join(os.path.dirname(os.path.abspath(path)), rho_cube))
        except Exception as err:
            print(f"  ! atoms not readable ({err}) - views on the world axes",
                  file=sys.stderr)
            atoms = None
    rot, axis_label, tilt = view_matrices(atoms)
    axis_note = f"\\nsigma axis: {axis_label}"
    if tilt is not None and tilt > 15.0:
        # If the bond is tilted, pi/edge follow the plane and sigma follows
        # the bond - then they are no longer three views of the same frame.
        axis_note += (f" ({tilt:.0f} degrees out of the best-fit plane; "
                      f"pi/edge follow the plane)")

    for key, value in {
        "@@STAMP@@": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "@@SOURCES@@": sources,
        "@@RHO_CUBE@@": rho_cube,
        "@@ESP_CUBE@@": esp_cube,
        "@@VOL_RHO@@": "0",
        "@@VOL_ESP@@": "1",
        "@@ISO@@": f"{iso:g}",
        "@@RANGE@@": f"{esp_range:.4f}",
        "@@RANGE_NEG@@": f"{-esp_range:.4f}",
        "@@KJ@@": f"{-esp_range * HARTREE_TO_KJ:.0f} .. "
                  f"{esp_range * HARTREE_TO_KJ:+.0f}",
        "@@OPACITY@@": f"{opacity:.6f}",
        "@@SCALE@@": scale if scale == "auto" else f"{scale:g}",
        "@@FILL@@": f"{fill:g}",
        "@@COLORSCALE@@": COLOR_SCALE["rainbow" if rainbow else "redblue"],
        "@@RAMP_STOPS@@": ramp_stops("rainbow" if rainbow else None),
        "@@STATS@@": note + axis_note,
        "@@ROT_PI@@": rot["pi"],
        "@@ROT_EDGE@@": rot["edge"],
        "@@ROT_SIGMA@@": rot["sigma"],
        "@@AXIS_LABEL@@": axis_label,
    }.items():
        text = text.replace(key, value)

    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


def classify(quantity):
    q = (quantity or "").lower()
    return "esp" if "potential" in q else "density" if "density" in q else None


# ----------------------------------------------------------------------------
# Main program
# ----------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        description="Convert Turbomole pointval grids to Gaussian cube and "
                    "write a VMD scene.",
        epilog="Example: python xyzToCubeToVMDVis.py "
               "--struct brombenzol_aro_opti.mol td.xyz tp.xyz")
    p.add_argument("grids", nargs="+", help="td.xyz tp.xyz (or td.cube "
                                            "tp.cube with --tcl-only)")
    p.add_argument("--struct", "-s", help="structure file for the cube header: "
                                          ".xyz, .mol, .sdf")
    p.add_argument("--struct-unit", choices=["angstrom", "bohr"],
                   default="angstrom", help="relevant for .xyz only")
    p.add_argument("--outdir", "-o", help="default: next to the input")
    p.add_argument("--stride", type=int, default=1,
                   help="every n-th grid point (2 => 8x smaller)")
    p.add_argument("--quiet", "-q", action="store_true")

    g = p.add_argument_group("VMD scene (esp.tcl)")
    g.add_argument("--tcl-only", action="store_true",
                   help="only write esp.tcl again, leave the cubes untouched "
                        "(fractions of a second instead of minutes)")
    g.add_argument("--no-vmd", action="store_true", help="the cube files only")
    g.add_argument("--esp-range", default="auto",
                   help="half width of the colour scale in a.u., or 'auto' "
                        "(default): from V_S,min/V_S,max on the shell")
    g.add_argument("--iso", type=float, default=0.001, help="isovalue (a.u.)")
    g.add_argument("--opacity", type=float, default=0.50,
                   help="opacity of the isosurface, 0..1")
    g.add_argument("--scale", default="auto",
                   help="zoom: a number, or 'auto' (default) from the size of "
                        "the molecule and the window height")
    g.add_argument("--fill", type=float, default=0.85,
                   help="fraction of the window height for the molecule with "
                        "--scale auto (default 0.85)")
    g.add_argument("--rainbow", action="store_true",
                   help="rainbow ramp instead of red-white-blue. Writes "
                        "esp_rainbow.tcl, so that the standard scene is "
                        "kept.")
    args = p.parse_args(argv)

    verbose = not args.quiet
    if args.tcl_only and args.no_vmd:
        p.error("--tcl-only and --no-vmd together leave nothing to do.")
    if not args.tcl_only and not args.struct:
        p.error("--struct is needed: the cube header requires an atom block, "
                "and the pointval file contains no atoms.")

    if verbose:
        print("=" * 66)
        print("xyzToCubeToVMDVis.py - pointval -> Cube -> VMD")
        print("=" * 66)

    outdir = args.outdir or os.path.dirname(os.path.abspath(args.grids[0]))
    cubes, stats = {}, None

    if args.tcl_only:
        for gpath in args.grids:
            base = os.path.splitext(os.path.basename(gpath))[0]
            cube = os.path.join(args.outdir or
                                os.path.dirname(os.path.abspath(gpath)),
                                base + ".cube")
            if not os.path.exists(cube):
                p.error(f"{cube} is missing - run once without --tcl-only.")
            with open(cube, "r", encoding="utf-8", errors="replace") as fh:
                fh.readline()
                kind = classify(fh.readline())
            if kind:
                cubes[kind] = cube
                if kind == "esp":
                    stats = read_stats(cube)
            if verbose:
                print(f"[1] {os.path.basename(cube)} -> {kind}")
    else:
        atoms = read_structure(args.struct, unit=args.struct_unit)
        if verbose:
            print(f"[1] structure: {args.struct} -> {len(atoms)} atoms "
                  f"({args.struct_unit} -> Bohr)")

        grids = {}
        for gpath in args.grids:
            if verbose:
                print(f"[2] grid file: {gpath}")
            info, data = read_values(gpath, verbose=verbose)
            if verbose:
                n = [info["grid"][i][2] for i in range(3)]
                print(f"    {n[0]} x {n[1]} x {n[2]}, "
                      f"delta = {info['grid'][0][1]} Bohr, "
                      f"'{info['quantity'] or 'unknown'}', "
                      f"{data.min():+.4g} .. {data.max():+.4g}")
            grids[classify(info["quantity"]) or gpath] = (info, data, gpath)

        if "density" in grids and "esp" in grids:
            stats = shell_range(grids["density"][1], grids["esp"][1], args.iso)
            if verbose and stats:
                v0, v1, amp, npts = stats
                print(f"[3] shell rho = {args.iso:g} a.u.: {npts:,} points")
                print(f"    V_S,min = {v0:+.5f} a.u. "
                      f"({v0 * HARTREE_TO_KJ:+7.1f} kJ/(mol*e))")
                print(f"    V_S,max = {v1:+.5f} a.u. "
                      f"({v1 * HARTREE_TO_KJ:+7.1f} kJ/(mol*e))")
                print(f"    -> colour scale +/- {amp:.4f} a.u.")

        for kind, (info, data, gpath) in grids.items():
            base = os.path.splitext(os.path.basename(gpath))[0]
            od = args.outdir or os.path.dirname(os.path.abspath(gpath))
            os.makedirs(od, exist_ok=True)
            out = os.path.join(od, base + ".cube")
            comment = f"{info['quantity'] or base} - from {os.path.basename(gpath)}"
            if kind == "esp" and stats:
                comment += stats_comment(stats, args.iso)
            shape = write_cube(out, info, data, atoms, args.stride, comment)
            if kind in ("density", "esp"):
                cubes[kind] = out
            if verbose:
                print(f"    -> {out}  ({'x'.join(map(str, shape))}, "
                      f"{os.path.getsize(out) / 1024 ** 2:.1f} MB)")

    if args.no_vmd:
        return 0

    if "esp" not in cubes or "density" not in cubes:
        print("    ! density or potential missing - no esp.tcl. The density "
              "provides the surface, the potential the colour.",
              file=sys.stderr)
        return 1

    if str(args.esp_range).lower() == "auto":
        if stats:
            rng = stats[2]
        else:
            rng = 0.035
            print("    ! colour scale not determinable - using +/- 0.035 a.u.",
                  file=sys.stderr)
    else:
        rng = float(args.esp_range)

    # Its own file name, otherwise a rainbow run overwrites the scene of the
    # red-white-blue run - the images are kept apart for the same reason.
    tcl = os.path.join(outdir, "esp_rainbow.tcl" if args.rainbow else "esp.tcl")
    write_vmd_script(tcl, os.path.basename(cubes["density"]),
                     os.path.basename(cubes["esp"]), rng, stats, iso=args.iso,
                     opacity=args.opacity,
                     scale=(args.scale if str(args.scale) == "auto"
                            else float(args.scale)),
                     fill=args.fill,
                     rainbow=args.rainbow,
                     sources=", ".join(os.path.basename(g) for g in args.grids))
    if verbose:
        print(f"[4] VMD scene: {tcl}   (colour scale +/- {rng:.4f} a.u.)")
        print(f"    start with:  vmd -e {os.path.basename(tcl)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

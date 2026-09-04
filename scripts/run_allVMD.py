#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_allVMD.py - batch run of the VMD pipeline

    python run_allVMD.py                 # self test on reference/
    python run_allVMD.py --root ../sandbox

Without --root the script runs on the reference dataset that ships with the
repository. That is the self test: it writes to images_check/ and
summary_check_<time>_<date>.csv so that the committed reference images stay
untouched, and afterwards you compare your own numbers with those in
reference/.

It walks a directory tree and does everything for each molecule folder in one
go: pointval -> cube, write esp.tcl, render the three views with Tachyon, plus
the colour bar and settings.txt. At the end there is a dated
summary_<HH-MM>_<DD-MM-YYYY>.csv with the statistics and the render parameters; a
--rainbow run writes its own summary_rainbow_<time>_<date>.csv, so no run overwrites
the summary of another.

A molecule folder is any directory holding

    td.xyz  + tp.xyz     (Turbomole output, gets converted)
or  td.cube + tp.cube    (already converted)

plus a structure file (.mol, .sdf or .xyz).

Order of work: analyse first, then render ONCE
----------------------------------------------
V_S,min and V_S,max are in the cube files, nothing has to be rendered for
them. So the script first collects the statistics of all molecules and then
renders exactly one pass.

    --esp-range auto     every molecule on its own scale (default)
    --esp-range 0.035    a fixed value
    --esp-range common   a common scale, straight from this run

The images are comparable only with ONE scale for all of them. The PyMOL
pipeline renders twice for that (--two-pass) - not here: at the end of every
run the console shows the smallest value that covers all molecules, together
with a ready-made command line. You read it, decide for yourself which scale
the figure should show, and start the run that counts. One image costs one to
two minutes in VMD; an automatic first set would be waste.

For a direct comparison with the PyMOL pipeline, take the value from its
summary.csv.

Not included
------------
No sigma hole, no trilinear interpolation along the C-X axis - that lives in
the sister project and is not implemented a second time here. V_S,min and
V_S,max are, though, and from the same grid points in the band
rho = iso +/- 12 % as over there: the same calculation, the same numbers. The
PyMOL values have to be quoted only for the sigma hole.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import fnmatch
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import xyzToCubeToVMDVis as conv          # noqa: E402
import render_espVMD as rend              # noqa: E402
from constants import HARTREE_TO_KJ       # noqa: E402

STRUCT_EXT = (".mol", ".sdf", ".xyz")
GRID_NAMES = ("td", "tp")
# reference/ sits next to scripts/, not in the current working directory.
# Without --root the run is therefore the self test, no matter where it was
# called from.
DEFAULT_ROOT = os.path.normpath(os.path.join(_HERE, "..", "reference"))


# ----------------------------------------------------------------------------
# Finding molecule folders
# ----------------------------------------------------------------------------

def discover(root):
    """All molecule folders below ``root``, sorted by name."""
    found = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames
                       if d not in ("images", "images_check", "__pycache__",
                                    ".git", "_to_delete")]
        names = set(filenames)
        has_raw = {"td.xyz", "tp.xyz"} <= names
        has_cubes = {"td.cube", "tp.cube"} <= names
        if not (has_raw or has_cubes):
            continue
        struct = find_structure(dirpath,
                                exclude=[os.path.join(dirpath, "td.xyz"),
                                         os.path.join(dirpath, "tp.xyz")])
        if struct is None and not has_cubes:
            print(f"  ! {dirpath}: grids found, but no structure file "
                  f"({'/'.join(STRUCT_EXT)}) - skipped")
            continue
        found.append({"dir": dirpath, "struct": struct,
                      "has_raw": has_raw, "has_cubes": has_cubes})
    return sorted(found, key=lambda e: e["dir"])


def find_structure(folder, exclude):
    """The structure file in ``folder``, or None if there is none.

    Aborts when the folder holds more than one candidate that is not simply
    the same structure in a second format.
    """
    excl = {os.path.abspath(p) for p in exclude if p}
    found = []
    for name in sorted(os.listdir(folder)):
        stem, ext = os.path.splitext(name)
        if ext.lower() not in STRUCT_EXT:
            continue
        path = os.path.join(folder, name)
        if os.path.abspath(path) in excl:
            continue
        if stem in GRID_NAMES:
            continue                       # td.xyz / tp.xyz are data
        found.append((stem, ext.lower(), path))

    if not found:
        return None

    # One structure in several formats is fine - brombenzol_aro_opti.mol next
    # to brombenzol_aro_opti.xyz is the same geometry twice. Anything else is
    # ambiguous, and the script must not pick for you: it would take the first
    # in sort order, and a stale left-over file sorting first silently pairs
    # the wrong geometry with the grids. That produces a molecule floating
    # beside its own isosurface, and every number in the run - the extrema,
    # the atom labels, the sigma hole - refers to atoms that are not where the
    # density says they are.
    stems = {s.lower() for s, _, _ in found}
    exts = [e for _, e, _ in found]
    if len(stems) > 1 or len(exts) != len(set(exts)):
        listing = "\n".join(f"      {os.path.basename(p)}" for _, _, p in found)
        raise SystemExit(
            f"{folder}: located {len(found)} structure files, keep the correct "
            f"one to proceed.\n{listing}\n"
            f"    They are not one structure in two formats, so the script "
            f"cannot tell which geometry belongs to the grids.")

    # .mol/.sdf before .xyz: they carry bond information, and a bare .xyz
    # easily collides with the Turbomole grid files of the same extension.
    for ext in STRUCT_EXT:
        for _, e, path in found:
            if e == ext:
                return path
    return found[0][2]


def cube_dims(path):
    """Grid dimensions from the cube header, without reading the data."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for _ in range(3):
            fh.readline()
        return "x".join(fh.readline().split()[0] for _ in range(3))


# ----------------------------------------------------------------------------
# Step 1: cubes and statistics
# ----------------------------------------------------------------------------

def ensure_cubes(entry, args):
    """Provide the cubes and return (V_S,min, V_S,max, scale, points).

    Three routes, in this order:
      1. Cubes missing or --force-convert: build them from the pointval grids,
         which produces the statistics along the way anyway.
      2. Cubes present and carrying the statistics in the header: take them
         over directly.
      3. Cubes present but without statistics (older runs): read the cubes and
         recompute - about 15 s instead of minutes for a reconversion.
    """
    folder = entry["dir"]
    cubes = {t: os.path.join(folder, f"{t}.cube") for t in GRID_NAMES}
    have = all(os.path.exists(c) for c in cubes.values())

    if have and not args.force_convert:
        stats = conv.read_stats(cubes["tp"])
        if stats:
            print("    cubes present, statistics from the cube header")
            return cubes, stats
        print("    cubes present but without statistics - recomputing them")
        density = conv.read_cube(cubes["td"])
        esp = conv.read_cube(cubes["tp"])
        cube_atoms, cube_origin, cube_voxel = conv.read_cube_frame(cubes["td"])
        conv.check_alignment(density, cube_atoms, cube_origin, cube_voxel,
                             label=os.path.basename(folder))
        return cubes, conv.shell_range(density, esp, args.iso)

    if entry["struct"] is None:
        raise SystemExit(f"{folder}: the structure file is missing, without "
                         f"it no cube header can be written.")
    atoms = conv.read_structure(entry["struct"], unit=args.struct_unit)
    grids = {}
    for tag in GRID_NAMES:
        raw = os.path.join(folder, f"{tag}.xyz")
        if not os.path.exists(raw):
            raise SystemExit(f"{folder}: neither {tag}.cube nor {tag}.xyz")
        print(f"    converting {tag}.xyz -> {tag}.cube (stride {args.stride})")
        grids[tag] = conv.read_values(raw, verbose=False)

    # The structure file is married to the grid right here - so this is where
    # a mismatch has to be caught, before minutes of rendering go into it.
    info_d = grids["td"][0]
    origin_d = info_d["origin"] + sum(info_d["grid"][i][0] * info_d["vectors"][i]
                                      for i in range(3))
    delta_d = np.array([info_d["grid"][i][1] for i in range(3)])
    conv.check_alignment(grids["td"][1], atoms, origin_d, delta_d,
                         label=os.path.basename(folder))

    stats = conv.shell_range(grids["td"][1], grids["tp"][1], args.iso)
    for tag in GRID_NAMES:
        info, data = grids[tag]
        comment = f"{info['quantity'] or tag} - from {tag}.xyz"
        if tag == "tp" and stats:
            comment += conv.stats_comment(stats, args.iso)
        conv.write_cube(cubes[tag], info, data, atoms, args.stride, comment)
    return cubes, stats


# ----------------------------------------------------------------------------
# summary.csv
# ----------------------------------------------------------------------------

FIELDS = ["molecule", "structure", "grid", "iso_au", "shell_points",
          "VS_min_au", "VS_max_au", "VS_min_kJ", "VS_max_kJ",
          "esp_range_used_au", "esp_range_mode", "colormap", "color_scale",
          "opacity", "resolution_px", "renderer",
          "backgrounds", "views"]


def summary_name(is_reference, rainbow, when=None):
    """summary[_check][_rainbow]_HH-MM_DD-MM-YYYY.csv

    Time-stamped on purpose, so that a later run does not silently overwrite
    the summary of an earlier one. The date alone was not enough: two runs on
    the same day - the usual case while a parameter is being settled - landed
    on the same name again. The minute is the coarsest unit that separates
    them in practice.

    The _rainbow part matters just as much: the scene and the images already
    carry that suffix, and without it here the CSV was the one file a
    --rainbow run clobbered - a run over one molecule replaced the summary of
    the whole set, and the loss was invisible until someone opened the file.

    The same name is built in the sister project (run_all.py), so summaries
    from the two pipelines can be filed next to each other.
    """
    stamp = (when or datetime.datetime.now()).strftime("%H-%M_%d-%m-%Y")
    parts = ["summary"]
    if is_reference:
        parts.append("check")
    if rainbow:
        parts.append("rainbow")
    parts.append(stamp)
    return "_".join(parts) + ".csv"


def write_summary(path, rows, common=None, advice=None):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for r in rows:
            s = r.get("stats")
            w.writerow({
                "molecule": r["prefix"],
                "structure": (os.path.basename(r["struct"]) if r["struct"]
                              else ""),
                "grid": r.get("grid", ""),
                "iso_au": r["iso"],
                "shell_points": s[3] if s else "",
                "VS_min_au": f"{s[0]:.5f}" if s else "",
                "VS_max_au": f"{s[1]:.5f}" if s else "",
                "VS_min_kJ": f"{s[0] * HARTREE_TO_KJ:.1f}" if s else "",
                "VS_max_kJ": f"{s[1] * HARTREE_TO_KJ:.1f}" if s else "",
                "esp_range_used_au": f"{r['rng']:.4f}",
                "esp_range_mode": r["mode"],
                "colormap": r["colormap"],
                "color_scale": r["color_scale"],
                "opacity": r["opacity"],
                "resolution_px": (f"{r['size'][0]}x{r['size'][1]}"
                                  if r.get("size") else ""),
                "renderer": r.get("renderer", ""),
                "backgrounds": " ".join(r["backgrounds"]),
                # Per view, what it was actually produced with - on a Tachyon
                # crash the fallback kicks in, and that must not disappear
                # from the record.
                "views": ";".join(f"{v}:{r['made'][v]}" for v in rend.VIEWS
                                  if v in r.get("made", {})),
            })
        if common is not None:
            fh.write(f"# common colour scale for all molecules: "
                     f"+/- {common:.4f} a.u.\n")
        if advice is not None:
            # In the file as well, not only in the console: whoever opens the
            # CSV later would otherwise not see that the images sit on
            # different scales, nor what a comparable set would cost.
            fh.write(f"# each molecule on its own scale; a comparable set "
                     f"needs --esp-range "
                     f"{advice:.4f}\n")


# ----------------------------------------------------------------------------
# Recommendation for the next run
# ----------------------------------------------------------------------------

def again_with(raw_argv, value):
    """The same call once more, only with a fixed --esp-range."""
    keep, skip = [], False
    for a in raw_argv:
        if skip:
            skip = False
            continue
        if a == "--esp-range":
            skip = True
            continue
        if a.startswith("--esp-range="):
            continue
        if a == "--no-render":
            # The suggestion aims at an image set - offering it with
            # --no-render would be the one line guaranteed to render nothing.
            continue
        keep.append(f'"{a}"' if " " in a else a)
    return " ".join(["python run_allVMD.py"] + keep
                    + [f"--esp-range {value:.4f}"])


def advise(mode, common, needed, rows, raw_argv):
    """Say which scale would give a comparable set - and how to set it.

    Instead of a second render pass (which is how the PyMOL pipeline does it,
    with --two-pass) there is only the suggestion here. Which scale is
    sensible depends on what the figure is meant to show: the smallest one
    that clips nothing is rarely the most informative - one molecule with a
    very negative oxygen otherwise flattens the contrast of all the others.
    That decision belongs to the user, not to the script.
    """
    if needed is None:
        return
    used = sorted({r["rng"] for r in rows})
    print()
    print("-" * 70)
    if mode == "auto":
        if len(used) == 1:
            print(f"All molecules ended up on +/- {used[0]:.4f} a.u. "
                  f"anyway -")
            print("so the images are already comparable. To fix it in place:")
        else:
            print("Every molecule has its own scale - "
                  f"{', '.join(f'{v:.4f}' for v in used)} a.u.")
            print("The images may NOT be laid side by side like this.")
            print(f"Smallest value that covers all of them: +/- {needed:.4f} a.u.")
            print("To render a comparable set:")
        print()
        print(f"    {again_with(raw_argv, needed)}")
    elif common is not None and common + 1e-12 < needed:
        low = [r["prefix"] for r in rows
               if r["stats"] and r["stats"][2] > common + 1e-12]
        print(f"! The fixed scale +/- {common:.4f} a.u. clips.")
        print(f"  {', '.join(low)} needs up to +/- {needed:.4f} a.u.; "
              f"everything above")
        print("  runs into saturation and can no longer be told apart in the "
              "image.")
        print("  On purpose? Then all is well. Otherwise:")
        print()
        print(f"    {again_with(raw_argv, needed)}")
    else:
        print(f"Scale +/- {common:.4f} a.u. covers all molecules "
              f"(+/- {needed:.4f} would be needed).")
    print("-" * 70)


# ----------------------------------------------------------------------------
# Main program
# ----------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        description="Batch conversion and rendering of the VMD ESP pipeline.")
    p.add_argument("--root", default=DEFAULT_ROOT,
                   help="directory tree holding the molecule folders "
                        "(default: the repository's reference/ - the self "
                        "test)")
    p.add_argument("--only", nargs="+", metavar="NAME",
                   help="these folders only, wildcards allowed: "
                        "--only paracetamol '*benzol'")

    g = p.add_argument_group("Conversion")
    g.add_argument("--stride", type=int, default=1,
                   help="every n-th grid point (2 => 8x smaller cubes). "
                        "Applies only when building, not to existing cubes.")
    g.add_argument("--struct-unit", choices=["angstrom", "bohr"],
                   default="angstrom")
    g.add_argument("--force-convert", action="store_true",
                   help="write the cubes again, even if they already exist")

    g = p.add_argument_group("Scene")
    g.add_argument("--esp-range", default="auto",
                   help="'auto' (default, every molecule on its own scale), "
                        "a fixed value in a.u. - the console prints a "
                        "suggestion for it at the end of every run - or "
                        "'common' (a common scale, straight from this run)")
    g.add_argument("--iso", type=float, default=0.001, help="isovalue in a.u.")
    g.add_argument("--opacity", type=float, default=0.50,
                   help="opacity of the isosurface, 0..1")
    g.add_argument("--scale", default="auto", help="zoom, a number or 'auto'")
    g.add_argument("--fill", type=float, default=0.85,
                   help="fraction of the image height for the molecule with --scale auto")
    g.add_argument("--stick-size", type=float,
                   default=conv.STICK_SIZE_DEFAULT,
                   help="bond radius of the Licorice skeleton in Angstrom "
                        "(default: %(default)g)")
    g.add_argument("--rainbow", action="store_true",
                   help="rainbow ramp instead of red-white-blue. Writes "
                        "esp_rainbow.tcl and an image set of its own, "
                        "<molecule>_rainbow_*; the standard set is kept.")

    g = p.add_argument_group("Rendering")
    g.add_argument("--no-render", action="store_true",
                   help="only convert and write esp.tcl")
    g.add_argument("--vmd", help="path to vmd.exe, if it is not on the PATH")
    g.add_argument("--res", default="1600x1280", help="image size")
    g.add_argument("--backgrounds", nargs="+", default=["white"],
                   metavar="COLOUR",
                   help="background colours, e.g. white black. From two "
                        "colours on, the files get a suffix: "
                        "<molecule>_pi_black.png")
    g.add_argument("--keep-tga", action="store_true")
    g.add_argument("--dpi", type=int, default=300, help="resolution of the bar")
    g.add_argument("--images-dir", default=None,
                   help="target folder inside each molecule folder (default: "
                        "'images', 'images_check' in the self test)")
    g.add_argument("--summary", default=None,
                   help="path of the CSV (default: "
                        "<root>/summary_HH-MM_DD-MM-YYYY.csv, with _check "
                        "on the self test and _rainbow with --rainbow)")
    args = p.parse_args(argv)
    # For the recommendation at the end: the same call, only with a different
    # scale.
    raw_argv = list(sys.argv[1:]) if argv is None else list(argv)

    # The self test must not overwrite anything committed: its own image
    # folder, its own scene file, its own CSV. Otherwise there is no way to
    # tell afterwards whether the reference is still the reference.
    is_reference = os.path.abspath(args.root) == os.path.abspath(DEFAULT_ROOT)
    # Its own scene name per ramp, otherwise a rainbow run overwrites the
    # scene of the red-white-blue run - the images are kept apart for the same
    # reason (<molecule>_rainbow_*).
    scene_name = ("esp_check" if is_reference else "esp") \
        + ("_rainbow" if args.rainbow else "") + ".tcl"
    if args.images_dir is None:
        args.images_dir = "images_check" if is_reference else "images"
    # For the summary only - the ramp itself is fixed, see COLOR_SCALE in
    # xyzToCubeToVMDVis.py.
    color_scale = conv.COLOR_SCALE["rainbow" if args.rainbow else "redblue"]

    # Make it absolute before the first chdir: rendering runs inside the
    # molecule folder, where a relative path would point somewhere else.
    if args.vmd:
        args.vmd = os.path.abspath(args.vmd)
    args.root = os.path.abspath(args.root)

    print("=" * 70)
    print("run_allVMD.py - batch run of the VMD ESP pipeline")
    print("=" * 70)
    if is_reference:
        print("  Reference dataset (self test).")
        print(f"  Writes to '{args.images_dir}/', {scene_name} and "
              f"{os.path.basename(summary_name(True, args.rainbow))},")
        print("  so that the committed reference files stay untouched.")
        print("  Afterwards compare your numbers with reference/summary.csv")
        print("  and your images with reference/*/images/.")

    entries = discover(args.root)
    if args.only:
        keep = [e for e in entries
                if any(fnmatch.fnmatch(os.path.basename(e["dir"]).lower(),
                                       pat.lower()) for pat in args.only)]
        if len(entries) - len(keep):
            print(f"  --only: {len(entries) - len(keep)} folder(s) skipped")
        entries = keep
    if not entries:
        raise SystemExit(f"No molecule folders found under '{args.root}'"
                         + (" (with --only)." if args.only else "."))
    print(f"{len(entries)} molecule folder(s):")
    for e in entries:
        s = os.path.basename(e["struct"]) if e["struct"] else "cubes only"
        print(f"  - {e['dir']}  ({s})")

    # --- Step 1: cubes and statistics ---------------------------------------
    print("\n" + "-" * 70)
    print("Step 1: cubes and statistics")
    print("-" * 70)
    # One unusable molecule must not cost the other eight their run: it is
    # skipped, named again at the end, and the exit code says something went
    # wrong. Every other error still stops the run - a mismatched structure is
    # the one failure that is both survivable and worth continuing past.
    ok_entries, skipped_bad = [], []
    for e in entries:
        name = os.path.basename(os.path.normpath(e["dir"]))
        print(f"\n[{name}]")
        try:
            e["cubes"], e["stats"] = ensure_cubes(e, args)
        except conv.StructureGridMismatch as err:
            print(f"    ! {err}")
            print(f"    ! {name} skipped - no scene and no images written.")
            skipped_bad.append(name)
            continue
        e["grid"] = cube_dims(e["cubes"]["td"])
        if e["stats"]:
            v0, v1, amp, n = e["stats"]
            print(f"    V_S,min = {v0:+.5f}   V_S,max = {v1:+.5f} a.u.  "
                  f"({n} points)  -> +/- {amp:.4f}")
        else:
            print("    ! no shell found - statistics missing")
        ok_entries.append(e)
    entries = ok_entries
    if not entries:
        raise SystemExit("No molecule could be processed.")

    # --- Fix the colour scale -----------------------------------------------
    # needed = the scale that covers ALL molecules of this run. It is always
    # computed, in auto mode too: the recommendation at the end comes out of
    # it. That is why this script needs no second render pass - the suggestion
    # is in the console, the user decides, and the next call with
    # --esp-range <value> is the pass that counts.
    ranges = [e["stats"][2] for e in entries if e["stats"]]
    needed = max(ranges) if ranges else None
    common = needed
    mode = str(args.esp_range).lower()
    if mode == "auto":
        print("\nColour scale per molecule from its own statistics (auto)")
    elif mode == "common":
        if common is None:
            raise SystemExit("No statistics - --esp-range needs a fixed "
                             "value.")
        print(f"\nCommon colour scale: +/- {common:.4f} a.u. "
              f"(largest value of {len(ranges)} molecule(s))")
    else:
        try:
            common = float(args.esp_range)
        except ValueError:
            raise SystemExit(f"--esp-range: 'auto', 'common' or a number, "
                             f"not '{args.esp_range}'.")
        print(f"\nFixed colour scale: +/- {common:.4f} a.u.")

    # --- Step 2: scene and images -------------------------------------------
    print("\n" + "-" * 70)
    print(f"Step 2: write {scene_name}" +
          ("" if args.no_render else " and render"))
    print("-" * 70)

    rows = []
    for e in entries:
        name = os.path.basename(os.path.normpath(e["dir"]))
        print(f"\n[{name}]")
        if mode == "auto":
            rng = e["stats"][2] if e["stats"] else 0.035
        else:
            rng = common

        tcl = os.path.join(e["dir"], scene_name)
        conv.write_vmd_script(
            tcl,
            rho_cube=os.path.basename(e["cubes"]["td"]),
            esp_cube=os.path.basename(e["cubes"]["tp"]),
            esp_range=rng, stats=e["stats"], iso=args.iso,
            opacity=args.opacity,
            scale=(args.scale if str(args.scale) == "auto"
                   else float(args.scale)),
            fill=args.fill, rainbow=args.rainbow,
            stick_size=args.stick_size,
            sources=", ".join(os.path.basename(c)
                              for c in e["cubes"].values()))
        print(f"    -> {tcl}   (colour scale +/- {rng:.4f} a.u.)")

        row = {"prefix": name, "struct": e["struct"], "grid": e["grid"],
               "stats": e["stats"], "rng": rng, "iso": args.iso, "mode": mode,
               "colormap": "rainbow" if args.rainbow else "redblue",
               "color_scale": color_scale, "opacity": args.opacity,
               "backgrounds": args.backgrounds,
               "made": {}, "size": None, "renderer": ""}

        if not args.no_render:
            # render_espVMD works in the current directory - esp.tcl, the
            # cubes and images/ all sit next to each other there.
            cwd = os.getcwd()
            try:
                os.chdir(e["dir"])
                res = rend.render_all(
                    outdir=args.images_dir, prefix=name, iso=args.iso,
                    rng=rng, stats=e["stats"], vmd=args.vmd, res=args.res,
                    backgrounds=args.backgrounds,
                    scene=scene_name, rainbow=args.rainbow,
                    keep_tga=args.keep_tga, dpi=args.dpi, verbose=True)
            finally:
                os.chdir(cwd)
            row.update(made=res["made"], size=res["size"],
                       renderer=res["renderer"])
        rows.append(row)

    summary = args.summary or os.path.join(
        args.root, summary_name(is_reference, args.rainbow))
    write_summary(summary, rows,
                  common=None if mode == "auto" else common,
                  advice=needed if mode == "auto" else None)

    # --- Wrap-up ------------------------------------------------------------
    print("\n" + "-" * 70)
    print(f"{'Molecule':<22}{'V_S,min':>10}{'V_S,max':>10}{'Scale':>9}"
          f"  Views")
    print("-" * 70)
    for r in rows:
        s = r["stats"]
        v0 = f"{s[0]:>+10.4f}" if s else f"{'-':>10}"
        v1 = f"{s[1]:>+10.4f}" if s else f"{'-':>10}"
        views = ", ".join(v for v in rend.VIEWS if v in r["made"]) or "-"
        print(f"{r['prefix']:<22}{v0}{v1}{r['rng']:>9.4f}  {views}")
    print("-" * 70)
    if mode != "auto" and common is not None:
        print(f"Common scale: +/- {common:.4f} a.u.")
    print(f"CSV: {summary}")
    if is_reference:
        print("Now compare with the committed reference/summary_*.csv "
              "and reference/*/images/.")
    advise(mode, common, needed, rows, raw_argv)
    if skipped_bad:
        print(f"\n! skipped, structure and grid do not match: "
              f"{', '.join(skipped_bad)}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

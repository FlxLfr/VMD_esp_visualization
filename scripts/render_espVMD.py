#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
render_espVMD.py - the standard image set from a finished VMD scene

    python render_espVMD.py            # render, convert, colour bar
    python render_espVMD.py --no-vmd   # only TGA -> PNG and the colour bar

Start it from the molecule folder, where esp.tcl, td.cube and tp.cube live.
The result goes to images/:

    <molecule>_pi.png  _edge.png  _sigma.png  _colorbar.png  _settings.txt

The same set of files as in the PyMOL pipeline, on purpose: only that way can
the two image sets be laid side by side and really compared.

Division of labour: render_esp.tcl does the VMD part (views, Tachyon), this
script calls it, converts VMD's TGA to PNG and produces the colour bar with
matplotlib. The bar deliberately does not come from VMD - PyMOL cannot do it
either and makes it with matplotlib in render_esp.py as well.

Needs matplotlib and Pillow.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import time
import datetime

VIEWS = ("pi", "edge", "sigma")
STATS_RE = re.compile(r"V_S,min=(\S+)\s+V_S,max=(\S+)\s+range=(\S+)\s+iso=(\S+)")
HARTREE_TO_KJ = 2625.4996


def find_vmd(explicit=None):
    """Find vmd.exe: argument, PATH, VMDDIR, then the usual paths.

    The Windows installer sets VMDDIR to the installation directory - that is
    the most reliable hit, because it holds regardless of where the program
    was installed. The installer does NOT extend the PATH.
    """
    if explicit:
        return explicit
    found = shutil.which("vmd") or shutil.which("vmd.exe")
    if found:
        return found

    bases = []
    if os.environ.get("VMDDIR"):
        bases.append(os.environ["VMDDIR"])
    bases += [r"C:\Program Files\VMD",
              r"C:\Program Files (x86)\VMD",
              r"C:\Program Files\University of Illinois\VMD",
              r"C:\Program Files (x86)\University of Illinois\VMD",
              "/usr/local/bin", "/opt/vmd/bin"]
    for base in bases:
        for name in ("vmd.exe", "vmd"):
            cand = os.path.join(base, name)
            if os.path.exists(cand):
                return cand
    return None


def run_vmd(vmd, outdir, prefix, views, w, h, ao, opaque, label,
            snapshot=0, scene="esp.tcl", bg="white", suffix="",
            verbose=True):
    """One VMD run for the given views. Returns (rc, renderer).

    A VMD crash leaves no error text behind, only missing files - which is why
    the complete log is written out and appended, so that the run before it
    can still be read as well.

    """
    tcl = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "render_esp.tcl").replace("\\", "/")
    with open("_render_opts.tcl", "w", encoding="utf-8") as fh:
        fh.write(f"set ESP_RES {{{w} {h}}}\n"
                 f"set ESP_QUIT 1\n"
                 f"set ESP_AO {ao}\n"
                 "set ESP_BG {" + bg + "}\n"
                 "set ESP_SUFFIX {" + suffix + "}\n"
                 f"set ESP_OPAQUE {opaque}\n"
                 f"set ESP_SNAPSHOT {snapshot}\n"
                 f"set ESP_VIEWS {{{' '.join(views)}}}\n"
                 # Scene file: normally esp.tcl, in the self test
                 # esp_check.tcl - that leaves the committed scene untouched.
                 "set ESP_SCENE {" + scene + "}\n"
                 # The target folder MUST be passed on: otherwise the Tcl
                 # script writes to images/ while this side looks in outdir/ -
                 # the images are produced but count as failed.
                 # Forward slashes, so that Tcl does not read Windows
                 # backslashes as escapes.
                 "set ESP_OUTDIR {" + outdir.replace("\\", "/") + "}\n"
                 "set ESP_PREFIX {" + prefix + "}\n"
                 # Braces: Tcl substitutes nothing inside them, and spaces in
                 # the path ("2. Semester") do no harm.
                 "source {" + tcl + "}\n")
    renderer = "TachyonInternal"
    try:
        out = subprocess.run([vmd, "-e", "_render_opts.tcl"],
                             capture_output=True, text=True, timeout=1800)
        text = (out.stdout or "") + (out.stderr or "")
        with open(os.path.join(outdir, "_vmd.log"), "a",
                  encoding="utf-8", errors="replace") as fh:
            fh.write(f"\n===== Pass: {label} ({', '.join(views)}) =====\n")
            fh.write(text)
        for line in text.splitlines():
            if line.startswith(("->", "!", "==", "Renderer:", "Window:",
                                "Target:")):
                if verbose or line.startswith("!"):
                    print("   ", line)
                if line.startswith("Renderer:"):
                    renderer = line.split(":", 1)[1].strip()
        return out.returncode, renderer
    except subprocess.TimeoutExpired:
        sys.exit("VMD did not respond within 30 minutes.")
    finally:
        if os.path.exists("_render_opts.tcl"):
            os.remove("_render_opts.tcl")


def read_scene(scene="esp.tcl"):
    """Isovalue, colour scale and the statistics from the scene resp. tp.cube."""
    iso, rng = 0.001, None
    if os.path.exists(scene):
        text = open(scene, encoding="utf-8", errors="replace").read()
        m = re.search(r"^set ISO\s+(\S+)", text, re.M)
        if m:
            iso = float(m.group(1))
        m = re.search(r"^set RANGE\s+(\S+)", text, re.M)
        if m:
            rng = float(m.group(1))
    stats = None
    if os.path.exists("tp.cube"):
        m = STATS_RE.search(open("tp.cube", encoding="utf-8",
                                 errors="replace").readline())
        if m:
            stats = (float(m.group(1)), float(m.group(2)), float(m.group(3)))
    return iso, rng, stats


def tga_to_png(outdir, prefix, keep_tga=False, suffix=""):
    """VMD's Tachyon writes TGA - convert to PNG for storage."""
    from PIL import Image
    done = []
    for view in VIEWS:
        tga = os.path.join(outdir, f"{prefix}_{view}{suffix}.tga")
        if not os.path.exists(tga):
            continue
        png = os.path.join(outdir, f"{prefix}_{view}{suffix}.png")
        try:
            with Image.open(tga) as im:
                im.convert("RGB").save(png)
        except Exception as err:
            # A VMD that crashed mid-write leaves half a file behind. That
            # should not take the rest of the run down with it.
            print(f"    ! {tga} unreadable ({err}) - skipped",
                  file=sys.stderr)
            continue
        if not keep_tga:
            os.remove(tga)
        with Image.open(png) as im:
            done.append((png, im.size))
    return done


def colorbar(path, rng, dpi=300, rainbow=False):
    """Horizontal colour bar as a PNG of its own.

    Red-white-blue in the same order as VMD's RWB colour scale, so that the
    bar matches the images - and matches the PyMOL variant, which uses the
    same ramp.

    The rainbow bar shows the five anchor colours red - yellow - green - cyan
    - blue, the same as in the PyMOL pipeline. VMD's built-in scales know only
    three; the scene therefore reprograms the colour table itself (esp_ramp in
    esp_template.tcl), so that bar and image agree.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.colorbar import ColorbarBase

    # Exactly the anchor colours of the PyMOL pipeline (render_esp.py,
    # RAMP_HEX) - and the same ones esp_ramp writes into the 1024 entries of
    # VMD's colour table. Bar and image therefore show the same ramp, and each
    # project shows the same one as the other.
    RAINBOW = ["#d40000", "#f0e000", "#00a000", "#00c8d4", "#0030d4"]
    cmap = LinearSegmentedColormap.from_list(
        "esp", RAINBOW if rainbow else ["red", "white", "blue"])
    fig = plt.figure(figsize=(6.0, 1.0))
    ax = fig.add_axes([0.06, 0.42, 0.88, 0.30])
    cb = ColorbarBase(ax, cmap=cmap, norm=Normalize(-rng, rng),
                      orientation="horizontal")
    cb.set_label("ESP / a.u.", fontsize=11)
    cb.set_ticks([-rng, -rng / 2, 0.0, rng / 2, rng])
    cb.ax.tick_params(labelsize=10)
    fig.savefig(path, dpi=dpi, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def settings(path, prefix, iso, rng, stats, size, renderer, made=None,
             ao=False, rainbow=False, backgrounds=("white",)):
    lines = [
        "Render parameters (written by render_espVMD.py)",
        "=" * 55,
        f"Molecule          : {prefix}",
        f"Density cube      : td.cube",
        f"ESP cube          : tp.cube",
        f"Isovalue rho      : {iso:g} a.u.",
    ]
    if stats:
        # 3-tuple from read_scene, 4-tuple (with the point count) from the
        # conversion
        vmin, vmax = stats[0], stats[1]
        npts = stats[3] if len(stats) > 3 else None
        lines += [
            f"V_S,min           : {vmin:+.5f} a.u. "
            f"({vmin * HARTREE_TO_KJ:+.1f} kJ/(mol*e))",
            f"V_S,max           : {vmax:+.5f} a.u. "
            f"({vmax * HARTREE_TO_KJ:+.1f} kJ/(mol*e))",
            (f"Shell points      : {npts}" if npts else ""),
            "                    (grid points in the band rho = iso +/- 12 %,",
            "                     the same calculation as in the PyMOL",
            "                     pipeline and the same numbers. The sigma",
            "                     hole is determined only there, not here.)",
        ]
    lines += [
        f"Colour scale      : {-rng:+.4f} .. {rng:+.4f} a.u.",
        (f"Colour ramp       : RGB (rainbow: red negative, green zero, "
         f"blue positive)" if rainbow else
         f"Colour ramp       : RWB (red negative, white zero, blue positive)"),
        f"Background        : {', '.join(backgrounds)}",
        f"Ambient occl.     : {'on' if ao else 'off'}",
        f"Image size        : {size[0]} x {size[1]} px" if size
        else "Image size        : unknown",
        f"Renderer          : {renderer}",
        f"Projection        : orthoscopic",
    ]
    if made:
        for view in VIEWS:
            lines.append(f"View {view:<13}: {made.get(view, 'not rendered')}")
    else:
        lines.append("Views             : pi, edge, sigma")
    lines += [
        f"Written           : {datetime.datetime.now():%Y-%m-%d %H:%M}",
    ]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(ln for ln in lines if ln) + "\n")


def render_all(outdir="images", prefix=None, iso=None, rng=None, stats=None,
               vmd=None, res="1600x1280", ao=False, backgrounds=None,
               keep_tga=False, dpi=300, no_vmd=False, scene=None,
               rainbow=False, verbose=True):
    """Image set for the CURRENT directory. Expects the scene and the cubes.

    Returns a dict with made (view -> what it was rendered with), renderer,
    size, rng, iso, stats and the paths. Used by main() and by run_allVMD.py -
    the retry logic should exist only once.
    """
    prefix = prefix or os.path.basename(os.path.abspath(os.getcwd()))
    # Its own name suffix and its own scene, otherwise a rainbow run
    # overwrites the red-white-blue image set of the same molecule.
    # prefix stays the molecule name, out_prefix names the files.
    scene = scene or ("esp_rainbow.tcl" if rainbow else "esp.tcl")
    out_prefix = f"{prefix}_rainbow" if rainbow else prefix
    os.makedirs(outdir, exist_ok=True)
    scene_iso, scene_rng, scene_stats = read_scene(scene)
    iso = scene_iso if iso is None else iso
    rng = rng if rng is not None else (scene_rng if scene_rng is not None
                                       else (scene_stats[2] if scene_stats
                                             else 0.035))
    stats = stats if stats is not None else scene_stats
    renderer = "TachyonInternal"

    # Ambient occlusion is OFF by default, see --ao. If it stays on, it is
    # the first pass.
    #
    # In the axial view Tachyon likes to die on the combination of many
    # transparent layers and ambient occlusion. Rather than turning everything
    # down at once: full quality first, then only the missing views with less.
    # What was produced by what is recorded in settings.txt.
    passes = ([("AO + transparency", dict(ao=1, opaque=0, snapshot=0))]
              if ao else []) + [
        ("transparency", dict(ao=0, opaque=0, snapshot=0)),
        # Tachyon bails out in the axial view on the many transparent layers
        # crossed. The window capture knows no recursion and keeps the
        # transparency - only after that, opaque.
        ("transparency, window capture", dict(ao=0, opaque=0, snapshot=1)),
        ("opaque", dict(ao=0, opaque=1, snapshot=0))]

    # Every background colour is an image set of its own with its own retry
    # logic: that Tachyon gets through on white says nothing about whether it
    # does on black too. Only with more than one colour do the files get a
    # name suffix, otherwise the standard set would suddenly be called
    # <molecule>_pi_white.png.
    backgrounds = list(backgrounds) if backgrounds else ["white"]
    tag = {bg: (f"_{bg}" if len(backgrounds) > 1 else "") for bg in backgrounds}
    made_per_bg, done = {}, []

    if not no_vmd:
        vmd = find_vmd(vmd)
        if not vmd:
            sys.exit("vmd not found. Pass it with --vmd <path>, or add the "
                     "VMD folder to the PATH.")
        try:
            w, h = res.lower().split("x")
        except ValueError:
            sys.exit("--res expects e.g. 1600x1280")
        if verbose:
            print(f"[1] VMD: {vmd}")

        for bg in backgrounds:
            made = {}
            made_per_bg[bg] = made
            if len(backgrounds) > 1 and verbose:
                print(f"    background {bg}")
            for label, opt in passes:
                todo = [v for v in VIEWS if v not in made]
                if not todo:
                    break
                if label != passes[0][0] and verbose:
                    print(f"    second attempt for {', '.join(todo)}: {label}")
                # Time stamp before the run: a TGA from an earlier, crashed
                # pass would otherwise still be lying there and would be
                # counted as a success.
                t0 = time.time() - 1.0
                rc, used = run_vmd(vmd, outdir, out_prefix, todo, w, h,
                                       opt["ao"], opt["opaque"], label,
                                       snapshot=opt["snapshot"], scene=scene,
                                       bg=bg, suffix=tag[bg], verbose=verbose)
                # The renderer of the FIRST pass is what the log records;
                # what a caught-up view was produced with is told by its own
                # line.
                if bg == backgrounds[0] and label == passes[0][0]:
                    renderer = used
                for v in todo:
                    tga = os.path.join(outdir,
                                       f"{out_prefix}_{v}{tag[bg]}.tga")
                    if os.path.exists(tga) and os.path.getmtime(tga) >= t0:
                        made[v] = label
                missing = [v for v in VIEWS if v not in made]
                if missing and label == passes[-1][0]:
                    print(f"    ! not rendered: {', '.join(missing)} "
                          f"({bg}, VMD return code {rc}). Log: "
                          f"{os.path.join(outdir, '_vmd.log')}",
                          file=sys.stderr)

    if verbose:
        print("[2] TGA -> PNG")
    for bg in backgrounds:
        done += tga_to_png(outdir, out_prefix, keep_tga=keep_tga,
                           suffix=tag[bg])
    for png, size in done:
        if verbose:
            print(f"    -> {png}  ({size[0]} x {size[1]} px)")
    if not done and not no_vmd:
        print("    ! no TGA files found - did VMD render?",
              file=sys.stderr)

    # For the summary the first background counts; the line in settings.txt
    # names all of them.
    made = made_per_bg.get(backgrounds[0], {})

    cb = os.path.join(outdir, f"{out_prefix}_colorbar.png")
    colorbar(cb, rng, dpi=dpi, rainbow=rainbow)
    if verbose:
        print(f"[3] colour bar -> {cb}  (+/- {rng:.4f} a.u.)")

    st = os.path.join(outdir, f"{out_prefix}_settings.txt")
    size = done[0][1] if done else None
    settings(st, prefix, iso, rng, stats, size, renderer, made,
             ao=ao, rainbow=rainbow, backgrounds=backgrounds)
    if verbose:
        print(f"[4] -> {st}")

    return {"prefix": prefix, "made": made, "renderer": renderer, "size": size,
            "rng": rng, "iso": iso, "stats": stats, "outdir": outdir,
            "images": [p for p, _ in done], "colorbar": cb, "settings": st}


def main(argv=None):
    p = argparse.ArgumentParser(
        description="The three standard views and the colour bar into images/")
    p.add_argument("--no-vmd", action="store_true",
                   help="do not call VMD, only convert TGA and build the bar")
    p.add_argument("--vmd", help="path to vmd.exe, if it is not on the PATH")
    p.add_argument("--res", default="1600x1280", help="image size (window)")
    p.add_argument("--ao", action="store_true",
                   help="ambient occlusion on. Off by default: it lays soft "
                        "shadows into the hollows, but also lines the sticks "
                        "as grey doubles on the isosurface. PyMOL renders "
                        "without it too.")
    p.add_argument("--backgrounds", nargs="+", default=["white"],
                   metavar="COLOUR",
                   help="background colours, e.g. white black. From two "
                        "colours on, the files get a suffix: "
                        "<molecule>_pi_black.png")
    p.add_argument("--keep-tga", action="store_true")
    p.add_argument("--dpi", type=int, default=300, help="resolution of the bar")
    p.add_argument("--outdir", default="images")
    p.add_argument("--rainbow", action="store_true",
                   help="rainbow ramp: renders esp_rainbow.tcl and writes "
                        "<prefix>_rainbow_*, the standard set is kept")
    p.add_argument("--scene", default=None,
                   help="scene file in the molecule folder (default esp.tcl, "
                        "with --rainbow esp_rainbow.tcl; the self test uses "
                        "esp_check.tcl)")
    args = p.parse_args(argv)
    if args.scene is None:
        args.scene = "esp_rainbow.tcl" if args.rainbow else "esp.tcl"

    if not os.path.exists(args.scene):
        # The scene comes from the first step, not from here. If the cubes
        # are already there, --tcl-only brings it back in a fraction of a
        # second.
        here = os.path.dirname(os.path.abspath(__file__))
        conv = os.path.join(here, "xyzToCubeToVMDVis.py")
        if os.path.exists("td.cube") and os.path.exists("tp.cube"):
            sys.exit(f"{args.scene} is missing, but the cubes are there. "
                     f"Write it again with:\n"
                     f"    python {conv} td.cube tp.cube --tcl-only")
        sys.exit(f"{args.scene} not found - start from the molecule folder.\n"
                 f"Convert first:\n"
                 f"    python {conv} --struct <structure> td.xyz tp.xyz")

    render_all(outdir=args.outdir, vmd=args.vmd, res=args.res,
               ao=args.ao, backgrounds=args.backgrounds, scene=args.scene,
               rainbow=args.rainbow,
               keep_tga=args.keep_tga, dpi=args.dpi, no_vmd=args.no_vmd)
    return 0


if __name__ == "__main__":
    sys.exit(main())

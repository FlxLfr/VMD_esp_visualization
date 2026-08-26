#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
render_esp.py - Standardbildersatz aus einer fertigen VMD-Szene

    python render_esp.py            # rendert, konvertiert, Farbskala
    python render_esp.py --no-vmd   # nur TGA -> PNG und Farbskala

Aus dem Molekuelordner starten, dort wo esp.tcl, td.cube und tp.cube liegen.
Ergebnis in images/:

    <molekuel>_pi.png  _edge.png  _sigma.png  _colorbar.png  _settings.txt

Derselbe Dateisatz wie in der PyMOL-Pipeline, mit Absicht: nur so lassen sich
die beiden Bildersaetze nebeneinanderlegen und wirklich vergleichen.

Aufgabenteilung: render_esp.tcl macht den VMD-Teil (Ansichten, Tachyon), dieses
Skript ruft es auf, wandelt VMDs TGA nach PNG und erzeugt die Farbskala mit
matplotlib. Die Skala kommt bewusst nicht aus VMD - PyMOL kann das genausowenig
und macht sie in render_esp.py ebenfalls mit matplotlib.

Braucht matplotlib und Pillow.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import datetime

VIEWS = ("pi", "edge", "sigma")
STATS_RE = re.compile(r"V_S,min=(\S+)\s+V_S,max=(\S+)\s+range=(\S+)\s+iso=(\S+)")
HARTREE_TO_KJ = 2625.4996


def find_vmd(explicit=None):
    """vmd.exe finden: Argument, PATH, VMDDIR, dann die ueblichen Pfade.

    VMDDIR setzt der Windows-Installer auf das Installationsverzeichnis - das
    ist der zuverlaessigste Treffer, weil es unabhaengig davon stimmt, wohin
    installiert wurde. Der Installer erweitert den PATH dagegen NICHT.
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
            shadows=0, snapshot=0, scene="esp.tcl", verbose=True):
    """Ein VMD-Lauf fuer die angegebenen Ansichten. Rueckgabe: (rc, renderer).

    Ein Absturz von VMD hinterlaesst keinen Fehlertext, nur fehlende Dateien -
    deshalb wird das komplette Protokoll mitgeschrieben und angehaengt, damit
    auch der Lauf davor noch nachlesbar ist.
    """
    tcl = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "render_esp.tcl").replace("\\", "/")
    with open("_render_opts.tcl", "w", encoding="utf-8") as fh:
        fh.write(f"set ESP_RES {{{w} {h}}}\n"
                 f"set ESP_QUIT 1\n"
                 f"set ESP_AO {ao}\n"
                 f"set ESP_SHADOWS {shadows}\n"
                 f"set ESP_OPAQUE {opaque}\n"
                 f"set ESP_SNAPSHOT {snapshot}\n"
                 f"set ESP_VIEWS {{{' '.join(views)}}}\n"
                 # Szenendatei: im Normalfall esp.tcl, im Selbsttest
                 # esp_check.tcl - so bleibt die committete Szene unberuehrt.
                 "set ESP_SCENE {" + scene + "}\n"
                 # Geschweifte Klammern: Tcl substituiert darin nichts und
                 # Leerzeichen im Pfad ("2. Semester") stoeren nicht.
                 "source {" + tcl + "}\n")
    renderer = "TachyonInternal"
    try:
        out = subprocess.run([vmd, "-e", "_render_opts.tcl"],
                             capture_output=True, text=True, timeout=1800)
        text = (out.stdout or "") + (out.stderr or "")
        with open(os.path.join(outdir, "_vmd.log"), "a",
                  encoding="utf-8", errors="replace") as fh:
            fh.write(f"\n===== Durchlauf: {label} ({', '.join(views)}) =====\n")
            fh.write(text)
        for line in text.splitlines():
            if line.startswith(("->", "!", "==", "Renderer:", "Fenster:")):
                if verbose or line.startswith("!"):
                    print("   ", line)
                if line.startswith("Renderer:"):
                    renderer = line.split(":", 1)[1].strip()
        return out.returncode, renderer
    except subprocess.TimeoutExpired:
        sys.exit("VMD hat nicht innerhalb von 30 Minuten geantwortet.")
    finally:
        if os.path.exists("_render_opts.tcl"):
            os.remove("_render_opts.tcl")


def read_scene(scene="esp.tcl"):
    """Isowert, Farbskala und die Kennzahlen aus der Szene bzw. tp.cube."""
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


def tga_to_png(outdir, prefix, keep_tga=False):
    """VMDs Tachyon schreibt TGA - fuer die Ablage nach PNG wandeln."""
    from PIL import Image
    done = []
    for view in VIEWS:
        tga = os.path.join(outdir, f"{prefix}_{view}.tga")
        if not os.path.exists(tga):
            continue
        png = os.path.join(outdir, f"{prefix}_{view}.png")
        try:
            with Image.open(tga) as im:
                im.convert("RGB").save(png)
        except Exception as err:
            # Ein mitten im Schreiben abgestuerztes VMD hinterlaesst eine halbe
            # Datei. Die soll den Rest des Durchlaufs nicht mitnehmen.
            print(f"    ! {tga} unlesbar ({err}) - uebersprungen",
                  file=sys.stderr)
            continue
        if not keep_tga:
            os.remove(tga)
        with Image.open(png) as im:
            done.append((png, im.size))
    return done


def colorbar(path, rng, dpi=300):
    """Waagerechte Farbskala als eigenes PNG.

    Rot-weiss-blau in derselben Reihenfolge wie VMDs Farbskala RWB, damit der
    Balken zu den Bildern passt - und zur PyMOL-Variante, die dieselbe Rampe
    benutzt.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.colorbar import ColorbarBase

    cmap = LinearSegmentedColormap.from_list("rwb", ["red", "white", "blue"])
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
             shadows=False, ao=False):
    lines = [
        "Renderparameter (erzeugt von render_espVMD.py)",
        "=" * 55,
        f"Molekuel          : {prefix}",
        f"Dichte-Cube       : td.cube",
        f"ESP-Cube          : tp.cube",
        f"Isowert rho       : {iso:g} a.u.",
    ]
    if stats:
        # 3er-Tupel aus read_scene, 4er (mit Punktzahl) aus der Konvertierung
        vmin, vmax = stats[0], stats[1]
        npts = stats[3] if len(stats) > 3 else None
        lines += [
            f"V_S,min           : {vmin:+.5f} a.u. "
            f"({vmin * HARTREE_TO_KJ:+.1f} kJ/(mol*e))",
            f"V_S,max           : {vmax:+.5f} a.u. "
            f"({vmax * HARTREE_TO_KJ:+.1f} kJ/(mol*e))",
            (f"Schalenpunkte     : {npts}" if npts else ""),
            "                    (Gitterpunkte auf der Schale, nicht "
            "interpoliert -",
            "                     zum Zitieren die Werte der PyMOL-Pipeline "
            "nehmen)",
        ]
    lines += [
        f"Farbskala         : {-rng:+.4f} .. {rng:+.4f} a.u.",
        f"Farbrampe         : RWB (rot negativ, blau positiv)",
        f"Schlagschatten    : {'an' if shadows else 'aus'}",
        f"Umgebungsverdeck. : {'an' if ao else 'aus'}",
        f"Bildgroesse       : {size[0]} x {size[1]} px" if size
        else "Bildgroesse       : unbekannt",
        f"Renderer          : {renderer}",
        f"Projektion        : orthoskopisch",
    ]
    if made:
        for view in VIEWS:
            lines.append(f"Ansicht {view:<10}: {made.get(view, 'nicht gerendert')}")
    else:
        lines.append("Ansichten         : pi, edge, sigma")
    lines += [
        f"Erzeugt           : {datetime.datetime.now():%Y-%m-%d %H:%M}",
    ]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(ln for ln in lines if ln) + "\n")


def render_all(outdir="images", prefix=None, iso=None, rng=None, stats=None,
               vmd=None, res="1600x1280", ao=False, shadows=False,
               keep_tga=False, dpi=300, no_vmd=False, scene="esp.tcl",
               verbose=True):
    """Bildersatz fuer das AKTUELLE Verzeichnis. Erwartet die Szene und die Cubes.

    Rueckgabe: dict mit made (Ansicht -> womit gerendert), renderer, size, rng,
    iso, stats und den Pfaden. Von main() und von run_allVMD.py benutzt - die
    Wiederholungslogik soll es nur einmal geben.
    """
    prefix = prefix or os.path.basename(os.path.abspath(os.getcwd()))
    os.makedirs(outdir, exist_ok=True)
    scene_iso, scene_rng, scene_stats = read_scene(scene)
    iso = scene_iso if iso is None else iso
    rng = rng if rng is not None else (scene_rng if scene_rng is not None
                                       else (scene_stats[2] if scene_stats
                                             else 0.035))
    stats = stats if stats is not None else scene_stats
    renderer = "TachyonInternal"

    # Umgebungsverdeckung ist standardmaessig AUS, siehe --ao. Bleibt sie an,
    # ist sie der erste Durchgang.
    #
    # Tachyon stirbt in der Achsenansicht gern an der Kombination aus vielen
    # transparenten Lagen und Umgebungsverdeckung. Statt gleich alles
    # herunterzuschrauben: erst in voller Qualitaet, dann nur die fehlenden
    # Ansichten mit weniger. Was womit entstanden ist, steht in settings.txt.
    passes = ([("AO + Transparenz", dict(ao=1, opaque=0, snapshot=0))]
              if ao else []) + [
        ("Transparenz", dict(ao=0, opaque=0, snapshot=0)),
        # Tachyon steigt in der Achsenansicht an den vielen durchquerten
        # transparenten Lagen aus. Der Fenstermitschnitt kennt keine
        # Rekursion und behaelt die Transparenz - erst danach opak.
        ("Transparenz, Fensterbild", dict(ao=0, opaque=0, snapshot=1)),
        ("opak", dict(ao=0, opaque=1, snapshot=0))]
    made = {}

    if not no_vmd:
        vmd = find_vmd(vmd)
        if not vmd:
            sys.exit("vmd nicht gefunden. Mit --vmd <pfad> angeben, oder den "
                     "VMD-Ordner in den PATH aufnehmen (README Abschnitt 8).")
        try:
            w, h = res.lower().split("x")
        except ValueError:
            sys.exit("--res erwartet z.B. 1600x1280")
        if verbose:
            print(f"[1] VMD: {vmd}")

        for label, opt in passes:
            todo = [v for v in VIEWS if v not in made]
            if not todo:
                break
            if label != passes[0][0] and verbose:
                print(f"    Zweiter Anlauf fuer {', '.join(todo)}: {label}")
            rc, used = run_vmd(vmd, outdir, prefix, todo, w, h,
                                   opt["ao"], opt["opaque"], label,
                                   shadows=1 if shadows else 0,
                                   snapshot=opt["snapshot"], scene=scene,
                                   verbose=verbose)
            # Der Renderer des ERSTEN Durchgangs steht im Protokoll; womit eine
            # nachgeholte Ansicht entstanden ist, sagt ihre eigene Zeile.
            if label == passes[0][0]:
                renderer = used
            for v in todo:
                if os.path.exists(os.path.join(outdir, f"{prefix}_{v}.tga")):
                    made[v] = label
            missing = [v for v in VIEWS if v not in made]
            if missing and label == passes[-1][0]:
                print(f"    ! nicht gerendert: {', '.join(missing)} "
                      f"(VMD-Rueckgabewert {rc}). Protokoll: "
                      f"{os.path.join(outdir, '_vmd.log')}", file=sys.stderr)

    if verbose:
        print("[2] TGA -> PNG")
    done = tga_to_png(outdir, prefix, keep_tga=keep_tga)
    for png, size in done:
        if verbose:
            print(f"    -> {png}  ({size[0]} x {size[1]} px)")
    if not done and not no_vmd:
        print("    ! keine TGA-Dateien gefunden - hat VMD gerendert?",
              file=sys.stderr)

    cb = os.path.join(outdir, f"{prefix}_colorbar.png")
    colorbar(cb, rng, dpi=dpi)
    if verbose:
        print(f"[3] Farbskala -> {cb}  (+/- {rng:.4f} a.u.)")

    st = os.path.join(outdir, f"{prefix}_settings.txt")
    size = done[0][1] if done else None
    settings(st, prefix, iso, rng, stats, size, renderer, made,
             shadows=shadows, ao=ao)
    if verbose:
        print(f"[4] -> {st}")

    return {"prefix": prefix, "made": made, "renderer": renderer, "size": size,
            "rng": rng, "iso": iso, "stats": stats, "outdir": outdir,
            "images": [p for p, _ in done], "colorbar": cb, "settings": st}


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Die drei Standardansichten und die Farbskala nach images/")
    p.add_argument("--no-vmd", action="store_true",
                   help="VMD nicht aufrufen, nur TGA wandeln und Skala bauen")
    p.add_argument("--vmd", help="Pfad zu vmd.exe, falls nicht im PATH")
    p.add_argument("--res", default="1600x1280", help="Bildgroesse (Fenster)")
    p.add_argument("--ao", action="store_true",
                   help="Umgebungsverdeckung an. Standard aus: sie legt weiche "
                        "Schatten in die Vertiefungen, saeumt dabei aber auch "
                        "die Staebchen als graue Doppelgaenger auf der "
                        "Isoflaeche. PyMOL rendert ebenfalls ohne.")
    p.add_argument("--shadows", action="store_true",
                   help="Schlagschatten an. Standard aus: sie werfen die "
                        "Staebchen als graue Kapseln auf die Isoflaeche, was "
                        "wie ein Datenartefakt aussieht.")
    p.add_argument("--keep-tga", action="store_true")
    p.add_argument("--dpi", type=int, default=300, help="Aufloesung der Skala")
    p.add_argument("--outdir", default="images")
    p.add_argument("--scene", default="esp.tcl",
                   help="Szenendatei im Molekuelordner (Standard esp.tcl; "
                        "der Selbsttest benutzt esp_check.tcl)")
    args = p.parse_args(argv)

    if not os.path.exists(args.scene):
        # Die Szene kommt aus dem ersten Schritt, nicht von hier. Wenn die
        # Cubes schon da sind, ist sie mit --tcl-only in Sekundenbruchteilen
        # zurueck.
        here = os.path.dirname(os.path.abspath(__file__))
        conv = os.path.join(here, "xyzToCubeToVMDVis.py")
        if os.path.exists("td.cube") and os.path.exists("tp.cube"):
            sys.exit(f"{args.scene} fehlt, die Cubes sind aber da. Neu "
                     f"schreiben mit:\n"
                     f"    python {conv} td.cube tp.cube --tcl-only")
        sys.exit(f"{args.scene} nicht gefunden - aus dem Molekuelordner "
                 f"starten.\nErst konvertieren:\n"
                 f"    python {conv} --struct <struktur> td.xyz tp.xyz")

    render_all(outdir=args.outdir, vmd=args.vmd, res=args.res,
               ao=args.ao, shadows=args.shadows, scene=args.scene,
               keep_tga=args.keep_tga, dpi=args.dpi, no_vmd=args.no_vmd)
    return 0


if __name__ == "__main__":
    sys.exit(main())

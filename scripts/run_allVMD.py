#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_allVMD.py - Stapellauf der VMD-Pipeline

    python run_allVMD.py                 # Selbsttest auf reference/
    python run_allVMD.py --root ../sandbox

Ohne --root laeuft das Skript auf dem mitgelieferten Referenzdatensatz. Das ist
der Selbsttest: er schreibt nach images_check/ und summary_check.csv, damit die
committeten Referenzbilder unangetastet bleiben, und man vergleicht danach die
eigenen Zahlen mit denen in reference/.

Durchsucht einen Verzeichnisbaum und macht fuer jeden Molekuelordner alles in
einem Rutsch: pointval -> Cube, esp.tcl schreiben, die drei Ansichten mit
Tachyon rendern, Farbskala und settings.txt dazu. Am Ende steht eine
summary.csv mit den Kennzahlen und den Renderparametern.

Ein Molekuelordner ist jedes Verzeichnis mit

    td.xyz  + tp.xyz     (Turbomole-Ausgabe, wird konvertiert)
oder td.cube + tp.cube   (schon konvertiert)

plus einer Strukturdatei (.mol, .sdf oder .xyz).

Ablauf: erst analysieren, dann EINMAL rendern
--------------------------------------------
V_S,min und V_S,max stehen in den Cube-Dateien, dafuer muss nichts gerendert
werden. Das Skript sammelt also zuerst die Kennzahlen aller Molekuele und
rendert danach genau einen Durchgang.

    --esp-range auto     jedes Molekuel auf seiner eigenen Skala (Standard)
    --esp-range 0.035    fester Wert
    --esp-range common   gemeinsame Skala, sofort aus diesem Lauf

Vergleichbar sind die Bilder nur mit EINER Skala fuer alle. Dafuer rendert die
PyMOL-Pipeline zweimal (--two-pass) - hier nicht: am Ende jedes Laufs steht der
kleinste Wert, der alle Molekuele abdeckt, samt fertiger Befehlszeile in der
Konsole. Man liest ihn, entscheidet selbst, welche Skala das Bild zeigen soll,
und startet den Lauf, der zaehlt. Ein Bild kostet in VMD ein bis zwei Minuten;
ein automatischer erster Satz waere Ausschuss.

Fuer den direkten Vergleich mit der PyMOL-Pipeline den Wert aus deren
summary.csv uebernehmen.

Nicht enthalten
---------------
Kein sigma-Loch, keine trilineare Interpolation auf die Isoflaeche - das steht
im Schwesterprojekt und wird hier nicht ein zweites Mal implementiert. Die
Kennzahlen hier stammen von Gitterpunkten nahe der Schale und dienen der
Farbskala; zitiert werden die Werte der PyMOL-Pipeline.
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import xyzToCubeToVMDVis as conv          # noqa: E402
import render_espVMD as rend              # noqa: E402
from constants import HARTREE_TO_KJ       # noqa: E402

STRUCT_EXT = (".mol", ".sdf", ".xyz")
GRID_NAMES = ("td", "tp")
# reference/ liegt neben scripts/, nicht im aktuellen Arbeitsverzeichnis.
# Ohne --root ist der Lauf damit der Selbsttest, egal von wo aufgerufen.
DEFAULT_ROOT = os.path.normpath(os.path.join(_HERE, "..", "reference"))


# ----------------------------------------------------------------------------
# Molekuelordner finden
# ----------------------------------------------------------------------------

def find_structure(folder, exclude):
    """Erste Strukturdatei im Ordner, die keine Gitterdatei ist."""
    excl = {os.path.abspath(p) for p in exclude if p}
    # .mol/.sdf vor .xyz: sie tragen Bindungsinformation, und ein nacktes .xyz
    # kollidiert leicht mit den Turbomole-Gitterdateien derselben Endung.
    for ext in STRUCT_EXT:
        for name in sorted(os.listdir(folder)):
            if not name.lower().endswith(ext):
                continue
            path = os.path.join(folder, name)
            if os.path.abspath(path) in excl:
                continue
            if os.path.splitext(name)[0] in GRID_NAMES:
                continue                   # td.xyz / tp.xyz sind Daten
            return path
    return None


def discover(root):
    """Alle Molekuelordner unterhalb von ``root``, nach Namen sortiert."""
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
            print(f"  ! {dirpath}: Gitter gefunden, aber keine Strukturdatei "
                  f"({'/'.join(STRUCT_EXT)}) - uebersprungen")
            continue
        found.append({"dir": dirpath, "struct": struct,
                      "has_raw": has_raw, "has_cubes": has_cubes})
    return sorted(found, key=lambda e: e["dir"])


def cube_dims(path):
    """Gitterabmessungen aus dem Cube-Kopf, ohne die Daten zu lesen."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for _ in range(3):
            fh.readline()
        return "x".join(fh.readline().split()[0] for _ in range(3))


# ----------------------------------------------------------------------------
# Schritt 1: Cubes und Kennzahlen
# ----------------------------------------------------------------------------

def ensure_cubes(entry, args):
    """Cubes bereitstellen und (V_S,min, V_S,max, Skala, Punkte) zurueckgeben.

    Drei Wege, in dieser Reihenfolge:
      1. Cubes fehlen oder --force-convert: aus den pointval-Gittern erzeugen,
         dabei fallen die Kennzahlen ohnehin an.
      2. Cubes da und tragen die Kennzahlen im Kopf: direkt uebernehmen.
      3. Cubes da, aber ohne Kennzahlen (aeltere Laeufe): Cubes lesen und
         nachrechnen - rund 15 s statt Minuten fuer eine Neukonvertierung.
    """
    folder = entry["dir"]
    cubes = {t: os.path.join(folder, f"{t}.cube") for t in GRID_NAMES}
    have = all(os.path.exists(c) for c in cubes.values())

    if have and not args.force_convert:
        stats = conv.read_stats(cubes["tp"])
        if stats:
            print("    Cubes vorhanden, Kennzahlen aus dem Cube-Kopf")
            return cubes, stats
        print("    Cubes vorhanden, aber ohne Kennzahlen - werden nachgerechnet")
        density = conv.read_cube(cubes["td"])
        esp = conv.read_cube(cubes["tp"])
        return cubes, conv.shell_range(density, esp, args.iso)

    if entry["struct"] is None:
        raise SystemExit(f"{folder}: Strukturdatei fehlt, ohne sie laesst sich "
                         f"kein Cube-Header schreiben.")
    atoms = conv.read_structure(entry["struct"], unit=args.struct_unit)
    grids = {}
    for tag in GRID_NAMES:
        raw = os.path.join(folder, f"{tag}.xyz")
        if not os.path.exists(raw):
            raise SystemExit(f"{folder}: weder {tag}.cube noch {tag}.xyz")
        print(f"    konvertiere {tag}.xyz -> {tag}.cube (stride {args.stride})")
        grids[tag] = conv.read_values(raw, verbose=False)

    stats = conv.shell_range(grids["td"][1], grids["tp"][1], args.iso)
    for tag in GRID_NAMES:
        info, data = grids[tag]
        comment = f"{info['quantity'] or tag} - aus {tag}.xyz"
        if tag == "tp" and stats:
            comment += conv.stats_comment(stats, args.iso)
        conv.write_cube(cubes[tag], info, data, atoms, args.stride, comment)
    return cubes, stats


# ----------------------------------------------------------------------------
# summary.csv
# ----------------------------------------------------------------------------

FIELDS = ["molecule", "structure", "grid", "iso_au", "shell_points",
          "VS_min_au", "VS_max_au", "VS_min_kJ", "VS_max_kJ",
          "esp_range_used_au", "esp_range_mode", "color_scale", "opacity",
          "resolution_px", "renderer", "ambient_occlusion", "backgrounds", "views"]


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
                "color_scale": r["color_scale"],
                "opacity": r["opacity"],
                "resolution_px": (f"{r['size'][0]}x{r['size'][1]}"
                                  if r.get("size") else ""),
                "renderer": r.get("renderer", ""),
                "ambient_occlusion": "on" if r["ao"] else "off",
                "backgrounds": " ".join(r["backgrounds"]),
                # Pro Ansicht, womit sie tatsaechlich entstanden ist - bei
                # einem Tachyon-Absturz greift der Fallback, und das darf im
                # Protokoll nicht verschwinden.
                "views": ";".join(f"{v}:{r['made'][v]}" for v in rend.VIEWS
                                  if v in r.get("made", {})),
            })
        if common is not None:
            fh.write(f"# gemeinsame Farbskala fuer alle Molekuele: "
                     f"+/- {common:.4f} a.u.\n")
        if advice is not None:
            # Auch in der Datei, nicht nur in der Konsole: wer die CSV spaeter
            # wieder aufmacht, sieht sonst nicht, dass die Bilder auf
            # verschiedenen Skalen liegen und was ein gemeinsamer Satz kostet.
            fh.write(f"# jedes Molekuel auf seiner eigenen Skala; ein "
                     f"vergleichbarer Satz braucht --esp-range "
                     f"{advice:.4f}\n")


# ----------------------------------------------------------------------------
# Empfehlung fuer den naechsten Lauf
# ----------------------------------------------------------------------------

def again_with(raw_argv, value):
    """Derselbe Aufruf noch einmal, nur mit festem --esp-range."""
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
            # Der Vorschlag zielt auf einen Bildersatz - ihn mit --no-render
            # anzubieten waere die eine Zeile, die garantiert nichts rendert.
            continue
        keep.append(f'"{a}"' if " " in a else a)
    return " ".join(["python run_allVMD.py"] + keep
                    + [f"--esp-range {value:.4f}"])


def advise(mode, common, needed, rows, raw_argv):
    """Sagen, welche Skala einen vergleichbaren Satz ergaebe - und wie man sie setzt.

    Statt eines zweiten Renderdurchgangs (so macht es die PyMOL-Pipeline mit
    --two-pass) steht hier nur der Vorschlag. Welche Skala vernuenftig ist,
    haengt davon ab, was das Bild zeigen soll: die kleinste, die nichts
    abschneidet, ist selten die aussagekraeftigste - ein Molekuel mit einem
    sehr negativen Sauerstoff drueckt sonst allen anderen die Kontraste flach.
    Diese Entscheidung gehoert dem Anwender, nicht dem Skript.
    """
    if needed is None:
        return
    used = sorted({r["rng"] for r in rows})
    print()
    print("-" * 70)
    if mode == "auto":
        if len(used) == 1:
            print(f"Alle Molekuele sind ohnehin auf +/- {used[0]:.4f} a.u. "
                  f"gelandet -")
            print("die Bilder sind also schon vergleichbar. Zum Festschreiben:")
        else:
            print("Jedes Molekuel hat seine eigene Skala - "
                  f"{', '.join(f'{v:.4f}' for v in used)} a.u.")
            print("Nebeneinanderlegen darf man die Bilder so NICHT.")
            print(f"Kleinster Wert, der alle abdeckt: +/- {needed:.4f} a.u.")
            print("Zum Rendern eines vergleichbaren Satzes:")
        print()
        print(f"    {again_with(raw_argv, needed)}")
    elif common is not None and common + 1e-12 < needed:
        low = [r["prefix"] for r in rows
               if r["stats"] and r["stats"][2] > common + 1e-12]
        print(f"! Die feste Skala +/- {common:.4f} a.u. schneidet ab.")
        print(f"  {', '.join(low)} braucht bis +/- {needed:.4f} a.u.; "
              f"alles darueber")
        print("  laeuft in die Saettigung und ist im Bild nicht mehr "
              "unterscheidbar.")
        print("  Absichtlich? Dann ist alles gut. Sonst:")
        print()
        print(f"    {again_with(raw_argv, needed)}")
    else:
        print(f"Skala +/- {common:.4f} a.u. deckt alle Molekuele ab "
              f"(noetig waeren +/- {needed:.4f}).")
    print("-" * 70)


# ----------------------------------------------------------------------------
# Hauptprogramm
# ----------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        description="Stapelkonvertierung und -rendering der VMD-ESP-Pipeline.")
    p.add_argument("--root", default=DEFAULT_ROOT,
                   help="Verzeichnisbaum mit den Molekuelordnern "
                        "(Standard: reference/ des Repositoriums - der "
                        "Selbsttest)")
    p.add_argument("--only", nargs="+", metavar="NAME",
                   help="nur diese Ordner, Platzhalter erlaubt: "
                        "--only paracetamol '*benzol'")

    g = p.add_argument_group("Konvertierung")
    g.add_argument("--stride", type=int, default=1,
                   help="jeden n-ten Gitterpunkt (2 => 8x kleinere Cubes). "
                        "Wirkt nur beim Erzeugen, nicht auf vorhandene Cubes.")
    g.add_argument("--struct-unit", choices=["angstrom", "bohr"],
                   default="angstrom")
    g.add_argument("--force-convert", action="store_true",
                   help="Cubes neu schreiben, auch wenn sie schon da sind")

    g = p.add_argument_group("Szene")
    g.add_argument("--esp-range", default="auto",
                   help="'auto' (Standard, jedes Molekuel auf seiner eigenen "
                        "Skala), ein fester Wert in a.u. - am Ende jedes Laufs "
                        "steht ein Vorschlag dafuer in der Konsole - oder "
                        "'common' (gemeinsame Skala, sofort aus diesem Lauf)")
    g.add_argument("--iso", type=float, default=0.001, help="Isowert in a.u.")
    g.add_argument("--opacity", type=float, default=0.50,
                   help="Deckkraft der Isoflaeche 0..1")
    g.add_argument("--scale", default="auto", help="Zoom, Zahl oder 'auto'")
    g.add_argument("--fill", type=float, default=0.85,
                   help="Anteil der Bildhoehe fuer das Molekuel bei --scale auto")
    g.add_argument("--color-scale", default=None,
                   help="VMD-Farbskala, Standard RWB (rot negativ); "
                        "mit --rainbow RGB")
    g.add_argument("--rainbow", action="store_true",
                   help="Regenbogenrampe statt rot-weiss-blau. Schreibt "
                        "esp_rainbow.tcl und einen eigenen Bildersatz "
                        "<molekuel>_rainbow_*, der Standardsatz bleibt "
                        "erhalten.")

    g = p.add_argument_group("Rendern")
    g.add_argument("--no-render", action="store_true",
                   help="nur konvertieren und esp.tcl schreiben")
    g.add_argument("--vmd", help="Pfad zu vmd.exe, falls nicht im PATH")
    g.add_argument("--res", default="1600x1280", help="Bildgroesse")
    g.add_argument("--ao", action="store_true",
                   help="Umgebungsverdeckung an. Standard aus: sie legt weiche "
                        "Schatten in die Vertiefungen, saeumt dabei aber auch "
                        "die Staebchen als graue Doppelgaenger auf der "
                        "Isoflaeche. PyMOL rendert ebenfalls ohne.")
    g.add_argument("--backgrounds", nargs="+", default=["white"],
                   metavar="FARBE",
                   help="Hintergrundfarben, z.B. white black. Ab zwei Farben "
                        "bekommen die Dateien einen Zusatz: "
                        "<molekuel>_pi_black.png")
    g.add_argument("--keep-tga", action="store_true")
    g.add_argument("--dpi", type=int, default=300, help="Aufloesung der Skala")
    g.add_argument("--images-dir", default=None,
                   help="Zielordner in jedem Molekuelordner (Standard: "
                        "'images', beim Selbsttest 'images_check')")
    g.add_argument("--summary", default=None,
                   help="Pfad der CSV (Standard: <root>/summary.csv)")
    args = p.parse_args(argv)
    # Fuer die Empfehlung am Ende: derselbe Aufruf, nur mit anderer Skala.
    raw_argv = list(sys.argv[1:]) if argv is None else list(argv)

    # Der Selbsttest darf nichts Committetes ueberschreiben: eigene Bildordner,
    # eigene Szenendatei, eigene CSV. Sonst laesst sich hinterher nicht mehr
    # sagen, ob die Referenz noch die Referenz ist.
    is_reference = os.path.abspath(args.root) == os.path.abspath(DEFAULT_ROOT)
    # Eigener Szenenname pro Rampe, sonst ueberschreibt ein Regenbogenlauf die
    # Szene des rot-weiss-blauen Laufs - die Bilder liegen aus demselben Grund
    # getrennt (<molekuel>_rainbow_*).
    scene_name = ("esp_check" if is_reference else "esp") \
        + ("_rainbow" if args.rainbow else "") + ".tcl"
    if args.images_dir is None:
        args.images_dir = "images_check" if is_reference else "images"
    # RGB ist VMDs Regenbogen: rot - gruen - blau, dieselbe Richtung wie die
    # Rampe der PyMOL-Pipeline. Ein ausdrueckliches --color-scale gewinnt.
    if args.color_scale is None:
        args.color_scale = "RGB" if args.rainbow else "RWB"

    # Vor dem ersten chdir absolut machen: das Rendern laeuft im Molekuelordner,
    # ein relativer Pfad zeigte dort woandershin.
    if args.vmd:
        args.vmd = os.path.abspath(args.vmd)
    args.root = os.path.abspath(args.root)

    print("=" * 70)
    print("run_allVMD.py - Stapellauf der VMD-ESP-Pipeline")
    print("=" * 70)
    if is_reference:
        print("  Referenzdatensatz (Selbsttest).")
        print(f"  Schreibt nach '{args.images_dir}/', {scene_name} und "
              f"summary_check.csv,")
        print("  damit die committeten Referenzdateien unberuehrt bleiben.")
        print("  Danach die eigenen Zahlen mit reference/summary.csv und die")
        print("  Bilder mit reference/*/images/ vergleichen.")

    entries = discover(args.root)
    if args.only:
        keep = [e for e in entries
                if any(fnmatch.fnmatch(os.path.basename(e["dir"]).lower(),
                                       pat.lower()) for pat in args.only)]
        if len(entries) - len(keep):
            print(f"  --only: {len(entries) - len(keep)} Ordner uebersprungen")
        entries = keep
    if not entries:
        raise SystemExit(f"Keine Molekuelordner unter '{args.root}' gefunden"
                         + (" (mit --only)." if args.only else "."))
    print(f"{len(entries)} Molekuelordner:")
    for e in entries:
        s = os.path.basename(e["struct"]) if e["struct"] else "nur Cubes"
        print(f"  - {e['dir']}  ({s})")

    # --- Schritt 1: Cubes und Kennzahlen ------------------------------------
    print("\n" + "-" * 70)
    print("Schritt 1: Cubes und Kennzahlen")
    print("-" * 70)
    for e in entries:
        print(f"\n[{os.path.basename(os.path.normpath(e['dir']))}]")
        e["cubes"], e["stats"] = ensure_cubes(e, args)
        e["grid"] = cube_dims(e["cubes"]["td"])
        if e["stats"]:
            v0, v1, amp, n = e["stats"]
            print(f"    V_S,min = {v0:+.5f}   V_S,max = {v1:+.5f} a.u.  "
                  f"({n} Punkte)  -> +/- {amp:.4f}")
        else:
            print("    ! keine Schale gefunden - Kennzahlen fehlen")

    # --- Farbskala festlegen ------------------------------------------------
    # needed = die Skala, die ALLE Molekuele dieses Laufs abdeckt. Sie wird
    # immer berechnet, auch im Modus auto: aus ihr entsteht die Empfehlung am
    # Ende. Deshalb braucht dieses Skript keinen zweiten Renderdurchgang -
    # der Vorschlag steht in der Konsole, der Anwender entscheidet, und der
    # naechste Aufruf mit --esp-range <wert> ist der Durchgang, der zaehlt.
    ranges = [e["stats"][2] for e in entries if e["stats"]]
    needed = max(ranges) if ranges else None
    common = needed
    mode = str(args.esp_range).lower()
    if mode == "auto":
        print("\nFarbskala je Molekuel aus den eigenen Kennzahlen (auto)")
    elif mode == "common":
        if common is None:
            raise SystemExit("Keine Kennzahlen - --esp-range braucht einen "
                             "festen Wert.")
        print(f"\nGemeinsame Farbskala: +/- {common:.4f} a.u. "
              f"(groesster Wert von {len(ranges)} Molekuel(en))")
    else:
        try:
            common = float(args.esp_range)
        except ValueError:
            raise SystemExit(f"--esp-range: 'auto', 'common' oder eine Zahl, "
                             f"nicht '{args.esp_range}'.")
        print(f"\nFeste Farbskala: +/- {common:.4f} a.u.")

    # --- Schritt 2: Szene und Bilder ---------------------------------------
    print("\n" + "-" * 70)
    print(f"Schritt 2: {scene_name} schreiben" +
          ("" if args.no_render else " und rendern"))
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
            fill=args.fill, colorscale=args.color_scale,
            rainbow=args.rainbow,
            sources=", ".join(os.path.basename(c)
                              for c in e["cubes"].values()))
        print(f"    -> {tcl}   (Farbskala +/- {rng:.4f} a.u.)")

        row = {"prefix": name, "struct": e["struct"], "grid": e["grid"],
               "stats": e["stats"], "rng": rng, "iso": args.iso, "mode": mode,
               "color_scale": args.color_scale, "opacity": args.opacity,
               "ao": args.ao, "backgrounds": args.backgrounds,
               "made": {}, "size": None, "renderer": ""}

        if not args.no_render:
            # render_espVMD arbeitet im aktuellen Verzeichnis - esp.tcl, die
            # Cubes und images/ liegen dort alle nebeneinander.
            cwd = os.getcwd()
            try:
                os.chdir(e["dir"])
                res = rend.render_all(
                    outdir=args.images_dir, prefix=name, iso=args.iso,
                    rng=rng, stats=e["stats"], vmd=args.vmd, res=args.res,
                    ao=args.ao, backgrounds=args.backgrounds,
                    scene=scene_name, rainbow=args.rainbow,
                    keep_tga=args.keep_tga, dpi=args.dpi, verbose=True)
            finally:
                os.chdir(cwd)
            row.update(made=res["made"], size=res["size"],
                       renderer=res["renderer"])
        rows.append(row)

    summary = args.summary or os.path.join(
        args.root, "summary_check.csv" if is_reference else "summary.csv")
    write_summary(summary, rows,
                  common=None if mode == "auto" else common,
                  advice=needed if mode == "auto" else None)

    # --- Abschluss ----------------------------------------------------------
    print("\n" + "-" * 70)
    print(f"{'Molekuel':<22}{'V_S,min':>10}{'V_S,max':>10}{'Skala':>9}"
          f"  Ansichten")
    print("-" * 70)
    for r in rows:
        s = r["stats"]
        v0 = f"{s[0]:>+10.4f}" if s else f"{'-':>10}"
        v1 = f"{s[1]:>+10.4f}" if s else f"{'-':>10}"
        views = ", ".join(v for v in rend.VIEWS if v in r["made"]) or "-"
        print(f"{r['prefix']:<22}{v0}{v1}{r['rng']:>9.4f}  {views}")
    print("-" * 70)
    if mode != "auto" and common is not None:
        print(f"Gemeinsame Skala: +/- {common:.4f} a.u.")
    print(f"CSV: {summary}")
    if is_reference:
        print("Vergleiche jetzt mit reference/summary.csv und "
              "reference/*/images/.")
    advise(mode, common, needed, rows, raw_argv)
    return 0


if __name__ == "__main__":
    sys.exit(main())

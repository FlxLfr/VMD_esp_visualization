#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_reference.py - kleinen Referenzdatensatz aus einem vollen Molekuelordner

    python tools/make_reference.py sandbox/brombenzol --name brombenzol

Schreibt nach reference/<name>/ ein ausgeduenntes td.xyz / tp.xyz im
Turbomole-pointval-Format plus die Strukturdatei. Das ist der Datensatz fuer
den Selbsttest: `python run_allVMD.py` ohne Argumente laeuft darauf.

Warum pointval und nicht Cube
-----------------------------
Der Selbsttest soll die ganze Kette pruefen, also auch xyzToCubeToVMDVis.py -
Einheitenumrechnung und Indexvertauschung sind die Stellen, an denen am ehesten
etwas kaputtgeht. Faendet er fertige Cubes vor, wuerde er genau diesen Schritt
ueberspringen.

Was ausgeduennt wird
--------------------
Zwei Schritte, beide noetig:

  1. **Zuschneiden.** Das volle Gitter ist eine 30-Bohr-Box, das Molekuel und
     seine rho=0.001-Flaeche fuellen davon nur die Mitte. Der Ausschnitt wird
     aus der DICHTE bestimmt: Huellquader aller Punkte mit rho > iso/2, plus
     Rand. Damit ist garantiert, dass die Isoflaeche vollstaendig enthalten ist
     und die Bilder nicht angeschnitten aussehen.
  2. **Ausduennen.** Nur jeder n-te Punkt je Achse.

Beides zusammen bringt 1.25 GB auf wenige MB. Der Datensatz ist bewusst zu grob
fuer einen zitierfaehigen V_S-Wert - er beantwortet die Frage "laeuft die
Installation und kommen die dokumentierten Zahlen heraus", nicht "wie gross ist
das sigma-Loch".

Dieselben Parameter wie im Schwesterprojekt
-------------------------------------------
Die PyMOL-Pipeline baut ihren Referenzsatz mit denselben Vorgaben (brombenzol,
--stride 5, --margin 2.5). Nur so stehen in beiden Repositorien Zahlen, die man
nebeneinanderlegen darf: weicht die Gitteraufloesung ab, weichen V_S,min und
V_S,max um ein bis drei Prozent ab, und der Vergleich misst dann die
Ausduennung statt der Viewer.
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
    """Gaussian-Cube -> (werte[i1,i2,i3], origin, delta) in Bohr."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        fh.readline()
        fh.readline()
        parts = fh.readline().split()
        natoms, origin = abs(int(parts[0])), [float(v) for v in parts[1:4]]
        dims, delta = [], []
        for i in range(3):
            p = fh.readline().split()
            dims.append(int(p[0]))
            delta.append(float(p[1 + i]))     # Diagonalelement
        for _ in range(natoms):
            fh.readline()
        # np.fromstring parst in C und legt kein Zwischenergebnis in Python an.
        # Ein volles Gitter hat 15.8 Millionen Werte; ueber .split() waere das
        # eine Liste aus 15.8 Millionen str-Objekten, gut ein Gigabyte nur fuer
        # den Zwischenschritt.
        values = np.fromstring(fh.read(), dtype=np.float32, sep=" ")
    if values.size != dims[0] * dims[1] * dims[2]:
        raise ValueError(f"{path}: {values.size} Werte, erwartet "
                         f"{dims[0] * dims[1] * dims[2]}")
    return values.reshape(dims), np.array(origin), np.array(delta)


def load(folder, tag, verbose=True):
    """Cube bevorzugen (Sekunden), sonst die pointval-Datei (Minuten)."""
    cube = os.path.join(folder, f"{tag}.cube")
    if os.path.exists(cube):
        if verbose:
            print(f"    lese {tag}.cube")
        return read_cube(cube)
    raw = os.path.join(folder, f"{tag}.xyz")
    if not os.path.exists(raw):
        raise SystemExit(f"{folder}: weder {tag}.cube noch {tag}.xyz")
    if verbose:
        print(f"    lese {tag}.xyz (pointval, das dauert)")
    import xyzToCubeToVMDVis as conv
    info, data = conv.read_values(raw, verbose=verbose)
    origin = np.array([info["grid"][i][0] for i in range(3)])
    delta = np.array([info["grid"][i][1] for i in range(3)])
    return data, origin, delta


def write_pointval(path, data, origin, delta, quantity, title="101"):
    """Turbomole-pointval schreiben: x laeuft am schnellsten."""
    n1, n2, n3 = data.shape
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        # #origin bleibt bei null, der Versatz steht in den #grid-Zeilen.
        # Turbomole schreibt es genauso, und xyzToCubeToVMDVis.py addiert
        # beide: stuende der Versatz an beiden Stellen, laege das Gitter
        # doppelt so weit vom Ursprung weg wie die Atome - der klassische
        # "Molekuel schwebt neben seiner Oberflaeche"-Fehler.
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
        for k in range(n3):                      # z aussen ...
            for j in range(n2):                  # ... x innen
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
        description="Ausgeduennten Referenzdatensatz fuer den Selbsttest bauen.")
    p.add_argument("source", help="Molekuelordner mit td/tp und Struktur")
    p.add_argument("--name", help="Name in reference/ (Standard: Ordnername)")
    p.add_argument("--outdir", default=None,
                   help="Standard: reference/ des Repositoriums")
    p.add_argument("--stride", type=int, default=5,
                   help="jeder n-te Gitterpunkt je Achse (Standard 5)")
    p.add_argument("--margin", type=float, default=2.5,
                   help="Rand um die Isoflaeche in Bohr (Standard 2.5)")
    p.add_argument("--iso", type=float, default=0.001,
                   help="Isowert, der vollstaendig enthalten sein muss")
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
        raise SystemExit("td und tp liegen nicht auf demselben Gitter.")
    print(f"[1] Gitter {dens.shape[0]}x{dens.shape[1]}x{dens.shape[2]}, "
          f"delta {delta[0]:.4f} Bohr, Ursprung {origin[0]:.2f}")

    # Ausschnitt aus der Dichte: alles, was die Isoflaeche braucht.
    mask = dens > args.iso / 2.0
    if not mask.any():
        raise SystemExit(f"Keine Punkte mit rho > {args.iso/2:g} gefunden.")
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
    print(f"[2] Ausschnitt {sub_d.shape[0]}x{sub_d.shape[1]}x{sub_d.shape[2]}, "
          f"delta {new_delta[0]:.4f} Bohr "
          f"({sub_d.size:,} statt {dens.size:,} Punkte)")
    print(f"    rho {sub_d.min():.2e} .. {sub_d.max():.4g} | "
          f"ESP {sub_e.min():+.4g} .. {sub_e.max():+.4g}")

    for tag, arr in (("td", sub_d), ("tp", sub_e)):
        out = os.path.join(dest, f"{tag}.xyz")
        write_pointval(out, arr, new_origin, new_delta, QUANTITY[tag])
        print(f"[3] -> {out}  ({os.path.getsize(out)/1024**2:.1f} MB)")

    # Strukturdatei mitnehmen - ohne sie findet run_allVMD.py den Ordner nicht.
    for fn in sorted(os.listdir(src)):
        stem, ext = os.path.splitext(fn)
        if ext.lower() in (".mol", ".sdf", ".xyz") and stem not in GRIDS:
            shutil.copy(os.path.join(src, fn), os.path.join(dest, fn))
            print(f"[4] -> {os.path.join(dest, fn)}")
            break
    else:
        print("    ! keine Strukturdatei gefunden - von Hand nachlegen",
              file=sys.stderr)

    print("\nSelbsttest:  python run_allVMD.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())

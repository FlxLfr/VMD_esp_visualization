#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
xyzToCubeToVMDVis.py - Turbomole pointval -> Gaussian Cube -> VMD-Szene

    python xyzToCubeToVMDVis.py --struct brombenzol_aro_opti.mol td.xyz tp.xyz
    vmd -e esp.tcl

Schreibt td.cube, tp.cube und ein esp.tcl, das beide Cubes laedt, die
Isoflaeche der Elektronendichte erzeugt und sie nach dem elektrostatischen
Potential einfaerbt. Die Szene selbst steht in esp_template.tcl daneben.

Umfang: Fuer die Farbskala werden V_S,min und V_S,max auf der rho=iso-Schale
bestimmt - ohne die Skala kein sinnvolles Bild. Genommen werden dafuer die
Gitterpunkte im Band rho = iso +/- 12 %, mit denselben Konstanten wie
xyzToCube.py im Schwesterprojekt Pymol_esp_visualization. Es ist derselbe
Rechenweg, und die Zahlen stimmen entsprechend ueberein - genau das prueft der
Selbsttest.

Bewusst NICHT enthalten ist das sigma-Loch. Es zu finden verlangt trilineare
Interpolation entlang der C-X-Achse, und dieser Code steht drueben,
parameterstudiert und dokumentiert; ein zweites, unabhaengig gepflegtes
Exemplar davon waere genau die Art Duplikat, die unbemerkt auseinanderlaeuft.
Fuer sigma-Loch-Werte also in die PyMOL-Pipeline schauen. V_S,min und V_S,max
sind hier dieselben.

Zwei Stolpersteine im Datenformat, die das Skript abfaengt:
  1. Achsenreihenfolge - Turbomole variiert x am schnellsten, Cube z.
  2. Einheiten - das Gitter steht in Bohr, die Strukturdatei meist in Angstrom.

Die Strukturdatei liefert nur den Atomblock des Cube-Headers (die pointval-
Datei enthaelt keine Atome). Danach braucht VMD sie nicht mehr.

Nur numpy wird benoetigt.
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

CHUNK_BYTES = 1 << 24          # 16 MB Lesepuffer
TEMPLATE = "esp_template.tcl"  # neben diesem Skript

# Schalendicke relativ zum Isowert. Duenn genug, dass die Werte wirklich von
# der Isoflaeche stammen; ist sie dadurch zu schwach besetzt, einmal aufweiten.
SHELL_TOL, SHELL_TOL_WIDE, SHELL_MIN_POINTS = 0.12, 0.30, 50


# ----------------------------------------------------------------------------
# Strukturdatei -> Atome in Bohr
# ----------------------------------------------------------------------------

def _z(sym, path):
    key = sym.strip().capitalize().upper()
    if key in SYMBOL_TO_Z:
        return SYMBOL_TO_Z[key]
    if re.fullmatch(r"\d+", sym.strip()):
        return int(sym)
    raise ValueError(f"Unbekanntes Element '{sym}' in {path}.")


def _read_xyz(lines, path):
    """``Symbol x y z``, mit oder ohne die zwei Kopfzeilen.

    Die Turbomole-Strukturdateien hier haben keine - deshalb lassen sie sich
    auch nicht direkt in VMD laden ("Unable to load molecule"). Ueber den
    Cube-Umweg ist das egal, die Atome landen im Cube-Header.
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
    """MDL-Molfile V2000/V3000. Koordinaten stehen VOR dem Symbol, in Angstrom."""
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
        raise ValueError(f"{path}: zu kurz fuer ein Molfile.")
    try:
        natoms = int(lines[3][0:3])          # Zaehlzeile
    except ValueError:
        raise ValueError(f"{path}: Zaehlzeile nicht lesbar: {lines[3].strip()!r}")

    atoms = []
    for line in lines[4:4 + natoms]:
        p = line.split()
        if len(p) < 4:
            raise ValueError(f"{path}: Atomzeile unvollstaendig: {line.strip()!r}")
        atoms.append((_z(p[3], path), float(p[0]), float(p[1]), float(p[2])))
    return atoms


def read_structure(path, unit="angstrom"):
    """.xyz / .mol / .sdf -> Liste von (Z, x, y, z) in Bohr."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    if not lines:
        raise ValueError(f"{path} ist leer.")

    if os.path.splitext(path)[1].lower() in (".mol", ".sdf", ".sd"):
        for i, ln in enumerate(lines):
            if ln.startswith("$$$$"):        # SD-File: nur der erste Datensatz
                lines = lines[:i]
                break
        atoms, unit = _read_molfile(lines, path), "angstrom"
    else:
        atoms = _read_xyz(lines, path)

    if not atoms:
        raise ValueError(f"Keine Atome in {path} gefunden.")
    if unit == "angstrom":
        b = BOHR_PER_ANGSTROM
        atoms = [(z, x * b, y * b, zz * b) for (z, x, y, zz) in atoms]
    elif unit != "bohr":
        raise ValueError("unit muss 'angstrom' oder 'bohr' sein")
    return atoms


# ----------------------------------------------------------------------------
# Turbomole-Gitterdatei
# ----------------------------------------------------------------------------

def parse_header(fh):
    """Liest den ``#``-Kopf. Rueckgabe: (info, bereits gelesene erste Datenzeile)."""
    info = {"origin": np.zeros(3), "vectors": np.eye(3),
            "grid": [None] * 3, "title": "", "quantity": ""}
    while True:
        line = fh.readline()
        if not line:
            raise ValueError("Datei endet im Header - keine Daten gefunden.")
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
                raise ValueError(f"Gitterzeile nicht lesbar: {line!r}")
            info["grid"][int(body[4]) - 1] = (float(m.group(1)),
                                              float(m.group(2)), int(m.group(3)))
        elif low.startswith("title"):
            info["title"] = body
        elif ("potential" in low or "density" in low) \
                and not low.startswith("cartesian"):
            info["quantity"] = body

    if any(g is None for g in info["grid"]):
        raise ValueError("Header unvollstaendig: #grid1/#grid2/#grid3 fehlen.")
    return info, line


def read_values(path, verbose=True):
    """4. Spalte einer pointval-Datei als float32-Gitter [i1, i2, i3].

    Blockweise statt zeilenweise - fuer die 1.25-GB-Dateien rund eine
    Groessenordnung schneller.
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
                raise ValueError(f"{path}: erwartete 4 Spalten pro Zeile.")
            arr = np.asarray(tokens, dtype=np.float32).reshape(-1, 4)[:, 3]
            if filled + arr.size > total:
                raise ValueError(f"{path}: mehr Datenpunkte als im Header.")
            values[filled:filled + arr.size] = arr
            filled += arr.size
            if verbose:
                sys.stdout.write(f"\r    lese {os.path.basename(path)}: "
                                 f"{100.0 * filled / total:5.1f} %")
                sys.stdout.flush()

        tokens = remainder.split()
        if tokens:
            if len(tokens) % 4:
                raise ValueError(f"{path}: letzte Zeile unvollstaendig.")
            arr = np.asarray(tokens, dtype=np.float32).reshape(-1, 4)[:, 3]
            values[filled:filled + arr.size] = arr
            filled += arr.size

    if verbose:
        sys.stdout.write(f"\r    lese {os.path.basename(path)}: 100.0 %  "
                         f"({filled:,} Punkte in {time.time() - t0:.1f} s)\n")
    if filled != total:
        raise ValueError(f"{path}: {filled} Werte gelesen, erwartet {total}.")

    # Turbomole: x am schnellsten -> [i3, i2, i1]. Cube: z am schnellsten.
    return info, np.ascontiguousarray(
        np.transpose(values.reshape(n3, n2, n1), (2, 1, 0)))


# ----------------------------------------------------------------------------
# Farbskala aus der Isoflaeche
# ----------------------------------------------------------------------------

def shell_range(density, esp, iso=0.001, step=0.005):
    """V_S,min / V_S,max auf der rho=iso-Schale und die symmetrische Skala.

    Den Bereich NICHT aus dem ganzen Gitter nehmen - dort dominieren die
    Kernsingularitaeten mit mehreren hundert a.u. Nur Gitterpunkte nahe der
    Isoflaeche zaehlen, das Ergebnis wird symmetrisch auf ein Vielfaches von
    ``step`` aufgerundet.

    Rueckgabe: (vmin, vmax, halbe_breite, anzahl_punkte) oder None.
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
    """Kennzahlen fuer die erste Cube-Zeile, damit --tcl-only sie wiederfindet."""
    vmin, vmax, amp, npts = stats
    return (f" | V_S,min={vmin:+.6f} V_S,max={vmax:+.6f} range={amp:.4f}"
            f" iso={iso:g} n={npts}")


def read_cube(path):
    """Gaussian-Cube -> float32-Gitter [i1, i2, i3].

    Nur so viel Parser, wie die Schalenstatistik braucht. Gebraucht wird das
    fuer Cubes, die noch keine Kennzahlen im Kopf tragen (aeltere Laeufe) - die
    ~15 s Lesen sind immer noch besser als Minuten Neukonvertierung.
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
        raise ValueError(f"{path}: {values.size} Werte, erwartet "
                         f"{dims[0] * dims[1] * dims[2]}.")
    return values.reshape(dims)


def read_cube_atoms(path):
    """Atomblock eines Cubes -> [(Z, x, y, z), ...], Laengen in Bohr.

    Der Cube traegt die Atome selbst - deshalb laesst sich die Orientierung
    auch mit --tcl-only bestimmen, ohne die Strukturdatei wiederzufinden.
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
            # Spalte 1 Ordnungszahl, Spalte 2 Kernladung, dann x y z
            atoms.append((int(p[0]), float(p[2]), float(p[3]), float(p[4])))
    return atoms


def read_stats(cube_path):
    """Kennzahlen aus der ersten Zeile eines frueher geschriebenen Cubes."""
    with open(cube_path, "r", encoding="utf-8", errors="replace") as fh:
        m = STATS_RE.search(fh.readline())
    if not m:
        return None
    return (float(m.group(1)), float(m.group(2)), float(m.group(3)),
            int(m.group(4)) if m.group(4) else 0)


# ----------------------------------------------------------------------------
# Orientierung der drei Ansichten
#
# Diese Rechnung ist aus der PyMOL-Pipeline uebernommen (render_esp.py,
# molecular_frame/view_matrix) und liefert absichtlich dieselben Achsen: nur
# dann zeigen beide Bildersaetze dasselbe und der Vergleich misst den Viewer.
#
# Vorher stand die Orientierung als feste Drehung im Tcl-Template - "pi, dann
# rotate x by -90". Das setzt voraus, dass das Molekuel planar in der xy-Ebene
# liegt und die C-X-Achse nach -y zeigt. Fuer die Halogenbenzole, wie Turbomole
# sie ablegt, stimmte das zufaellig. Fuer die substituierten Pyridine nicht:
# dort blickte die sigma-Ansicht am Loch vorbei. Auffallen tut das kaum, weil
# in jeder Richtung eine runde, eingefaerbte Flaeche zu sehen ist - man haelt
# die Gegenseite fuer das sigma-Loch.
# ----------------------------------------------------------------------------

HALOGENS = {9: "F", 17: "Cl", 35: "Br", 53: "I"}


def molecular_frame(atoms, halogen_index=None):
    """Reproduzierbares Molekuelkoordinatensystem aus der Geometrie.

    Rueckgabe: (normal, axis, sigma_axis, label)
      normal      Flaechennormale: Hauptachse der kleinsten Ausdehnung ueber
                  die Schweratome (Wasserstoff wackelt und traegt nichts bei)
      axis        die IN DIE EBENE PROJIZIERTE C->Halogen-Achse; sie spannt mit
                  normal ein Rechtssystem auf und richtet pi und edge aus
      sigma_axis  die ECHTE, unprojizierte C->Halogen-Achse
      label       woher die Achse kommt, fuer die Meldung im Skript

    Warum zwei Achsen: bei planaren Molekuelen sind sie identisch und die
    Projektion ist Rundungskosmetik. Steht die C-X-Bindung aus der
    Ausgleichsebene heraus, dreht sie die Achse aber wirklich weg - und die
    sigma-Ansicht blickte um genau diesen Winkel am Loch vorbei.

    Ohne Halogen gibt es kein sigma-Loch. Dann nimmt die Achse die groesste
    Traegheitsausdehnung, also die Laengsachse; die sigma-Ansicht ist dort der
    Blick von der Schmalseite und keine Aussage ueber ein Loch.
    """
    coords = np.array([[a[1], a[2], a[3]] for a in atoms], dtype=float)
    znums = np.array([a[0] for a in atoms])

    heavy = coords[znums > 1]
    if len(heavy) < 3:
        heavy = coords
    _, _, vt = np.linalg.svd(heavy - heavy.mean(axis=0), full_matrices=False)
    normal = vt[2]                                  # kleinste Ausdehnung
    long_axis = vt[0]                               # groesste Ausdehnung

    axis, label = None, "laengste Traegheitsachse (kein Halogen)"
    hal_idx = [i for i, z in enumerate(znums) if z in HALOGENS]
    if hal_idx:
        hi = halogen_index if halogen_index in hal_idx else hal_idx[0]
        carbons = [i for i, z in enumerate(znums) if z == 6]
        if carbons:
            d = np.linalg.norm(coords[carbons] - coords[hi], axis=1)
            ci = carbons[int(np.argmin(d))]
            axis = coords[hi] - coords[ci]          # C -> X, zeigt zum Loch
            label = (f"C{ci + 1}->{HALOGENS[int(znums[hi])]}{hi + 1}")

    if axis is None:
        axis = long_axis.copy()

    axis = axis / np.linalg.norm(axis)
    normal = normal / np.linalg.norm(normal)
    sigma_axis = axis.copy()

    axis = axis - normal * float(np.dot(axis, normal))
    if np.linalg.norm(axis) < 1e-6:                 # C-X steht senkrecht
        axis = long_axis
    axis = axis / np.linalg.norm(axis)

    return normal, axis, sigma_axis, label


def view_matrix(forward, up):
    """Rotationsmatrix fuer VMDs ``rotate_matrix``.

    ``forward`` zeigt vom Molekuel zur Kamera, ``up`` nach oben im Bild.

    VMD rechnet Bildkoordinaten als R * (Welt - Zentrum). Damit die Kamera
    entlang forward blickt, muss R genau diesen Vektor auf +z abbilden - die
    ZEILEN von R sind also die Kamerabasis. (PyMOLs set_view will dieselbe
    Matrix transponiert; deshalb steht dort ein .T und hier keines.)
    """
    z = np.asarray(forward, dtype=float)
    z = z / np.linalg.norm(z)
    up = np.asarray(up, dtype=float)
    up = up - z * float(np.dot(up, z))
    if np.linalg.norm(up) < 1e-8:                   # up parallel zu forward
        alt = np.array([0.0, 0.0, 1.0])
        if abs(float(np.dot(alt, z))) > 0.9:
            alt = np.array([1.0, 0.0, 0.0])
        up = alt - z * float(np.dot(alt, z))
    up = up / np.linalg.norm(up)
    right = np.cross(up, z)
    right = right / np.linalg.norm(right)
    return np.array([right, up, z])


def tcl_matrix(m):
    """3x3-Rotation -> 4x4-Matrix als Tcl-Liste, wie molinfo sie erwartet."""
    rows = []
    for i in range(3):
        rows.append("{%s 0.000000}" % " ".join(f"{v:.6f}" for v in m[i]))
    rows.append("{0.000000 0.000000 0.000000 1.000000}")
    return "{%s}" % " ".join(rows)


# Fallback ohne Atome: genau die Achsen, die die frueheren festen Drehungen
# ergaben (pi = Blick aus +z, edge = rotate y by 90, sigma = rotate x by -90).
# So aendert sich nichts, wenn die Geometrie ausnahmsweise fehlt.
_FALLBACK = {
    "pi":    ((0.0, 0.0, 1.0), (0.0, 1.0, 0.0)),
    "edge":  ((-1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
    "sigma": ((0.0, -1.0, 0.0), (0.0, 0.0, 1.0)),
}


def view_matrices(atoms):
    """Die drei Ansichten als Tcl-Matrizen plus eine Notiz zur sigma-Achse."""
    if not atoms:
        mats = {k: view_matrix(f, u) for k, (f, u) in _FALLBACK.items()}
        return {k: tcl_matrix(m) for k, m in mats.items()}, \
               "Weltachsen (keine Geometrie gefunden)", None

    normal, axis, sigma_axis, label = molecular_frame(atoms)
    mats = {
        # Blick senkrecht auf die Ebene, C-X-Achse zeigt nach unten
        "pi":    view_matrix(forward=normal, up=-axis),
        # Blick in der Ebene, senkrecht zur C-X-Achse; C-X-Achse waagerecht
        "edge":  view_matrix(forward=np.cross(normal, axis), up=normal),
        # Blick von aussen entlang der ECHTEN C-X-Achse auf das sigma-Loch
        "sigma": view_matrix(forward=sigma_axis, up=normal),
    }
    # Neigung der C-X-Achse gegen die Ausgleichsebene: ist sie gross, laufen
    # sigma und pi/edge auseinander, und das sollte im Protokoll stehen.
    tilt = abs(90.0 - math.degrees(math.acos(
        min(1.0, abs(float(np.dot(sigma_axis, normal)))))))
    return {k: tcl_matrix(m) for k, m in mats.items()}, label, tilt


# ----------------------------------------------------------------------------
# Cube schreiben
# ----------------------------------------------------------------------------

def write_cube(path, info, data, atoms, stride=1, comment=""):
    """Gaussian-Cube, alle Laengen in Bohr.

    Geschrieben wird erst nach <name>.part und am Ende umbenannt. Ein voller
    Cube ist 200 MB und braucht Minuten; wird der Lauf in dieser Zeit
    abgebrochen - Strg-C, geschlossenes Fenster, volle Platte -, bliebe sonst
    eine halbe Datei mit gueltigem Kopf liegen. Der naechste Lauf haelt die fuer
    fertig, ueberspringt die Konvertierung und rendert eine Isoflaeche mit
    fehlender hinterer Haelfte. Das Umbenennen ist der Moment, in dem die Datei
    entsteht - vorher gibt es sie unter ihrem Namen nicht.
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
        fh.write(f"{comment or 'erzeugt mit xyzToCubeToVMDVis.py'}\n")
        fh.write(f"{info.get('quantity') or 'volumetric data'} | "
                 f"{info.get('title', '')} | Einheiten: Bohr\n")
        fh.write(f"{len(atoms):5d} {origin[0]:12.6f} {origin[1]:12.6f} "
                 f"{origin[2]:12.6f}\n")
        for i in range(3):
            fh.write(f"{n[i]:5d} {voxel[i][0]:12.6f} {voxel[i][1]:12.6f} "
                     f"{voxel[i][2]:12.6f}\n")
        for (znum, x, y, z) in atoms:
            fh.write(f"{znum:5d} {float(znum):12.6f} {x:12.6f} {y:12.6f} "
                     f"{z:12.6f}\n")

        # Werte: z am schnellsten, 6 pro Zeile. Ein vorkompiliertes
        # Formatmuster pro z-Reihe ist deutlich schneller als eine Schleife
        # ueber alle Einzelwerte.
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
    # Erst hier existiert <name>.cube. os.replace ersetzt auch eine vorhandene
    # Datei und ist auf einem Dateisystem unteilbar - es gibt keinen Moment, in
    # dem der alte Cube weg und der neue noch nicht da ist.
    os.replace(tmp, path)
    return n


# ----------------------------------------------------------------------------
# VMD-Szene
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Farbrampen
#
# Dieselben Stuetzfarben wie in der PyMOL-Pipeline (render_esp.py, RAMP_HEX),
# damit beide Bildersaetze wirklich dieselbe Rampe zeigen und nicht nur
# dieselbe Richtung. VMDs eingebaute Skalen kennen nur drei Stuetzfarben; die
# fuenf des Regenbogens werden deshalb ueber esp_ramp in esp_template.tcl in
# die 1024 Eintraege der Farbtabelle geschrieben.
#
# Rot-weiss-blau bleibt vorerst VMDs eigenes RWB (reines Rot und Blau statt
# PyMOLs etwas dunklerem #d40000/#0030d4). Ein Wechsel wuerde jeden bereits
# gerenderten Standardsatz minimal veraendern.
# ----------------------------------------------------------------------------

RAMP_HEX = {
    "rainbow": ["#d40000", "#f0e000", "#00a000", "#00c8d4", "#0030d4"],
    "redblue": ["#d40000", "#ffffff", "#0030d4"],
}


def ramp_stops(name):
    """Stuetzfarben als Tcl-Liste {{r g b} {r g b} ...}, Werte 0..1."""
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
                     opacity=0.50, scale="auto", fill=0.85, colorscale="RWB",
                     sources="", atoms=None, rainbow=False):
    """Fuellt esp_template.tcl.

    Ueber @@PLATZHALTER@@ statt str.format oder string.Template: Tcl benutzt
    geschweifte Klammern und Dollarzeichen als Syntax, beide Mechanismen wuerden
    genau daran ersticken.

    rho_cube wird zuerst geladen und ist damit Volume 0, esp_cube Volume 1.
    """
    tpl_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), TEMPLATE)
    with open(tpl_path, "r", encoding="utf-8") as fh:
        text = fh.read()

    note = ""
    if stats:
        note = (f"\\nV_S,min = {stats[0]:+.5f}   V_S,max = {stats[1]:+.5f} a.u. "
                f"(Gitterpunkte auf der Schale)")

    # Orientierung aus der Geometrie. Ohne uebergebene Atome aus dem Cube
    # lesen - der traegt sie, auch bei --tcl-only.
    if atoms is None:
        try:
            atoms = read_cube_atoms(
                os.path.join(os.path.dirname(os.path.abspath(path)), rho_cube))
        except Exception as err:
            print(f"  ! Atome nicht lesbar ({err}) - Ansichten auf Weltachsen",
                  file=sys.stderr)
            atoms = None
    rot, axis_label, tilt = view_matrices(atoms)
    axis_note = f"\\nsigma-Achse: {axis_label}"
    if tilt is not None and tilt > 15.0:
        # Steht die Bindung schraeg, folgen pi/edge der Ebene und sigma der
        # Bindung - dann sind es nicht mehr drei Ansichten desselben Rahmens.
        axis_note += (f" ({tilt:.0f} Grad aus der Ausgleichsebene; "
                      f"pi/edge folgen der Ebene)")

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
        "@@COLORSCALE@@": colorscale,
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
# Hauptprogramm
# ----------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        description="Turbomole-pointval-Gitter nach Gaussian-Cube konvertieren "
                    "und eine VMD-Szene schreiben.",
        epilog="Beispiel: python xyzToCubeToVMDVis.py "
               "--struct brombenzol_aro_opti.mol td.xyz tp.xyz")
    p.add_argument("grids", nargs="+", help="td.xyz tp.xyz (oder td.cube "
                                            "tp.cube bei --tcl-only)")
    p.add_argument("--struct", "-s", help="Strukturdatei fuer den Cube-Header: "
                                          ".xyz, .mol, .sdf")
    p.add_argument("--struct-unit", choices=["angstrom", "bohr"],
                   default="angstrom", help="nur fuer .xyz relevant")
    p.add_argument("--outdir", "-o", help="Standard: neben der Eingabe")
    p.add_argument("--stride", type=int, default=1,
                   help="jeden n-ten Gitterpunkt (2 => 8x kleiner)")
    p.add_argument("--quiet", "-q", action="store_true")

    g = p.add_argument_group("VMD-Szene (esp.tcl)")
    g.add_argument("--tcl-only", action="store_true",
                   help="nur esp.tcl neu schreiben, Cubes unangetastet lassen "
                        "(Sekundenbruchteile statt Minuten)")
    g.add_argument("--no-vmd", action="store_true", help="nur die Cube-Dateien")
    g.add_argument("--esp-range", default="auto",
                   help="halbe Breite der Farbskala in a.u., oder 'auto' "
                        "(Standard): aus V_S,min/V_S,max auf der Schale")
    g.add_argument("--iso", type=float, default=0.001, help="Isowert (a.u.)")
    g.add_argument("--opacity", type=float, default=0.50,
                   help="Deckkraft der Isoflaeche 0..1")
    g.add_argument("--scale", default="auto",
                   help="Zoom: Zahl, oder 'auto' (Standard) aus der "
                        "Molekuelgroesse und der Fensterhoehe")
    g.add_argument("--fill", type=float, default=0.85,
                   help="Anteil der Fensterhoehe fuer das Molekuel bei "
                        "--scale auto (Standard 0.85)")
    g.add_argument("--color-scale", default=None,
                   help="VMD-Farbskala, Standard RWB (rot negativ); "
                        "mit --rainbow RGB")
    g.add_argument("--rainbow", action="store_true",
                   help="Regenbogenrampe statt rot-weiss-blau. Schreibt "
                        "esp_rainbow.tcl, damit die Standardszene erhalten "
                        "bleibt.")
    args = p.parse_args(argv)

    # RGB ist VMDs Regenbogen: kleinster Wert rot, Mitte gruen, groesster blau -
    # dieselbe Richtung wie die Rampe der PyMOL-Pipeline. Ein ausdrueckliches
    # --color-scale gewinnt, sonst waere die Option neben --rainbow wirkungslos.
    if args.color_scale is None:
        args.color_scale = "RGB" if args.rainbow else "RWB"

    verbose = not args.quiet
    if args.tcl_only and args.no_vmd:
        p.error("--tcl-only und --no-vmd zusammen ergeben nichts zu tun.")
    if not args.tcl_only and not args.struct:
        p.error("--struct wird gebraucht: der Cube-Header verlangt einen "
                "Atomblock, die pointval-Datei enthaelt keine Atome.")

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
                p.error(f"{cube} fehlt - einmal ohne --tcl-only laufen lassen.")
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
            print(f"[1] Struktur: {args.struct} -> {len(atoms)} Atome "
                  f"({args.struct_unit} -> Bohr)")

        grids = {}
        for gpath in args.grids:
            if verbose:
                print(f"[2] Gitterdatei: {gpath}")
            info, data = read_values(gpath, verbose=verbose)
            if verbose:
                n = [info["grid"][i][2] for i in range(3)]
                print(f"    {n[0]} x {n[1]} x {n[2]}, "
                      f"delta = {info['grid'][0][1]} Bohr, "
                      f"'{info['quantity'] or 'unbekannt'}', "
                      f"{data.min():+.4g} .. {data.max():+.4g}")
            grids[classify(info["quantity"]) or gpath] = (info, data, gpath)

        if "density" in grids and "esp" in grids:
            stats = shell_range(grids["density"][1], grids["esp"][1], args.iso)
            if verbose and stats:
                v0, v1, amp, npts = stats
                print(f"[3] Schale rho = {args.iso:g} a.u.: {npts:,} Punkte")
                print(f"    V_S,min = {v0:+.5f} a.u. "
                      f"({v0 * HARTREE_TO_KJ:+7.1f} kJ/(mol*e))")
                print(f"    V_S,max = {v1:+.5f} a.u. "
                      f"({v1 * HARTREE_TO_KJ:+7.1f} kJ/(mol*e))")
                print(f"    -> Farbskala +/- {amp:.4f} a.u.")

        for kind, (info, data, gpath) in grids.items():
            base = os.path.splitext(os.path.basename(gpath))[0]
            od = args.outdir or os.path.dirname(os.path.abspath(gpath))
            os.makedirs(od, exist_ok=True)
            out = os.path.join(od, base + ".cube")
            comment = f"{info['quantity'] or base} - aus {os.path.basename(gpath)}"
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
        print("    ! Dichte oder Potential fehlt - kein esp.tcl. Die Dichte "
              "liefert die Flaeche, das Potential die Farbe.", file=sys.stderr)
        return 1

    if str(args.esp_range).lower() == "auto":
        if stats:
            rng = stats[2]
        else:
            rng = 0.035
            print("    ! Farbskala nicht bestimmbar - benutze +/- 0.035 a.u.",
                  file=sys.stderr)
    else:
        rng = float(args.esp_range)

    # Eigener Dateiname, sonst ueberschreibt ein Regenbogenlauf die Szene des
    # rot-weiss-blauen Laufs - die Bilder liegen aus demselben Grund getrennt.
    tcl = os.path.join(outdir, "esp_rainbow.tcl" if args.rainbow else "esp.tcl")
    write_vmd_script(tcl, os.path.basename(cubes["density"]),
                     os.path.basename(cubes["esp"]), rng, stats, iso=args.iso,
                     opacity=args.opacity,
                     scale=(args.scale if str(args.scale) == "auto"
                            else float(args.scale)),
                     fill=args.fill,
                     colorscale=args.color_scale, rainbow=args.rainbow,
                     sources=", ".join(os.path.basename(g) for g in args.grids))
    if verbose:
        print(f"[4] VMD-Szene: {tcl}   (Farbskala +/- {rng:.4f} a.u.)")
        print(f"    Start mit:  vmd -e {os.path.basename(tcl)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

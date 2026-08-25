#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
xyzToCubeToVMDVis.py
====================

Turbomole-``pointval``-Gitterdateien (td.xyz, tp.xyz) -> Gaussian-Cube ->
fertige VMD-Szene (``esp.tcl``).

Was dieses Skript tut
---------------------
1. Es liest die beiden ASCII-Punktwolken, die Turbomole ausgibt, und schreibt
   sie als Gaussian-Cube-Dateien (``td.cube``, ``tp.cube``).
2. Es schreibt daneben ein ``esp.tcl``, das beide Cubes in VMD laedt, die
   Isoflaeche der Elektronendichte erzeugt und sie nach dem elektrostatischen
   Potential einfaerbt.

Was dieses Skript **nicht** tut
-------------------------------
Es rechnet nichts aus. Keine V_S,min / V_S,max, kein sigma-Loch, keine
Oberflaechenstatistik. Diese Groessen werden im Schwesterprojekt
``Pymol_esp_visualization`` bestimmt; hier waeren sie eine zweite, unabhaengig
gepflegte Implementierung derselben Zahlen - also genau die Sorte Duplikat, die
irgendwann leise auseinanderlaeuft. Dieses Projekt visualisiert, mehr nicht.

Eine Folge davon betrifft die Farbskala: Sie wird **nicht** aus den Daten
bestimmt, sondern als Parameter uebergeben (``--esp-range``, Standard 0.035 a.u.).
Das ist kein Mangel, sondern der Punkt: Nur wenn beide Pipelines nachweislich
dieselbe Skala benutzen, sind ihre Bilder ueberhaupt vergleichbar. Den Wert
liefert die PyMOL-Pipeline in ``<molekuel>_settings.txt``, Zeile "Farbskala".

Hintergrund zum Datenformat
---------------------------
Turbomole schreibt Volumendaten als reine ASCII-Punktwolke::

    #origin           0.000000      0.000000      0.000000
    #vector1          1.000000      0.000000      0.000000
    ...
    #grid1  start  -15.000000  delta    0.120000  points    251
    #plotdata
    # cartesian coordinates x,y,z and f(x,y,z)
          -15.00000000   -15.00000000   -15.00000000   -0.00054019
          ...

Jede Zeile enthaelt die vollen Koordinaten -> bei 251^3 Punkten sind das 1.25 GB.
Eine Cube-Datei speichert dieselbe Information mit implizitem Gitter (~200 MB).

Zwei Stolpersteine, die dieses Skript abfaengt:

1. **Achsenreihenfolge.** In der Turbomole-Datei laeuft *x* am schnellsten,
   im Cube-Format laeuft *z* am schnellsten. Ohne Umsortierung erhaelt man ein
   transponiertes, gespiegeltes Molekuel.
2. **Einheiten.** Das Gitter steht in Bohr, die Strukturdatei ueblicherweise in
   Angstrom. Das Skript rechnet die Atome standardmaessig um
   (``--struct-unit angstrom``).

Warum ueberhaupt eine Strukturdatei noetig ist: Der Cube-Header verlangt einen
Atomblock, und die pointval-Datei enthaelt keine Atome. Sobald der Cube
geschrieben ist, braucht VMD die Strukturdatei nicht mehr - Atome und Gitter
kommen dann aus derselben Datei und koennen nicht mehr gegeneinander verrutschen.

Benutzung
---------
::

    python xyzToCubeToVMDVis.py --struct brombenzol_aro_opti.mol td.xyz tp.xyz

Ergebnis: ``td.cube``, ``tp.cube``, ``esp.tcl``. Danach::

    vmd -e esp.tcl

Wenn das volle 251^3-Gitter zu traege ist, jeden zweiten Punkt nehmen::

    python xyzToCubeToVMDVis.py --struct brombenzol_aro_opti.mol td.xyz tp.xyz --stride 2

Als Strukturdatei werden ``.xyz`` (mit oder ohne Kopfzeilen), ``.mol`` und
``.sdf`` akzeptiert.

Nur numpy wird benoetigt.
"""

from __future__ import annotations

import argparse
import datetime
import os
import re
import sys
import time

import numpy as np

from constants import BOHR_PER_ANGSTROM, HARTREE_TO_KJ

# ----------------------------------------------------------------------------
# Konstanten
# ----------------------------------------------------------------------------

# Elementsymbole in der Reihenfolge der Ordnungszahl: die Position in der Liste
# IST Z-1, daraus wird unten SYMBOL_TO_Z gebaut.
ELEMENTS = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
]
SYMBOL_TO_Z = {sym.upper(): i + 1 for i, sym in enumerate(ELEMENTS)}

CHUNK_BYTES = 1 << 24                      # 16 MB Lesepuffer


# ----------------------------------------------------------------------------
# Strukturdatei einlesen
# ----------------------------------------------------------------------------

def _symbol_to_z(sym: str, path: str) -> int:
    """Elementsymbol (oder Ordnungszahl) -> Ordnungszahl."""
    key = sym.strip().capitalize().upper()
    if key in SYMBOL_TO_Z:
        return SYMBOL_TO_Z[key]
    if re.fullmatch(r"\d+", sym.strip()):
        return int(sym)                                # Ordnungszahl statt Symbol
    raise ValueError(
        f"Unbekanntes Element '{sym}' in {path}. "
        f"Bitte die Liste ELEMENTS im Skript ergaenzen."
    )


def _read_xyz(lines, path):
    """xyz-Format: Zeilen ``Symbol x y z``.

    Akzeptiert sowohl das Standardformat (Atomanzahl + Kommentarzeile + Atome)
    als auch eine nackte Koordinatenliste, wie die Turbomole-Rechnungen sie
    hier liefern. Nebenbei der Grund, warum diese Dateien sich nicht direkt in
    VMD laden lassen: VMDs xyz-Reader besteht auf den zwei Kopfzeilen und
    bricht sonst mit "Unable to load molecule" ab. Ueber den Cube-Umweg ist das
    egal - die Atome landen im Cube-Header.
    """
    start = 0
    first = lines[0].strip() if lines else ""
    if re.fullmatch(r"\d+", first):
        start = 2                                     # Anzahl + Kommentar ueberspringen

    atoms = []
    for line in lines[start:]:
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            x, y, z = (float(parts[1]), float(parts[2]), float(parts[3]))
        except ValueError:
            continue                                   # Kommentar-/Muellzeile
        atoms.append((_symbol_to_z(parts[0], path), x, y, z))
    return atoms


def _read_molfile(lines, path):
    """MDL-Molfile / SD-File (.mol, .sdf), V2000 und V3000.

    Aufbau V2000::

        Zeile 1   Titel
        Zeile 2   Programmzeile
        Zeile 3   Kommentar
        Zeile 4   Zaehlzeile:  " 12 12  0 ... V2000"
        dann      je Atom:  x  y  z  Symbol  ...      <- Koordinaten ZUERST
        dann      Bindungsblock

    Bei SD-Files wird nur der erste Datensatz gelesen (bis ``$$$$``).
    Molfile-Koordinaten sind per Definition in Angstrom.
    """
    atoms = []

    # --- V3000 ---------------------------------------------------------
    if any("V3000" in ln for ln in lines[:8]):
        inside = False
        for line in lines:
            s = line.strip()
            if s.startswith("M  V30 BEGIN ATOM"):
                inside = True
                continue
            if s.startswith("M  V30 END ATOM"):
                break
            if inside and s.startswith("M  V30"):
                # M  V30 <index> <symbol> <x> <y> <z> <aamap> ...
                parts = s.split()
                if len(parts) >= 7:
                    atoms.append((_symbol_to_z(parts[3], path),
                                  float(parts[4]), float(parts[5]),
                                  float(parts[6])))
        return atoms

    # --- V2000 ---------------------------------------------------------
    if len(lines) < 5:
        raise ValueError(f"{path}: zu kurz fuer ein Molfile.")

    counts = lines[3]
    try:
        natoms = int(counts[0:3])
    except ValueError:
        raise ValueError(
            f"{path}: Zaehlzeile (Zeile 4) nicht lesbar: {counts.strip()!r}")

    for line in lines[4:4 + natoms]:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"{path}: Atomzeile unvollstaendig: {line.strip()!r}")
        x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
        atoms.append((_symbol_to_z(parts[3], path), x, y, z))

    return atoms


def read_structure(path: str, unit: str = "angstrom"):
    """Liest eine Strukturdatei und liefert die Atome in **Bohr**.

    Unterstuetzt:

    ==========  =====================================================
    ``.xyz``    ``Symbol x y z``, mit oder ohne Kopfzeilen
    ``.mol``    MDL-Molfile V2000/V3000 (Koordinaten *vor* dem Symbol)
    ``.sdf``    SD-File, erster Datensatz
    ==========  =====================================================

    Rueckgabe: Liste von ``(Z, x, y, z)``.

    ``unit`` gilt nur fuer xyz-Dateien - Molfile-Koordinaten sind per
    Definition in Angstrom, dort wird die Angabe ignoriert.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()

    if not lines:
        raise ValueError(f"{path} ist leer.")

    ext = os.path.splitext(path)[1].lower()

    if ext in (".mol", ".sdf", ".sd"):
        # SD-File: nur der erste Datensatz
        for i, ln in enumerate(lines):
            if ln.startswith("$$$$"):
                lines = lines[:i]
                break
        atoms = _read_molfile(lines, path)
        file_unit = "angstrom"
    else:
        atoms = _read_xyz(lines, path)
        file_unit = unit

    if not atoms:
        raise ValueError(f"Keine Atome in {path} gefunden.")

    if file_unit == "angstrom":
        atoms = [(z, x * BOHR_PER_ANGSTROM,
                     y * BOHR_PER_ANGSTROM,
                     zz * BOHR_PER_ANGSTROM) for (z, x, y, zz) in atoms]
    elif file_unit != "bohr":
        raise ValueError("unit muss 'angstrom' oder 'bohr' sein")

    return atoms


# ----------------------------------------------------------------------------
# Turbomole-Gitterdatei einlesen
# ----------------------------------------------------------------------------

def parse_header(fh):
    """Liest den ``#``-Kopf einer pointval-Datei.

    Rueckgabe: (info-dict, erste_datenzeile). Die erste Datenzeile wurde bereits
    aus dem Stream gelesen und muss vom Aufrufer mitverarbeitet werden.
    """
    info = {
        "origin": np.zeros(3),
        "vectors": np.eye(3),
        "grid": [None, None, None],       # je (start, delta, points)
        "title": "",
        "quantity": "",
    }
    first_data_line = None

    while True:
        line = fh.readline()
        if not line:
            raise ValueError("Datei endet im Header - keine Daten gefunden.")
        if not line.startswith("#"):
            first_data_line = line
            break

        body = line[1:].strip()
        low = body.lower()

        if low.startswith("origin"):
            info["origin"] = np.array([float(v) for v in body.split()[1:4]])
        elif low.startswith("vector"):
            idx = int(body[6]) - 1
            info["vectors"][idx] = np.array([float(v) for v in body.split()[1:4]])
        elif low.startswith("grid"):
            idx = int(body[4]) - 1
            m = re.search(
                r"start\s+(\S+)\s+delta\s+(\S+)\s+points\s+(\d+)", body, re.I)
            if not m:
                raise ValueError(f"Gitterzeile nicht lesbar: {line!r}")
            info["grid"][idx] = (float(m.group(1)),
                                 float(m.group(2)),
                                 int(m.group(3)))
        elif low.startswith("title"):
            info["title"] = body
        elif low in ("density", "electrostatic potential", "plotdata") \
                or "potential" in low or "density" in low:
            if low != "plotdata" and not low.startswith("cartesian"):
                info["quantity"] = body

    if any(g is None for g in info["grid"]):
        raise ValueError("Header unvollstaendig: #grid1/#grid2/#grid3 fehlen.")

    return info, first_data_line


def read_values(path, verbose=True):
    """Liest die 4. Spalte einer pointval-Datei als float32-Array.

    Liest in grossen Bloecken statt Zeile fuer Zeile - fuer die 1.25-GB-Dateien
    ist das etwa eine Groessenordnung schneller.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        info, first_line = parse_header(fh)

        n1 = info["grid"][0][2]
        n2 = info["grid"][1][2]
        n3 = info["grid"][2][2]
        total = n1 * n2 * n3

        values = np.empty(total, dtype=np.float32)
        filled = 0
        t0 = time.time()

        remainder = first_line
        while True:
            chunk = fh.read(CHUNK_BYTES)
            if not chunk:
                break
            chunk = remainder + chunk
            cut = chunk.rfind("\n")
            if cut == -1:                      # extrem lange Zeile - weiterlesen
                remainder = chunk
                continue
            remainder = chunk[cut + 1:]
            block = chunk[:cut]

            tokens = block.split()
            if not tokens:
                continue
            if len(tokens) % 4 != 0:
                raise ValueError(
                    f"{path}: erwartete 4 Spalten pro Zeile, "
                    f"gefunden {len(tokens)} Werte in einem Block."
                )
            arr = np.asarray(tokens, dtype=np.float32).reshape(-1, 4)[:, 3]
            n = arr.size
            if filled + n > total:
                raise ValueError(
                    f"{path}: mehr Datenpunkte als der Header angibt "
                    f"({filled + n} > {total})."
                )
            values[filled:filled + n] = arr
            filled += n

            if verbose:
                pct = 100.0 * filled / total
                sys.stdout.write(f"\r    lese {os.path.basename(path)}: "
                                 f"{pct:5.1f} %")
                sys.stdout.flush()

        # letzte, unvollstaendig gepufferte Zeile
        tokens = remainder.split()
        if tokens:
            if len(tokens) % 4 != 0:
                raise ValueError(f"{path}: letzte Zeile unvollstaendig.")
            arr = np.asarray(tokens, dtype=np.float32).reshape(-1, 4)[:, 3]
            values[filled:filled + arr.size] = arr
            filled += arr.size

    if verbose:
        sys.stdout.write(f"\r    lese {os.path.basename(path)}: 100.0 %  "
                         f"({filled:,} Punkte in {time.time() - t0:.1f} s)\n")

    if filled != total:
        raise ValueError(
            f"{path}: {filled} Werte gelesen, laut Header erwartet {total}."
        )

    # Turbomole: x laeuft am schnellsten -> Speicherlayout ist [i3, i2, i1]
    data = values.reshape(n3, n2, n1)
    # Cube: z laeuft am schnellsten -> wir wollen [i1, i2, i3]
    data = np.ascontiguousarray(np.transpose(data, (2, 1, 0)))

    return info, data


# ----------------------------------------------------------------------------
# Cube schreiben
# ----------------------------------------------------------------------------

def write_cube(path, info, data, atoms, stride=1, comment=""):
    """Schreibt ein Gaussian-Cube. Alle Laengen in Bohr."""
    if stride > 1:
        data = data[::stride, ::stride, ::stride]

    n = data.shape
    starts = np.array([info["grid"][i][0] for i in range(3)])
    deltas = np.array([info["grid"][i][1] for i in range(3)])
    vecs = info["vectors"]

    # Ursprung des ersten Voxels im kartesischen Raum
    origin = info["origin"] + sum(starts[i] * vecs[i] for i in range(3))
    # Voxelvektoren (mit Stride skaliert)
    voxel = np.array([deltas[i] * stride * vecs[i] for i in range(3)])

    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"{comment or 'Cube erzeugt mit xyzToCubeToVMDVis.py'}\n")
        fh.write(f"{info.get('quantity', '') or 'volumetric data'} | "
                 f"{info.get('title', '')} | Einheiten: Bohr\n")
        fh.write(f"{len(atoms):5d} {origin[0]:12.6f} {origin[1]:12.6f} "
                 f"{origin[2]:12.6f}\n")
        for i in range(3):
            fh.write(f"{n[i]:5d} {voxel[i][0]:12.6f} {voxel[i][1]:12.6f} "
                     f"{voxel[i][2]:12.6f}\n")
        for (znum, x, y, z) in atoms:
            fh.write(f"{znum:5d} {float(znum):12.6f} {x:12.6f} {y:12.6f} "
                     f"{z:12.6f}\n")

        # Werte: z am schnellsten, 6 pro Zeile.
        # Ein vorkompiliertes Formatmuster pro z-Reihe ist deutlich schneller
        # als eine Schleife ueber alle 15.8 Mio. Einzelwerte.
        nz = n[2]
        n_full, rest = divmod(nz, 6)
        row_fmt = ("%13.5E" * 6 + "\n") * n_full
        if rest:
            row_fmt += "%13.5E" * rest + "\n"

        flat = np.ascontiguousarray(data).reshape(n[0] * n[1], nz)
        buf = []
        for idx in range(flat.shape[0]):
            buf.append(row_fmt % tuple(flat[idx].tolist()))
            if len(buf) >= 4096:                    # gebuendelt schreiben
                fh.write("".join(buf))
                buf.clear()
        if buf:
            fh.write("".join(buf))

    return n, origin, voxel


# ----------------------------------------------------------------------------
# VMD-Skript erzeugen
#
# Der Text wird ueber @@PLATZHALTER@@ gefuellt, nicht ueber str.format oder
# string.Template: Tcl benutzt geschweifte Klammern und Dollarzeichen als
# Syntax, und beide Mechanismen wuerden genau daran ersticken.
# ----------------------------------------------------------------------------

TCL_TEMPLATE = r"""# ==============================================================
# esp.tcl - elektrostatisches Potential auf der Dichte-Isoflaeche
#
# Erzeugt von xyzToCubeToVMDVis.py  (@@STAMP@@)
# Quelle: @@SOURCES@@
#
# Start:   vmd -e esp.tcl
# oder in der Tk Console:   source esp.tcl
#
# Danach verfuegbar:
#   esp_view pi | edge | sigma      Standardansichten
#   esp_iso <wert>                  Isowert der Dichteflaeche aendern
#   esp_range <halbe_breite>        Farbskala aendern (a.u.)
#   esp_opacity <0..1>              Deckkraft der Isoflaeche aendern
#   esp_colorbar [y] [breite] ...   Farbbalken neu zeichnen / verschieben
#   esp_colorbar_off                Farbbalken entfernen
#   esp_snapshot <name>             Bild rendern
# ==============================================================

# --- Robustheit ------------------------------------------------
#
# VMD-Versionen unterscheiden sich darin, welche display- und material-Optionen
# sie kennen. Ein unbekannter Befehl wirft in Tcl einen Fehler, und ein Fehler
# bricht "source" an Ort und Stelle ab: alles danach - Materialien, Prozeduren,
# Startansicht - existiert dann einfach nicht, ohne dass im Bild ein Hinweis
# darauf zu sehen waere. Genau das ist hier passiert (display depthsort).
#
# Deshalb laeuft alles Optionale ueber _try: schlaegt es fehl, wird es gemeldet
# und uebersprungen, aber die Szene steht.
proc _try {args} {
    if {[catch {uplevel 1 $args} err]} {
        puts "! uebersprungen: $args   ($err)"
        return 0
    }
    return 1
}

# --- 0) sauberer Start ----------------------------------------
mol delete all
display resetview

# --- 1) Volumendaten laden ------------------------------------
# Die Cube-Datei bringt die Atome selbst mit - in Bohr und aus derselben
# Quelle wie das Gitter. Eine separate Strukturdatei ist in VMD deshalb
# nicht noetig, und Atome und Gitter koennen nicht gegeneinander verrutschen.
set espmol [mol new @@RHO_CUBE@@ type cube waitfor all]
mol addfile @@ESP_CUBE@@ type cube waitfor all

# Die zweite Cube-Datei bringt einen identischen Koordinatensatz mit, den VMD
# als weiteres Frame anhaengt. Harmlos, aber verwirrend - weg damit.
if {[molinfo $espmol get numframes] > 1} {
    animate delete beg 1 end [expr {[molinfo $espmol get numframes] - 1}] $espmol
}

# Volume 0 = @@RHO_LABEL@@
# Volume 1 = @@ESP_LABEL@@
set VOL_RHO @@VOL_RHO@@
set VOL_ESP @@VOL_ESP@@

set ISO   @@ISO@@
set RANGE @@RANGE@@
set SCALE @@SCALE@@

# --- 2) Darstellungen -----------------------------------------
#
# DIE REIHENFOLGE IST NICHT KOSMETIK.
#
# VMD zeichnet die Reps in der Reihenfolge ihrer Nummer und schreibt dabei auch
# fuer transparente Flaechen in den Tiefenpuffer. Steht die Isoflaeche auf Rep 0,
# wird sie zuerst gezeichnet, und das anschliessend gezeichnete Geruest faellt
# ueberall dort aus dem Tiefentest, wo die Flaeche davor liegt - also ueberall.
# Ergebnis: eine sichtbar durchscheinende Flaeche, hinter der trotzdem nichts zu
# sehen ist, egal wie klein man die Deckkraft macht. Nur da, wo die nahe
# Clipping-Ebene die Flaeche aufschneidet, taucht das Geruest auf.
#
# Deshalb: erst das opake Geruest, dann die transparente Flaeche.
mol delrep 0 $espmol

set REP_MOL  0
set REP_SURF 1

# Rep 0: Molekuelgeruest. Bindungen erkennt VMD ueber Abstaende - fuer
#        organische Molekuele zuverlaessig, Bindungsordnungen braucht die
#        Darstellung nicht.
#
#        Kohlenstoff auf Grau, wie in der PyMOL-Pipeline (dort util.cnc mit
#        grey70). VMD faerbt C sonst cyan, was neben der roten Isoflaeche
#        aussieht wie ein zweites Signal.
_try color Name C gray
mol representation Licorice 0.150000 24.000000 24.000000
mol color Name
mol selection {all}
mol material Opaque
mol addrep $espmol

# Rep 1: Isoflaeche der Elektronendichte bei rho = ISO a.u.
#        (Politzer/Murray-Konvention fuer die "Molekueloberflaeche"),
#        eingefaerbt nach dem zweiten Volumendatensatz.
#
#        Das ist der Kern der ganzen Sache und das VMD-Aequivalent zu
#        ramp_new + set surface_color in PyMOL: VMD kennt kein Rampenobjekt,
#        sondern faerbt eine Isoflaeche direkt nach den Werten eines anderen
#        Gitters, das in DERSELBEN Molekuel-ID liegt.
#
#        Isosurface <iso> <volID> <show> <draw> <step> <size>
#          show 0 = nur Flaeche (1 = Flaeche+Box, 2 = nur Box)
#          draw 0 = solide      (1 = Drahtgitter, 2 = Punkte)
mol representation Isosurface $ISO $VOL_RHO 0 0 1 1
mol color Volume $VOL_ESP
mol selection {all}
mol material Transparent
mol addrep $espmol

# Farbskala: @@RANGE_NEG@@ .. +@@RANGE@@ a.u.  =  @@KJ_NEG@@ .. +@@KJ_POS@@ kJ/(mol*e)
# Rot = negativ, blau = positiv (RWB laeuft von Rot nach Blau).
# Der Wert ist bewusst ein Parameter und wird NICHT aus den Daten bestimmt -
# nur eine von aussen vorgegebene, identische Skala macht die Bilder aus
# PyMOL und VMD vergleichbar.
color scale method @@COLORSCALE@@
color scale midpoint 0.5
mol scaleminmax $espmol $REP_SURF @@RANGE_NEG@@ @@RANGE@@

# --- 3) Ansichten und Befehle ---------------------------------------------

# Auf die ATOME zentrieren, nicht auf die Gitterbox. VMDs resetview passt die
# Ansicht an die Ausdehnung aller Reps an, und die Isoflaeche traegt die
# Ausdehnung des gesamten Gitters (hier 30 x 30 x 30 Bohr). Das Molekuel sitzt
# darin fast nie mittig und landet sonst klein und schief am Bildrand.
proc esp_center {} {
    global espmol SCALE
    set sel [atomselect $espmol all]
    set c [measure center $sel]
    molinfo $espmol set center_matrix [list [transoffset [vecscale -1.0 $c]]]
    molinfo $espmol set rotate_matrix [list [transidentity]]
    molinfo $espmol set global_matrix [list [transidentity]]
    molinfo $espmol set scale_matrix  [list [transidentity]]
    $sel delete
    scale to $SCALE
}

# Die drei Standardansichten der PyMOL-Pipeline.
#
# ACHTUNG, Annahme: Die Rotationen gelten fuer ein planares, entlang y
# ausgerichtetes Molekuel in der xy-Ebene, mit dem Halogen bei y = 0 und dem
# Ring bei y > 0 - so, wie die Turbomole-Optimierung die Halogenbenzole hier
# ablegt. Fuer andere Molekuele von Hand drehen und esp_snapshot benutzen.
#
# Der Drehsinn bei sigma ist nicht beliebig: falsch herum zeigt das Bild die
# gegenueberliegende Seite des Molekuels, und weil dort ebenfalls eine runde,
# gefaerbte Flaeche sitzt, faellt der Fehler nicht auf.
proc esp_view {which} {
    esp_center
    switch -- $which {
        pi      { }
        edge    { rotate y by 90 }
        sigma   {
            # -90, nicht +90. Das Brom sitzt bei y = 0, der Ring bei y > 0;
            # das sigma-Loch liegt auf der Verlaengerung der C-Br-Achse, also
            # bei y < 0. Mit +90 kippt der RING zur Kamera und man schaut auf
            # das para-H statt auf das sigma-Loch - das Bild sieht plausibel
            # aus und zeigt die falsche Seite.
            rotate x by -90
            puts "Hinweis: In der Achsenansicht stapeln sich die transparenten"
            puts "         Lagen der Isoflaeche. Fuer ein sauberes Bild"
            puts "         esp_opacity 1.0 oder esp_snapshot (Tachyon)."
        }
        default { puts "esp_view: pi | edge | sigma" ; return }
    }
    # Der Balken haengt an Zentrierung und Skalierung des Molekuels und wird
    # von "mol fix" bewusst nicht mitbewegt - also nach jedem Ansichtswechsel
    # neu zeichnen, damit er nicht zurueckbleibt.
    global cbmol cbopts
    if {[info exists cbmol] && [lsearch [molinfo list] $cbmol] >= 0} {
        eval esp_colorbar $cbopts
    }
    puts "Ansicht: $which"
}

proc esp_iso {value} {
    global espmol VOL_RHO REP_SURF ISO
    set ISO $value
    mol modstyle $REP_SURF $espmol Isosurface $value $VOL_RHO 0 0 1 1
    puts "Isowert: rho = $value a.u."
}

proc esp_range {half} {
    global espmol REP_SURF RANGE cbmol
    set RANGE $half
    mol scaleminmax $espmol $REP_SURF [expr {-1.0 * $half}] $half
    if {[info exists cbmol] && [lsearch [molinfo list] $cbmol] >= 0} {
        esp_colorbar
    }
    puts "Farbskala: +/- $half a.u."
}

# --- Farbskala als Balken im Bild -----------------------------
#
# VMD hat keine fertige Legende. Der Balken wird als Grafikprimitive in eine
# EIGENE, leere Molekuel-ID gezeichnet und mit "mol fix" vom Maus-Transform
# abgekoppelt - sonst kippt er beim Drehen des Molekuels mit.
#
# Positioniert wird er NICHT in geratenen Bildschirmkoordinaten, sondern
# relativ zum Molekuel: in Angstroem, unterhalb von dessen Huellkugel, und mit
# derselben Zentrierung und Skalierung wie das Molekuel selbst. Damit landet er
# unabhaengig von Fenstergroesse und Molekuelgroesse immer knapp unter der
# Isoflaeche. Erst dadurch ist er ueberhaupt zuverlaessig im Bild.
#
#   gap        Abstand zur Huellkugel in Angstroem
#   barheight  Hoehe des Balkens in Angstroem
#   textsize   Schriftgroesse (bildschirmbezogen, nicht in Angstroem)
#
# Nach einem manuellen Zoom sitzt er nicht mehr passend - dann einfach
# esp_colorbar erneut aufrufen.
proc esp_colorbar {{gap 0.5} {barheight 0.7} {textsize 0.8}} {
    global espmol RANGE cbmol cbopts
    set cbopts [list $gap $barheight $textsize]

    set sel [atomselect $espmol all]
    set c  [measure center $sel]
    set mm [measure minmax $sel]
    $sel delete

    set ex [expr {[lindex [lindex $mm 1] 0] - [lindex [lindex $mm 0] 0]}]
    set ey [expr {[lindex [lindex $mm 1] 1] - [lindex [lindex $mm 0] 1]}]
    set ez [expr {[lindex [lindex $mm 1] 2] - [lindex [lindex $mm 0] 2]}]
    # Halbe Raumdiagonale der Atomhuelle, plus rund 1.8 A fuer die Isoflaeche.
    # Ueber die Huellkugel statt ueber die Box, damit der Balken auch nach dem
    # Drehen in die sigma- oder edge-Ansicht nicht ins Molekuel rutscht.
    set r [expr {0.5 * sqrt($ex*$ex + $ey*$ey + $ez*$ez) + 1.8}]

    set cx [lindex $c 0]
    set cy [lindex $c 1]
    set cz [lindex $c 2]
    set hw [expr {0.85 * $r}]
    set y0 [expr {$cy - $r - $gap}]
    set y1 [expr {$y0 + $barheight}]
    set x0 [expr {$cx - $hw}]

    if {[info exists cbmol] && [lsearch [molinfo list] $cbmol] >= 0} {
        graphics $cbmol delete all
    } else {
        mol new
        set cbmol [molinfo top]
        mol rename $cbmol colorbar
    }

    # Ohne Beleuchtung zeichnen: der Balken soll die Farben der Skala zeigen,
    # nicht eine schattierte Version davon.
    graphics $cbmol materials off

    # Die Farbskala belegt die Farb-IDs oberhalb der benannten Farben.
    set c0 [colorinfo num]
    set c1 [expr {[colorinfo max] - 1}]

    set nseg 64
    set dx [expr {2.0 * $hw / $nseg}]
    for {set i 0} {$i < $nseg} {incr i} {
        set cid [expr {int($c0 + double($i) / ($nseg - 1) * ($c1 - $c0))}]
        graphics $cbmol color $cid
        set xa [expr {$x0 + $i * $dx}]
        set xb [expr {$xa + $dx}]
        graphics $cbmol triangle "$xa $y0 $cz" "$xb $y0 $cz" "$xb $y1 $cz"
        graphics $cbmol triangle "$xa $y0 $cz" "$xb $y1 $cz" "$xa $y1 $cz"
    }

    graphics $cbmol color black
    set ty [expr {$y0 - 1.1 * $barheight}]
    graphics $cbmol text "$x0 $ty $cz" \
        [format "%.3f" [expr {-1.0 * $RANGE}]] size $textsize thickness 2
    graphics $cbmol text "[expr {$cx - 0.3}] $ty $cz" "0" \
        size $textsize thickness 2
    graphics $cbmol text "[expr {$cx + 0.45 * $hw}] $ty $cz" \
        [format "+%.3f" $RANGE] size $textsize thickness 2
    graphics $cbmol text "[expr {$cx - 1.2}] [expr {$y1 + 0.4 * $barheight}] $cz" \
        "ESP / a.u." size $textsize thickness 2

    # Gleiche Zentrierung und Skalierung wie das Molekuel, aber ohne dessen
    # Rotation - der Balken bleibt so beim Drehen waagerecht und an Ort und
    # Stelle. mol fix koppelt ihn zusaetzlich von der Maus ab; danach das
    # Molekuel wieder nach oben, sonst zielen Mausaktionen auf den Balken.
    molinfo $cbmol set center_matrix [list [transoffset [vecscale -1.0 $c]]]
    molinfo $cbmol set rotate_matrix [list [transidentity]]
    molinfo $cbmol set global_matrix [molinfo $espmol get global_matrix]
    molinfo $cbmol set scale_matrix  [molinfo $espmol get scale_matrix]
    mol fix $cbmol
    mol top $espmol
    puts "Farbbalken: Molekuel-ID $cbmol, y = [format %.2f $y0] A"
}

proc esp_colorbar_off {} {
    global cbmol
    if {[info exists cbmol] && [lsearch [molinfo list] $cbmol] >= 0} {
        mol delete $cbmol
        unset cbmol
    }
}

# Deckkraft der Isoflaeche live aendern, ohne die Cubes neu zu erzeugen.
# 0.30 = Standard, Geruest scheint durch. Ab etwa 0.5 verschwindet es (zwei
# Lagen!), unter etwa 0.15 wird die Farbe zu blass, um sie noch abzulesen.
proc esp_opacity {value} {
    material change opacity Transparent $value
    puts "Deckkraft: $value"
}

proc esp_snapshot {name} {
    # Strahlverfolgt, mit echter Transparenz und weichen Kanten.
    # Schnelle Alternative (exakt das Fensterbild):  render snapshot $name.png
    render TachyonInternal $name.tga
    puts "geschrieben: $name.tga"
}

# --- 4) Anzeige -----------------------------------------------
_try display projection Orthographic
_try display depthcue off
_try display shadows off
_try display culling off
_try axes location off
_try color Display Background white

# Nahe Clipping-Ebene fast auf null. Der Default (0.5) schneidet beim
# Hineinzoomen eine Scheibe aus der Isoflaeche heraus - man sieht dann ins
# Molekuel hinein und haelt es leicht faelschlich fuer Transparenz.
_try display nearclip set 0.010000

# Transparente Flaechen muessen von hinten nach vorne gezeichnet werden,
# sonst blendet VMD in willkuerlicher Reihenfolge und die Flaeche wirkt dicht.
_try display depthsort on

# Echte Transparenz statt Rasterpunkten. Der Default-Rendermode loest
# Transparenz durch Weglassen von Pixeln ("screen door") - das sieht aus wie
# ein kaputter Drucker. Ueber _try, weil aeltere Treiber kein GLSL koennen; dann
# bleibt es beim Default, und die finalen Bilder rendert ohnehin Tachyon.
_try display rendermode GLSL

# Deckkraft: Der Blick geht durch ZWEI Lagen der geschlossenen Isoflaeche
# (vorne rein, hinten raus). Bei einer Deckkraft a bleibt hinter beiden Lagen
# nur (1-a)^2 uebrig - bei a = 0.6 sind das 16 %, und das Geruest verschwindet,
# obwohl die Zahl nach "halb durchsichtig" klingt. Deshalb sind Werte um 0.3
# noetig, um in VMD den Eindruck zu bekommen, den PyMOL bei transparency 0.15
# liefert. Gleiche Zahl heisst hier nicht gleiches Bild.
_try material change opacity  Transparent @@OPACITY@@
_try material change diffuse  Transparent 0.750000

# Glanzlicht klein halten: ein kraeftiger Speculareffekt legt sich als milchiger
# Schleier ueber die Flaeche und frisst genau den Durchblick wieder auf, den die
# Deckkraft freigibt.
_try material change specular Transparent 0.100000
_try material change shininess Transparent 0.300000

# --- 5) Startansicht ------------------------------------------
esp_view pi
@@COLORBAR_CALL@@

puts ""
puts "esp.tcl geladen."
puts "  Isoflaeche  rho = $ISO a.u., eingefaerbt nach dem ESP"
puts "  Farbskala   +/- $RANGE a.u. (rot = negativ, blau = positiv)"
puts "  Deckkraft   @@OPACITY@@"
puts "  Befehle     esp_view pi|edge|sigma, esp_iso, esp_range, esp_opacity,"
puts "              esp_colorbar, esp_colorbar_off, esp_snapshot"
puts ""
"""


def classify_cube(path):
    """'density' / 'esp' / None - anhand der zweiten Kopfzeile des Cubes.

    Die schreibt dieses Skript selbst; sie traegt die ``#quantity``-Angabe aus
    der pointval-Datei. Damit laesst sich eine Szene aus fertigen Cubes neu
    bauen, ohne die 1.25-GB-Gitter noch einmal anzufassen (``--tcl-only``).
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        fh.readline()
        second = fh.readline().lower()
    if "potential" in second:
        return "esp"
    if "density" in second:
        return "density"
    return None


def write_vmd_script(path, rho_cube, esp_cube, rho_label, esp_label,
                     esp_range, iso=0.001, opacity=0.50, scale=0.12,
                     colorscale="RWB", sources="", colorbar=True):
    """Schreibt das ``esp.tcl``.

    ``rho_cube`` wird als erstes geladen und ist damit Volume 0, ``esp_cube``
    als zweites und damit Volume 1. Die Reihenfolge steht hier fest und wird
    nicht geraten - sie folgt aus der ``#quantity``-Zeile der pointval-Dateien.
    """
    repl = {
        "@@STAMP@@": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "@@SOURCES@@": sources,
        "@@RHO_CUBE@@": rho_cube,
        "@@ESP_CUBE@@": esp_cube,
        "@@RHO_LABEL@@": rho_label,
        "@@ESP_LABEL@@": esp_label,
        "@@VOL_RHO@@": "0",
        "@@VOL_ESP@@": "1",
        "@@ISO@@": f"{iso:g}",
        "@@RANGE@@": f"{esp_range:.4f}",
        "@@RANGE_NEG@@": f"{-esp_range:.4f}",
        "@@KJ_NEG@@": f"{-esp_range * HARTREE_TO_KJ:.0f}",
        "@@KJ_POS@@": f"{esp_range * HARTREE_TO_KJ:.0f}",
        "@@OPACITY@@": f"{opacity:.6f}",
        "@@SCALE@@": f"{scale:g}",
        "@@COLORSCALE@@": colorscale,
        "@@COLORBAR_CALL@@": (
            "# Fehler beim Farbbalken duerfen den Rest der Szene nicht\n"
            "# mitreissen - Tcl bricht ein source sonst an Ort und Stelle ab.\n"
            "if {[catch {esp_colorbar} cberr]} {\n"
            "    puts \"! Farbbalken nicht gezeichnet: $cberr\"\n"
            "}" if colorbar else
                              "# Farbbalken abgeschaltet (--no-colorbar); "
                              "mit  esp_colorbar  nachtraeglich einblenden"),
    }
    text = TCL_TEMPLATE
    for key, value in repl.items():
        text = text.replace(key, str(value))

    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


# ----------------------------------------------------------------------------
# Hauptprogramm
# ----------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        description="Turbomole-pointval-Gitterdateien (td.xyz, tp.xyz) nach "
                    "Gaussian-Cube konvertieren und eine fertige VMD-Szene "
                    "schreiben.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Beispiel:\n"
               "  python xyzToCubeToVMDVis.py --struct brombenzol_aro_opti.mol "
               "td.xyz tp.xyz\n"
               "  vmd -e esp.tcl\n",
    )
    p.add_argument("grids", nargs="+",
                   help="Turbomole-Gitterdateien (z.B. td.xyz tp.xyz)")
    p.add_argument("--struct", "-s", default=None,
                   help="Strukturdatei fuer den Cube-Header: .xyz, .mol "
                        "oder .sdf (entfaellt bei --tcl-only)")
    p.add_argument("--struct-unit", choices=["angstrom", "bohr"],
                   default="angstrom",
                   help="Einheit der Strukturdatei (Standard: angstrom); gilt "
                        "nur fuer .xyz - .mol/.sdf sind immer Angstrom")
    p.add_argument("--outdir", "-o", default=None,
                   help="Ausgabeverzeichnis (Standard: neben der Eingabe)")
    p.add_argument("--stride", type=int, default=1,
                   help="Nur jeden n-ten Gitterpunkt schreiben "
                        "(2 => 8x kleinere Datei, Standard: 1)")
    p.add_argument("--quiet", "-q", action="store_true")

    # Diese Optionen aendern NICHTS an den Cube-Dateien - sie beschreiben nur
    # die Szene, die nebenbei geschrieben wird.
    g = p.add_argument_group("VMD-Szene (esp.tcl)")
    g.add_argument("--no-vmd", action="store_true",
                   help="Nur die Cube-Dateien schreiben, kein esp.tcl")
    g.add_argument("--tcl-only", action="store_true",
                   help="Nur esp.tcl neu schreiben, die vorhandenen Cubes "
                        "unangetastet lassen. Fuer jede Aenderung an "
                        "Deckkraft, Farbskala oder Isowert - das Umrechnen "
                        "der 1.25-GB-Gitter dauert Minuten, das hier "
                        "Millisekunden. --struct wird dann nicht gebraucht.")
    g.add_argument("--esp-range", type=float, default=0.035,
                   help="Halbe Breite der ESP-Farbskala in a.u. "
                        "(Standard: 0.035). Wird bewusst NICHT aus den Daten "
                        "bestimmt: Den Wert liefert die PyMOL-Pipeline in "
                        "<molekuel>_settings.txt, Zeile 'Farbskala'. Nur eine "
                        "identische Skala macht die Bilder vergleichbar.")
    g.add_argument("--iso", type=float, default=0.001,
                   help="Isowert der Dichteflaeche in der Szene "
                        "(Standard: 0.001 a.u.)")
    g.add_argument("--opacity", type=float, default=0.50,
                   help="Deckkraft der Isoflaeche, 0..1 (Standard: 0.50). "
                        "Der Blick geht durch zwei Lagen der Flaeche, es "
                        "bleibt also (1-a)^2 uebrig - hoehere Werte geben "
                        "kraeftigere Farben, blassere ein deutlicheres "
                        "Geruest. In VMD live aenderbar mit esp_opacity.")
    g.add_argument("--no-colorbar", action="store_true",
                   help="Farbbalken beim Start nicht einblenden (in VMD "
                        "jederzeit mit esp_colorbar nachholbar)")
    g.add_argument("--scale", type=float, default=0.10,
                   help="Zoomfaktor der Startansicht (Standard: 0.10 - etwas "
                        "herausgezoomt, damit unter dem Molekuel Platz fuer "
                        "den Farbbalken bleibt). 0.08 bis 0.18 ist der "
                        "brauchbare Bereich.")
    g.add_argument("--color-scale", default="RWB",
                   help="VMD-Farbskala (Standard: RWB = rot-weiss-blau, also "
                        "rot negativ). BWR dreht sie um.")
    args = p.parse_args(argv)

    verbose = not args.quiet

    if args.tcl_only and args.no_vmd:
        p.error("--tcl-only und --no-vmd schliessen sich aus - dann bliebe "
                "nichts zu tun.")
    if not args.tcl_only and not args.struct:
        p.error("--struct wird gebraucht: der Cube-Header verlangt einen "
                "Atomblock, und die pointval-Datei enthaelt keine Atome.")

    if verbose:
        print("=" * 70)
        print("xyzToCubeToVMDVis.py - Turbomole pointval -> Cube -> VMD")
        print("=" * 70)

    written = {}
    labels = {}

    # --- Nur die Szene neu schreiben ------------------------------------
    if args.tcl_only:
        for gpath in args.grids:
            base = os.path.splitext(os.path.basename(gpath))[0]
            outdir_g = args.outdir or os.path.dirname(os.path.abspath(gpath))
            cube = os.path.join(outdir_g, base + ".cube")
            if not os.path.exists(cube):
                p.error(f"{cube} fehlt - ohne Cube-Datei kann --tcl-only "
                        f"nichts einordnen. Einmal ohne --tcl-only laufen "
                        f"lassen.")
            kind = classify_cube(cube)
            if kind is None:
                print(f"    ! {os.path.basename(cube)}: Kopfzeile nennt weder "
                      f"Dichte noch Potential - uebersprungen.", file=sys.stderr)
                continue
            written[kind] = cube
            labels[kind] = ("Elektronendichte" if kind == "density"
                            else "elektrostatisches Potential")
            if verbose:
                print(f"[1] gefunden: {os.path.basename(cube)} -> {kind}")
    else:
        atoms = read_structure(args.struct, unit=args.struct_unit)
        if verbose:
            print(f"[1] Struktur: {args.struct} -> {len(atoms)} Atome "
                  f"(eingelesen als {args.struct_unit}, gespeichert als Bohr)")

    for gpath in ([] if args.tcl_only else args.grids):
        if verbose:
            print(f"[2] Gitterdatei: {gpath}")
        info, data = read_values(gpath, verbose=verbose)

        n1, n2, n3 = (info["grid"][i][2] for i in range(3))
        if verbose:
            print(f"    Gitter {n1} x {n2} x {n3}, "
                  f"delta = {info['grid'][0][1]} Bohr, "
                  f"Groesse = '{info['quantity'] or 'unbekannt'}'")
            print(f"    Wertebereich: {data.min():+.6g} .. {data.max():+.6g}")

        base = os.path.splitext(os.path.basename(gpath))[0]
        outdir = args.outdir or os.path.dirname(os.path.abspath(gpath))
        os.makedirs(outdir, exist_ok=True)
        outpath = os.path.join(outdir, base + ".cube")

        shape, origin, voxel = write_cube(
            outpath, info, data, atoms, stride=args.stride,
            comment=f"{info['quantity'] or base} - konvertiert aus "
                    f"{os.path.basename(gpath)}",
        )
        if verbose:
            mb = os.path.getsize(outpath) / 1024 ** 2
            print(f"    -> {outpath}  ({shape[0]}x{shape[1]}x{shape[2]}, "
                  f"{mb:.1f} MB)")

        q = (info["quantity"] or "").lower()
        if "potential" in q:
            written["esp"] = outpath
            labels["esp"] = info["quantity"] or base
        elif "density" in q:
            written["density"] = outpath
            labels["density"] = info["quantity"] or base
        else:
            written.setdefault("other", []).append(outpath)
            if verbose:
                print(f"    ! '{info['quantity']}' ist weder Dichte noch "
                      f"Potential - fuer esp.tcl unbenutzt.")

        # Das Gitter kann mehrere hundert MB belegen; hier wird es nicht mehr
        # gebraucht, weil dieses Projekt bewusst nichts nachrechnet.
        del data

    if not args.no_vmd:
        outdir = args.outdir or os.path.dirname(os.path.abspath(args.grids[0]))
        tcl = os.path.join(outdir, "esp.tcl")

        if "esp" not in written or "density" not in written:
            fehlt = "Dichte (td)" if "density" not in written else "Potential (tp)"
            print(f"    ! {fehlt} fehlt - esp.tcl wird uebersprungen. "
                  f"Fuer die Szene werden beide Gitter gebraucht: die Dichte "
                  f"liefert die Flaeche, das Potential die Farbe.",
                  file=sys.stderr)
        else:
            write_vmd_script(
                tcl,
                rho_cube=os.path.basename(written["density"]),
                esp_cube=os.path.basename(written["esp"]),
                rho_label=labels.get("density", "Elektronendichte"),
                esp_label=labels.get("esp", "elektrostatisches Potential"),
                esp_range=args.esp_range,
                iso=args.iso,
                opacity=args.opacity,
                scale=args.scale,
                colorscale=args.color_scale,
                sources=", ".join(os.path.basename(g) for g in args.grids),
                colorbar=not args.no_colorbar,
            )
            if verbose:
                print(f"[3] VMD-Szene: {tcl}")
                print(f"    Farbskala:  +/- {args.esp_range:.4f} a.u. "
                      f"(Parameter, nicht aus den Daten bestimmt)")
                print(f"    Start mit:  vmd -e {os.path.basename(tcl)}")

    if verbose:
        print("Fertig.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ansi.py
=======

Minimale Farbunterstuetzung fuer die Konsolenausgabe. Keine Abhaengigkeiten.

Farben werden automatisch abgeschaltet, wenn die Ausgabe kein Terminal ist -
also beim Umleiten in eine Datei oder beim Weiterreichen durch eine Pipe.
Sonst landen Steuerzeichen wie ``\\x1b[32m`` in den Logdateien.

Zusaetzlich respektiert das Modul die Konvention ``NO_COLOR`` (siehe
https://no-color.org): ist die Umgebungsvariable gesetzt, bleibt alles farblos.
Mit ``FORCE_COLOR`` laesst sich das Gegenteil erzwingen.

Unter Windows muss die Verarbeitung von ANSI-Sequenzen einmal aktiviert werden
(``ENABLE_VIRTUAL_TERMINAL_PROCESSING``); das erledigt ``_enable_windows_vt``
ueber die Win32-API, ohne colorama als zusaetzliches Paket.
"""

from __future__ import annotations

import os
import re
import sys

RESET = "\033[0m"
BOLD = "\033[1m"
GREEN = "\033[32m"
CYAN = "\033[36m"
YELLOW = "\033[33m"

def _enable_windows_vt() -> bool:
    """Schaltet die ANSI-Verarbeitung der Windows-Konsole frei."""
    try:
        import ctypes

        kernel32 = ctypes.windll.kernel32
        handle = kernel32.GetStdHandle(-11)          # STD_OUTPUT_HANDLE
        mode = ctypes.c_uint32()
        if not kernel32.GetConsoleMode(handle, ctypes.byref(mode)):
            return False
        ENABLE_VIRTUAL_TERMINAL_PROCESSING = 0x0004
        return bool(kernel32.SetConsoleMode(
            handle, mode.value | ENABLE_VIRTUAL_TERMINAL_PROCESSING))
    except Exception:
        return False

def _supported() -> bool:
    if os.environ.get("NO_COLOR"):
        return False
    if os.environ.get("FORCE_COLOR"):
        return True
    if not (hasattr(sys.stdout, "isatty") and sys.stdout.isatty()):
        return False
    if os.name == "nt":
        return _enable_windows_vt()
    return True

ENABLED = _supported()

def disable():
    """Schaltet die Farbausgabe zur Laufzeit ab (fuer --no-color)."""
    global ENABLED
    ENABLED = False

def paint(text, code) -> str:
    """Faerbt ``text``, oder gibt ihn unveraendert zurueck."""
    if not ENABLED or not text:
        return text
    return f"{code}{text}{RESET}"

def green(text):
    return paint(text, GREEN)

def cyan(text):
    return paint(text, CYAN)

def yellow(text):
    return paint(text, YELLOW)

def bold(text):
    return paint(text, BOLD)

# ----------------------------------------------------------------------------
# Chemie-spezifisch
# ----------------------------------------------------------------------------

HALOGEN_SYMBOLS = ("F", "Cl", "Br", "I", "At")

_LABEL = re.compile(r"^([A-Za-z]+)(\d*)$")

def atom_label(label: str) -> str:
    """Faerbt das Elementsymbol eines Atomlabels, wenn es ein Halogen ist.

    ``"Cl12"`` -> tuerkises ``Cl`` plus normales ``12``. Die laufende Nummer
    bleibt ungefaerbt, damit das Symbol ins Auge springt und nicht der Index.
    Nicht-Halogene wie ``"H5"`` bleiben unveraendert.
    """
    if not label:
        return label
    m = _LABEL.match(label)
    if not m:
        return label
    symbol, number = m.group(1), m.group(2)
    if symbol in HALOGEN_SYMBOLS:
        return cyan(symbol) + number
    return label

def element(symbol: str) -> str:
    """Faerbt ein blankes Elementsymbol, wenn es ein Halogen ist."""
    return cyan(symbol) if symbol in HALOGEN_SYMBOLS else symbol

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
constants.py
============

Physikalische Konstanten der ESP-Pipeline, an EINER Stelle.

Warum ein eigenes Modul: der Umrechnungsfaktor Angstrom -> Bohr stand vorher
in jedem Skript einzeln. Solche Duplikate sind harmlos, bis eines davon
angefasst wird - dann rechnen zwei Teile derselben Pipeline mit verschiedenen
Werten, und das Ergebnis ist eine um Promille verschobene Geometrie, die man
im Bild nicht sieht und in den Zahlen fuer echt haelt.

Bewusst ohne jede Abhaengigkeit (auch ohne numpy), damit dieses Modul von
ueberall importiert werden kann, ohne eine Umgebung vorauszusetzen.
"""

from __future__ import annotations

# CODATA 2018: 1 Bohr = 0.529177210903 Angstrom
BOHR_PER_ANGSTROM = 1.8897259886
ANGSTROM_PER_BOHR = 1.0 / BOHR_PER_ANGSTROM

# Energieumrechnung fuer die ESP-Werte (Hartree pro Elementarladung)
HARTREE_TO_KCAL = 627.5095
HARTREE_TO_KJ = 2625.4996

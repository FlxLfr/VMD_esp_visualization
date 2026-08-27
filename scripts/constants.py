#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
constants.py
============

The physical constants of the ESP pipeline, in ONE place.

Why a module of its own: the Angstrom -> Bohr conversion factor used to sit in
every script separately. Duplicates like that are harmless until one of them is
touched - then two parts of the same pipeline compute with different values,
and the result is a geometry shifted by a fraction of a per mille, invisible in
the picture and taken for real in the numbers.

Deliberately without any dependency (numpy included), so that this module can
be imported from anywhere without presupposing an environment.
"""

from __future__ import annotations

# CODATA 2018: 1 Bohr = 0.529177210903 Angstrom
BOHR_PER_ANGSTROM = 1.8897259886
ANGSTROM_PER_BOHR = 1.0 / BOHR_PER_ANGSTROM

# Energy conversion for the ESP values (Hartree per elementary charge)
HARTREE_TO_KCAL = 627.5095
HARTREE_TO_KJ = 2625.4996

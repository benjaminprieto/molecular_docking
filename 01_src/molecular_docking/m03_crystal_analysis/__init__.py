"""
m03_crystal_analysis — Crystal Structure Interaction Analysis
===============================================================
Analyze protein-ligand interactions from crystallographic data
using PLIP (Protein-Ligand Interaction Profiler).

This module operates on the reference complex (receptor + co-crystallized
ligand) with correct protonation states. Results feed into downstream
modules (m06 Pharmit, docking validation, etc.).

Modules:
    plip_interaction_analysis  03a — PLIP interaction detection + JSON output

Location: 01_src/molecular_docking/m03_crystal_analysis/
Project: molecular_docking
Module: m03
Version: 1.0

Dependencies: plip (conda install -c conda-forge plip)
"""
__version__ = "1.0.0"
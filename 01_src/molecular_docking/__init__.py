"""
molecular_docking - Multi-Engine Docking Pipeline
===================================================
Generic pipeline for molecular docking with DOCK6 and GNINA.
Produces outputs compatible with dock2profile.

Modules:
    m00_preparation  - Parse molecules, prepare receptor, binding site (shared)
    m01_docking      - DOCK6 engine (antechamber, grids, docking, scores)
    m02_gnina        - GNINA engine (Vina + CNN scoring, docking, scores)
"""
__version__ = "2.1.0"
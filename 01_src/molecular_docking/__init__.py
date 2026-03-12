"""
molecular_docking - DOCK6 Docking Pipeline
============================================
Generic pipeline for molecular docking with DOCK6.
Produces outputs compatible with dock2profile.

Modules:
    m00_preparation  - Parse molecules, prepare receptor & ligands
    m01_docking      - Grid generation & DOCK6 execution
    m02_collection   - Score collection & Excel generation
"""
__version__ = "1.0.0"

"""
molecular_docking - Virtual Screening Pipeline
===============================================
Virtual screening pipeline with DOCK6 and GNINA engines.
Produces outputs compatible with dock2profile.

Modules:
    m00_preparation    - Parse molecules, prepare receptor, binding site (shared)
    m01_docking        - DOCK6 engine (antechamber, grids, docking, scores)
    m02_gnina          - GNINA engine (Vina + CNN scoring, docking, scores)
    m04_dock6_analysis - DOCK6 results analysis (ranking, footprints, modes, contacts)
    m05_gnina_analysis - GNINA results analysis (fragments, clustering, hotspots, export)
"""
__version__ = "3.0.0"
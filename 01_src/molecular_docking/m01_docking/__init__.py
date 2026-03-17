"""
m01_docking - DOCK6 Engine
============================
Complete DOCK6 docking pipeline.

Modules:
    antechamber_preparation  01a — SDF → mol2 (AM1-BCC charges, GAFF2)
    grid_generation          01b — DMS → spheres → grids
    dock6_runner             01c — dock6 per molecule
    score_collector          01d — parse scores → Excel
"""
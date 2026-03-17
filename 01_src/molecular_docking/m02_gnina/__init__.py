"""
m02_gnina - GNINA Docking Engine
===================================
GNINA = Vina + CNN scoring in a single run.
Produces: Vina affinity + CNN score + CNN affinity per molecule.

Modules:
    gnina_preparation     02a — binding box + inputs JSON (no format conversion needed)
    gnina_runner          02b — GNINA docking per molecule (dock or rescore mode)
    gnina_score_collector 02c — parse scores → CSV/Excel with CNN columns
"""

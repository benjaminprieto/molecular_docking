"""
m05_gnina_analysis — Pose Clustering & Fragment Analysis
======================================================
Multi-resolution analysis of docking poses by rigid fragment decomposition.

Pipeline:
    05a parse_and_fragment     — Parse 3D coords + identify rigid fragments
    05b fragment_clustering    — DBSCAN clustering per fragment per molecule
    05c score_decomposition    — Exact energy per fragment (atom_terms) + CNN regression
    05d binding_site_hotspots  — Cross-molecule hotspot detection
    05e structure_export       — mol2 representatives + ChimeraX scripts
    05f contact_mapping        — Receptor contacts per fragment cluster
    05g campaign_report        — Consolidated HTML report with figures and tables

Location: 01_src/molecular_docking/m05_gnina_analysis/
Project: molecular_docking
Module: m05
Version: 2.1
"""

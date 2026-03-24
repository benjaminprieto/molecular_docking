"""
m04_dock6_analysis - DOCK6 Post-Docking Analysis
===================================================
Residue-based analysis of DOCK6 docking results using footprint scoring.

DOCK6 analysis is RESIDUE-BASED (footprint), complementary to GNINA's
atom-level (fragment-based) analysis in m05.

Pipeline:
    01d footprint_rescore      — Re-score poses with fps_primary (in m01)
    04a score_ranking          — Grid_Score ranking + vdW/ES decomposition
    04b footprint_analysis     — Per-residue energy parsing + consensus (reads 01d)
    04c binding_modes          — Binding mode characterization per molecule
    04d contact_mapping        — Distance contacts + footprint cross-reference
    04e campaign_report        — HTML report with composite ranking

Reference: Balius et al. J Chem Inf Model 2011, 51(8):1942-56
"""
__version__ = "2.0.0"

"""
m04_dock6_analysis - DOCK6 Post-Docking Analysis
===================================================
Residue-based analysis of DOCK6 docking results using footprint scoring.

DOCK6 analysis is RESIDUE-BASED (footprint), complementary to GNINA's
atom-level (fragment-based) analysis in m05.

Modules:
    score_ranking          04a — Grid_Score ranking + vdW/ES decomposition
    footprint_rescoring    04b — Footprint re-scoring (calls dock6 via m01 utils)
    footprint_analysis     04b — Per-residue energy parsing + consensus
    binding_modes          04c — Binding mode characterization per molecule
    contact_mapping        04d — Distance contacts + footprint cross-reference
    campaign_report        04e — HTML report with composite ranking

Reference: Balius et al. J Chem Inf Model 2011, 51(8):1942-56
"""
__version__ = "1.0.0"

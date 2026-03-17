#!/bin/bash
# =============================================================================
# run_pipeline.sh — Run the complete docking pipeline
# =============================================================================
#
# Architecture (v2.1 — 2026-03-16):
#   00   Shared preparation (all engines)
#   01   DOCK6 engine (01a antechamber → 01b grids → 01c dock → 01d scores)
#   02   GNINA engine (02a prep → 02b dock+CNN → 02c scores)
#
# Usage:
#   bash run_pipeline.sh <campaign_id> [dock6|gnina|both]
#
# Examples:
#   bash run_pipeline.sh example_campaign              # both engines
#   bash run_pipeline.sh example_campaign dock6         # DOCK6 only
#   bash run_pipeline.sh example_campaign gnina         # GNINA only
#
# Critical modules (00a, 00b, 00c, 00d) stop the pipeline on failure.
# Batch modules (01a, 01c, 02b) continue even if some molecules fail.
#
# =============================================================================

CAMPAIGN_ID="${1:-}"
ENGINE="${2:-both}"

if [ -z "$CAMPAIGN_ID" ]; then
    echo "Usage: bash run_pipeline.sh <campaign_id> [dock6|gnina|both]"
    echo ""
    echo "Available campaigns:"
    ls -1 04_data/campaigns/ 2>/dev/null || echo "  (none)"
    exit 1
fi

CAMPAIGN="04_data/campaigns/${CAMPAIGN_ID}/campaign_config.yaml"

if [ ! -f "$CAMPAIGN" ]; then
    echo "ERROR: ${CAMPAIGN} not found"
    exit 1
fi

echo "============================================================"
echo "  molecular_docking — Full pipeline"
echo "  Campaign: ${CAMPAIGN_ID}"
echo "  Engine:   ${ENGINE}"
echo "  Started:  $(date '+%Y-%m-%d %H:%M')"
echo "============================================================"
echo ""

# =============================================================================
# 00 — SHARED PREPARATION
# =============================================================================

echo "========== 00: SHARED PREPARATION =========="
echo ""

echo "[00a] Molecule Parser"
python 02_scripts/00a_molecule_parser.py --config 03_configs/00a_molecule_parser.yaml --campaign "$CAMPAIGN"
if [ $? -ne 0 ]; then echo "FAILED at 00a"; exit 1; fi
echo ""

echo "[00b] Receptor Preparation"
python 02_scripts/00b_receptor_preparation.py --config 03_configs/00b_receptor_preparation.yaml --campaign "$CAMPAIGN"
if [ $? -ne 0 ]; then echo "FAILED at 00b"; exit 1; fi
echo ""

echo "[00c] Ionization Profiling"
python 02_scripts/00c_ionization_profiling.py --config 03_configs/00c_ionization_profiling.yaml --campaign "$CAMPAIGN"
if [ $? -ne 0 ]; then echo "FAILED at 00c"; exit 1; fi
echo ""

echo "[00d] Binding Site Definition"
python 02_scripts/00d_binding_site_definition.py --config 03_configs/00d_binding_site_definition.yaml --campaign "$CAMPAIGN"
if [ $? -ne 0 ]; then echo "FAILED at 00d"; exit 1; fi
echo ""

# =============================================================================
# 01 — DOCK6 ENGINE
# =============================================================================

if [ "$ENGINE" = "dock6" ] || [ "$ENGINE" = "both" ]; then
    echo "========== 01: DOCK6 ENGINE =========="
    echo ""

    echo "[01a] Antechamber (DOCK6-specific: SDF → mol2)"
    python 02_scripts/01a_antechamber_preparation.py --config 03_configs/01a_antechamber_preparation.yaml --campaign "$CAMPAIGN"
    # Partial failure OK
    echo ""

    echo "[01b] Grid Generation"
    python 02_scripts/01b_grid_generation.py --config 03_configs/01b_grid_generation.yaml --campaign "$CAMPAIGN"
    if [ $? -ne 0 ]; then echo "FAILED at 01b (grids)"; exit 1; fi
    echo ""

    echo "[01c] DOCK6 Docking"
    python 02_scripts/01c_dock6_run.py --config 03_configs/01c_dock6_run.yaml --campaign "$CAMPAIGN"
    # Partial failure OK
    echo ""

    echo "[01d] DOCK6 Score Collection"
    python 02_scripts/01d_score_collection.py --config 03_configs/01d_score_collection.yaml --campaign "$CAMPAIGN"
    if [ $? -ne 0 ]; then echo "FAILED at 01d (scores)"; exit 1; fi
    echo ""
fi

# =============================================================================
# 02 — GNINA ENGINE
# =============================================================================

if [ "$ENGINE" = "gnina" ] || [ "$ENGINE" = "both" ]; then
    echo "========== 02: GNINA ENGINE (Vina + CNN) =========="
    echo ""

    echo "[02a] GNINA Preparation (binding box + inputs JSON)"
    python 02_scripts/02a_gnina_preparation.py --config 03_configs/02a_gnina_preparation.yaml --campaign "$CAMPAIGN"
    if [ $? -ne 0 ]; then echo "FAILED at 02a (gnina prep)"; exit 1; fi
    echo ""

    echo "[02b] GNINA Docking (Vina + CNN scoring)"
    python 02_scripts/02b_gnina_runner.py --config 03_configs/02b_gnina_runner.yaml --campaign "$CAMPAIGN"
    # Partial failure OK
    echo ""

    echo "[02c] GNINA Score Collection"
    python 02_scripts/02c_gnina_score_collector.py --config 03_configs/02c_gnina_score_collector.yaml --campaign "$CAMPAIGN"
    if [ $? -ne 0 ]; then echo "FAILED at 02c (gnina scores)"; exit 1; fi
    echo ""
fi

# =============================================================================
# DONE
# =============================================================================

echo "============================================================"
echo "  Pipeline complete: ${CAMPAIGN_ID} (${ENGINE})"
echo "  Finished: $(date '+%Y-%m-%d %H:%M')"
echo "============================================================"
echo ""
echo "Results in: 05_results/${CAMPAIGN_ID}/"
if [ "$ENGINE" = "dock6" ] || [ "$ENGINE" = "both" ]; then
    echo "  DOCK6: 01d_score_collection/dock6_scores.xlsx"
fi
if [ "$ENGINE" = "gnina" ] || [ "$ENGINE" = "both" ]; then
    echo "  GNINA: 02c_gnina_scores/gnina_scores.xlsx"
fi
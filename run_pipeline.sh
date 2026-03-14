#!/bin/bash
# =============================================================================
# run_pipeline.sh — Run the complete docking pipeline
# =============================================================================
#
# Usage:
#   bash run_pipeline.sh <campaign_id>
#
# Example:
#   bash run_pipeline.sh example_campaign
#
# Requirements:
#   - conda activate molecular_docking_env
#   - DOCK6 and ChimeraX installed
#   - Campaign configured in 04_data/campaigns/{campaign_id}/
#
# =============================================================================

set -e

CAMPAIGN_ID="${1:-}"
if [ -z "$CAMPAIGN_ID" ]; then
    echo "Usage: bash run_pipeline.sh <campaign_id>"
    echo ""
    echo "Available campaigns:"
    ls -1 04_data/campaigns/ 2>/dev/null || echo "  (none)"
    exit 1
fi

CAMPAIGN="04_data/campaigns/${CAMPAIGN_ID}/campaign_config.yaml"

if [ ! -f "$CAMPAIGN" ]; then
    echo "ERROR: ${CAMPAIGN} not found"
    echo ""
    echo "Create a campaign first:"
    echo "  cp -r 04_data/campaigns/example_campaign 04_data/campaigns/${CAMPAIGN_ID}"
    exit 1
fi

echo "============================================================"
echo "  molecular_docking — Full pipeline"
echo "  Campaign: ${CAMPAIGN_ID}"
echo "  Started:  $(date '+%Y-%m-%d %H:%M')"
echo "============================================================"
echo ""

echo "[1/5] 00a — Molecule Parser"
python 02_scripts/00a_molecule_parser.py --config 03_configs/00a_molecule_parser.yaml --campaign "$CAMPAIGN"
echo ""

echo "[2/5] 00b — Receptor Preparation"
python 02_scripts/00b_receptor_preparation.py --config 03_configs/00b_receptor_preparation.yaml --campaign "$CAMPAIGN"
echo ""

echo "[3/5] 00c-00d — Ligand Preparation"
python 02_scripts/00c_ionization_profiling.py --config 03_configs/00c_ionization_profiling.yaml --campaign "$CAMPAIGN"
python 02_scripts/00d_antechamber_preparation.py --config 03_configs/00d_antechamber_preparation.yaml --campaign "$CAMPAIGN"
echo ""

echo "[4/5] 01a — Grid Generation"
python 02_scripts/01a_grid_generation.py --config 03_configs/01a_grid_generation.yaml --campaign "$CAMPAIGN"
echo ""

echo "[5/5] 01b — DOCK6 Docking"
python 02_scripts/01b_dock6_run.py --config 03_configs/01b_dock6_run.yaml --campaign "$CAMPAIGN"
echo ""

echo "============================================================"
echo "  Done — $(date '+%Y-%m-%d %H:%M')"
echo "  Results: 05_results/${CAMPAIGN_ID}/"
echo "============================================================"

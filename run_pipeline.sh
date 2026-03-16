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
# Critical modules (00a, 00b, 00c, 00e, 01a) stop the pipeline on failure.
# Batch modules (00d, 01b) continue even if some molecules fail.
#
# =============================================================================

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

# --- Critical steps (stop on failure) ---

echo "[1/7] 00a — Molecule Parser"
python 02_scripts/00a_molecule_parser.py --config 03_configs/00a_molecule_parser.yaml --campaign "$CAMPAIGN"
if [ $? -ne 0 ]; then echo "FAILED at 00a"; exit 1; fi
echo ""

echo "[2/7] 00b — Receptor Preparation"
python 02_scripts/00b_receptor_preparation.py --config 03_configs/00b_receptor_preparation.yaml --campaign "$CAMPAIGN"
if [ $? -ne 0 ]; then echo "FAILED at 00b"; exit 1; fi
echo ""

echo "[3/7] 00c-00d — Ligand Preparation"
python 02_scripts/00c_ionization_profiling.py --config 03_configs/00c_ionization_profiling.yaml --campaign "$CAMPAIGN"
if [ $? -ne 0 ]; then echo "FAILED at 00c"; exit 1; fi

# 00d: partial failure is OK (some molecules may fail antechamber)
python 02_scripts/00d_antechamber_preparation.py --config 03_configs/00d_antechamber_preparation.yaml --campaign "$CAMPAIGN"
echo ""

echo "[4/7] 00e — Binding Site Definition"
python 02_scripts/00e_binding_site_definition.py --config 03_configs/00e_binding_site_definition.yaml --campaign "$CAMPAIGN"
if [ $? -ne 0 ]; then echo "FAILED at 00e"; exit 1; fi
echo ""

echo "[5/7] 01a — Grid Generation"
python 02_scripts/01a_grid_generation.py --config 03_configs/01a_grid_generation.yaml --campaign "$CAMPAIGN"
if [ $? -ne 0 ]; then echo "FAILED at 01a"; exit 1; fi
echo ""

# 01b: partial failure is OK (some molecules may fail docking)
echo "[6/7] 01b — DOCK6 Docking"
python 02_scripts/01b_dock6_run.py --config 03_configs/01b_dock6_run.yaml --campaign "$CAMPAIGN"
echo ""

echo "[7/7] 02a — Score Collection"
python 02_scripts/02a_score_collection.py --config 03_configs/02a_score_collection.yaml --campaign "$CAMPAIGN"
if [ $? -ne 0 ]; then echo "FAILED at 02a"; exit 1; fi
echo ""

echo "============================================================"
echo "  Done — $(date '+%Y-%m-%d %H:%M')"
echo "  Results: 05_results/${CAMPAIGN_ID}/"
echo "============================================================"
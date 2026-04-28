#!/bin/bash
# =============================================================================
# run_pipeline.sh — Run the virtual screening pipeline
# =============================================================================
#
# Architecture (v7.0 — 2026-03-25):
#
#   Screening-only pipeline:
#     00a → 00b → 00c → 00d → 01a-01e → 02a-02c → 04a-04e → 05a-05g
#
#   Engines run first (01x, 02x), then all analysis (04x, 05x).
#
# Usage:
#   bash run_pipeline.sh <campaign_id> [engine] [start] [stop]
#
# Examples:
#   bash run_pipeline.sh XT1_hits                    # everything, both engines
#   bash run_pipeline.sh XT1_hits dock6              # DOCK6 only
#   bash run_pipeline.sh XT1_hits gnina              # GNINA only
#   bash run_pipeline.sh XT1_hits both 01c 01e       # both engines, range
#   bash run_pipeline.sh XT1_hits dock6 01c 01e      # DOCK6, range
#
# For nohup (background server runs):
#   # Full pipeline, both engines
#   nohup bash run_pipeline.sh XT1_hits both > logs/xt1_both.log 2>&1 &
#   nohup bash run_pipeline.sh SD1_pharmit_pH63_SD1_1 both > logs/SD1_pharmit_pH63_SD1_full.log 2>&1 &
#
#   # Single engine, full run
#   nohup bash run_pipeline.sh XT1_hits dock6 > logs/xt1_dock6.log 2>&1 &
#   nohup bash run_pipeline.sh XT1_hits gnina > logs/xt1_gnina.log 2>&1 &
#
#   # Resume after interruption (e.g. 01a done, continue from 01b)
#   nohup bash run_pipeline.sh XT1_hits dock6 01b 04e > logs/xt1_dock6_from01b.log 2>&1 &
#
#   # Analysis only (docking already finished)
#   nohup bash run_pipeline.sh XT1_hits dock6 04a 04e > logs/xt1_dock6_analysis.log 2>&1 &
#   nohup bash run_pipeline.sh XT1_hits gnina 05a 05g > logs/xt1_gnina_analysis.log 2>&1 &
#
#   # Monitor running jobs
#   jobs -l                    # list background jobs
#   tail -f logs/xt1_dock6.log # follow progress
#
# =============================================================================

set -e

CAMPAIGN_ID="${1:-}"
ENGINE="${2:-both}"
START="${3:-00a}"
STOP="${4:-07a}"

if [ -z "$CAMPAIGN_ID" ]; then
    echo "Usage:"
    echo "  bash run_pipeline.sh <campaign_id> [engine] [start] [stop]"
    echo ""
    echo "  engine = dock6 | gnina | both (default: both)"
    echo "  start/stop = module codes (default: 00a to 07a)"
    echo ""
    echo "Pipeline sequence:"
    echo "  00a → 00b → 00c → 00d → 01a-01e → 02a-02c → 04a-04e → 05a-05g → 07a"
    echo ""
    echo "Modules:"
    echo "  00a  Molecule Parser           00b  Receptor Preparation"
    echo "  00c  Ionization Profiling      00d  Binding Site Definition"
    echo "  01a  Antechamber (DOCK6)       01b  Grid Generation"
    echo "  01c  DOCK6 Docking             01d  Footprint Re-scoring"
    echo "  01e  DOCK6 Score Collection"
    echo "  02a  GNINA Preparation         02b  GNINA Docking"
    echo "  02c  GNINA Scores"
    echo "  04a  Score Ranking (DOCK6)     04b  Footprint Analysis"
    echo "  04c  Binding Modes             04d  Contact Mapping"
    echo "  04e  DOCK6 Campaign Report"
    echo "  05a  Parse & Fragment          05b  Fragment Clustering"
    echo "  05c  Score Decomposition       05d  Binding Site Hotspots"
    echo "  05e  Structure Export          05f  Contact Mapping"
    echo "  05g  Campaign Report"
    echo "  07a  Cross-Engine Comparison   (engine=both only)"
    echo ""
    echo "Examples:"
    echo "  bash run_pipeline.sh XT1_hits                    # everything"
    echo "  bash run_pipeline.sh XT1_hits dock6              # DOCK6 only"
    echo "  bash run_pipeline.sh XT1_hits both 01c 01e       # range"
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

# Receptor path (for 05d, 05e, 05f)
RECEPTOR="04_data/campaigns/${CAMPAIGN_ID}/receptor/XT1_6EJ7.pdb"

# Module ordering
declare -A MODULE_ORDER=(
    [00a]=1  [00b]=2  [00c]=3  [00d]=4
    [01a]=10 [01b]=11 [01c]=12 [01d]=13 [01e]=14
    [02a]=20 [02b]=21 [02c]=22
    [04a]=40 [04b]=41 [04c]=42 [04d]=43 [04e]=44
    [05a]=50 [05b]=51 [05c]=52 [05d]=53 [05e]=54 [05f]=55 [05g]=56
    [07a]=70
)

START_NUM=${MODULE_ORDER[$START]:-1}
STOP_NUM=${MODULE_ORDER[$STOP]:-70}

should_run() {
    local module=$1
    local num=${MODULE_ORDER[$module]:-0}
    [ $num -ge $START_NUM ] && [ $num -le $STOP_NUM ]
}

echo "============================================================"
echo "  molecular_docking — Pipeline v7.0"
echo "  Campaign: ${CAMPAIGN_ID}"
echo "  Engine:   ${ENGINE}"
echo "  Range:    ${START} → ${STOP}"
echo "  Started:  $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo ""

# =============================================================================
# 00 — SHARED PREPARATION
# =============================================================================

if should_run "00a"; then
    echo "[00a] Molecule Parser — $(date '+%H:%M:%S')"
    python 02_scripts/00a_molecule_parser.py --config 03_configs/00a_molecule_parser.yaml --campaign "$CAMPAIGN"
    echo ""
fi

if should_run "00b"; then
    echo "[00b] Receptor Preparation — $(date '+%H:%M:%S')"
    python 02_scripts/00b_receptor_preparation.py --config 03_configs/00b_receptor_preparation.yaml --campaign "$CAMPAIGN"
    echo ""
fi

if should_run "00c"; then
    echo "[00c] Ionization Profiling — $(date '+%H:%M:%S')"
    python 02_scripts/00c_ionization_profiling.py --config 03_configs/00c_ionization_profiling.yaml --campaign "$CAMPAIGN"
    echo ""
fi

if should_run "00d"; then
    echo "[00d] Binding Site Definition — $(date '+%H:%M:%S')"
    python 02_scripts/00d_binding_site_definition.py --config 03_configs/00d_binding_site_definition.yaml --campaign "$CAMPAIGN"
    echo ""
fi

# =============================================================================
# 01 — DOCK6 ENGINE
# =============================================================================

if [ "$ENGINE" = "dock6" ] || [ "$ENGINE" = "both" ]; then

    if should_run "01a"; then
        echo "[01a] Antechamber Preparation — $(date '+%H:%M:%S')"
        python 02_scripts/01a_antechamber_preparation.py --config 03_configs/01a_antechamber_preparation.yaml --campaign "$CAMPAIGN"
        echo ""
    fi

    if should_run "01b"; then
        echo "[01b] Grid Generation — $(date '+%H:%M:%S')"
        python 02_scripts/01b_grid_generation.py --config 03_configs/01b_grid_generation.yaml --campaign "$CAMPAIGN"
        echo ""
    fi

    if should_run "01c"; then
        echo "[01c] DOCK6 Docking — $(date '+%H:%M:%S')"
        python 02_scripts/01c_dock6_run.py --config 03_configs/01c_dock6_run.yaml --campaign "$CAMPAIGN"
        echo ""
    fi

    if should_run "01d"; then
        echo "[01d] Footprint Re-scoring — $(date '+%H:%M:%S')"
        python 02_scripts/01d_footprint_rescore.py --config 03_configs/01d_footprint_rescore.yaml --campaign "$CAMPAIGN"
        echo ""
    fi

    if should_run "01e"; then
        echo "[01e] DOCK6 Score Collection — $(date '+%H:%M:%S')"
        python 02_scripts/01e_score_collection.py --config 03_configs/01e_score_collection.yaml --campaign "$CAMPAIGN"
        echo ""
    fi
fi

# =============================================================================
# 02 — GNINA ENGINE
# =============================================================================

if [ "$ENGINE" = "gnina" ] || [ "$ENGINE" = "both" ]; then

    if should_run "02a"; then
        echo "[02a] GNINA Preparation — $(date '+%H:%M:%S')"
        python 02_scripts/02a_gnina_preparation.py --config 03_configs/02a_gnina_preparation.yaml --campaign "$CAMPAIGN"
        echo ""
    fi

    if should_run "02b"; then
        echo "[02b] GNINA Docking — $(date '+%H:%M:%S')"
        python 02_scripts/02b_gnina_runner.py --config 03_configs/02b_gnina_runner.yaml --campaign "$CAMPAIGN"
        echo ""
    fi

    if should_run "02c"; then
        echo "[02c] GNINA Score Collection — $(date '+%H:%M:%S')"
        python 02_scripts/02c_gnina_score_collector.py --config 03_configs/02c_gnina_score_collector.yaml --campaign "$CAMPAIGN"
        echo ""
    fi
fi

# =============================================================================
# 04 — DOCK6 ANALYSIS
# =============================================================================

if [ "$ENGINE" = "dock6" ] || [ "$ENGINE" = "both" ]; then

    if should_run "04a"; then
        echo "[04a] DOCK6 Score Ranking — $(date '+%H:%M:%S')"
        python 02_scripts/04a_score_ranking.py --config 03_configs/04a_score_ranking.yaml --campaign "$CAMPAIGN"
        echo ""
    fi

    if should_run "04b"; then
        echo "[04b] Footprint Analysis — $(date '+%H:%M:%S')"
        python 02_scripts/04b_footprint_analysis.py --config 03_configs/04b_footprint_analysis.yaml --campaign "$CAMPAIGN"
        echo ""
    fi

    if should_run "04c"; then
        echo "[04c] Binding Modes — $(date '+%H:%M:%S')"
        python 02_scripts/04c_binding_modes.py --config 03_configs/04c_binding_modes.yaml --campaign "$CAMPAIGN"
        echo ""
    fi

    if should_run "04d"; then
        echo "[04d] Contact Mapping — $(date '+%H:%M:%S')"
        python 02_scripts/04d_contact_mapping.py --config 03_configs/04d_contact_mapping.yaml --campaign "$CAMPAIGN"
        echo ""
    fi

    if should_run "04e"; then
        echo "[04e] DOCK6 Campaign Report — $(date '+%H:%M:%S')"
        python 02_scripts/04e_campaign_report.py --config 03_configs/04e_campaign_report.yaml --campaign "$CAMPAIGN"
        echo ""
    fi
fi

# =============================================================================
# 05 — GNINA ANALYSIS
# =============================================================================

# Determine engine for analysis
if [ "$ENGINE" = "both" ]; then
    ANALYSIS_ENGINE="gnina"
else
    ANALYSIS_ENGINE="$ENGINE"
fi

echo "========== m05: ${ANALYSIS_ENGINE^^} ANALYSIS =========="
echo ""

if should_run "05a"; then
    echo "[05a] Parse & Fragment — $(date '+%H:%M:%S')"
    python 02_scripts/05a_parse_and_fragment.py -c 03_configs/05a_parse_and_fragment.yaml --campaign "$CAMPAIGN" --engine "$ANALYSIS_ENGINE"
    echo ""
fi

if should_run "05b"; then
    echo "[05b] Fragment Clustering — $(date '+%H:%M:%S')"
    python 02_scripts/05b_fragment_clustering.py -c 03_configs/05b_fragment_clustering.yaml --campaign "$CAMPAIGN"
    echo ""
fi

if should_run "05c"; then
    echo "[05c] Score Decomposition — $(date '+%H:%M:%S')"
    python 02_scripts/05c_score_decomposition.py -c 03_configs/05c_score_decomposition.yaml --campaign "$CAMPAIGN"
    echo ""
fi

if should_run "05d"; then
    echo "[05d] Binding Site Hotspots — $(date '+%H:%M:%S')"
    python 02_scripts/05d_binding_site_hotspots.py -c 03_configs/05d_binding_site_hotspots.yaml --campaign "$CAMPAIGN"
    echo ""
fi

if should_run "05e"; then
    echo "[05e] Structure Export — $(date '+%H:%M:%S')"
    python 02_scripts/05e_structure_export.py -c 03_configs/05e_structure_export.yaml --campaign "$CAMPAIGN"
    echo ""
fi

if should_run "05f"; then
    echo "[05f] Contact Mapping — $(date '+%H:%M:%S')"
    python 02_scripts/05f_contact_mapping.py -c 03_configs/05f_contact_mapping.yaml --campaign "$CAMPAIGN" --receptor "$RECEPTOR"
    echo ""
fi

if should_run "05g"; then
    echo "[05g] Campaign Report — $(date '+%H:%M:%S')"
    python 02_scripts/05g_campaign_report.py -c 03_configs/05g_campaign_report.yaml --campaign "$CAMPAIGN"
    echo ""
fi

# =============================================================================
# 07 — CROSS-ENGINE COMPARISON
# =============================================================================

if [ "$ENGINE" = "both" ]; then
    if should_run "07a"; then
        echo "[07a] Cross-Engine Comparison — $(date '+%H:%M:%S')"
        python 02_scripts/07a_cross_engine_comparison.py \
            -c 03_configs/07a_cross_engine_comparison.yaml \
            --campaign "$CAMPAIGN"
        echo ""
    fi
fi

# =============================================================================
# DONE
# =============================================================================

echo "============================================================"
echo "  Screening complete"
echo "  Campaign: ${CAMPAIGN_ID}"
echo "  Engine:   ${ENGINE}"
echo "  Finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
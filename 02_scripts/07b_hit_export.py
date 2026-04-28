#!/usr/bin/env python3
"""
07b Hit Export - CLI
=======================
Export selected molecules (from Excel) as protonated SDF files
for downstream hit_validation.

Usage:
    python 02_scripts/07b_hit_export.py \
        --config 03_configs/07b_hit_export.yaml \
        --campaign 04_data/campaigns/UDX_pharmit_UDX_1_pH63/campaign_config.yaml

    # Override selection:
    python 02_scripts/07b_hit_export.py \
        --config 03_configs/07b_hit_export.yaml \
        --campaign 04_data/campaigns/UDX_pharmit_UDX_1_pH63/campaign_config.yaml \
        --selection selections/pharmacophore_v1.xlsx

Project: molecular_docking
Module: 07b
Version: 1.0
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from molecular_docking.m07_cross_engine.hit_export import run_hit_export

logger = logging.getLogger(__name__)


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def setup_log_file(log_path: Path, log_level: str = "INFO"):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(str(log_path), encoding="utf-8")
    fh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    fh.setFormatter(logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s"))
    logging.getLogger().addHandler(fh)


def main():
    parser = argparse.ArgumentParser(
        description="07b Hit Export — export selected molecules as SDF for hit_validation",
    )
    parser.add_argument("--config", "-c", type=str, required=True,
                        help="Module config YAML")
    parser.add_argument("--campaign", type=str, required=True,
                        help="Campaign config YAML")
    parser.add_argument("--selection", type=str, default=None,
                        help="Override selection Excel (relative to campaign dir or absolute)")
    parser.add_argument("--output", "-o", type=str, default=None,
                        help="Override output dir")
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    args = parser.parse_args()

    # =========================================================================
    # RESOLVE PARAMETERS
    # =========================================================================

    cc = load_yaml(args.campaign)
    campaign_dir = Path(args.campaign).parent
    campaign_id = cc.get("campaign_id", campaign_dir.name)

    # Module config
    mc = load_yaml(args.config)
    params = mc.get("parameters", {})

    log_level = args.log_level or params.get("log_level", "INFO")

    # --- Selection file ---
    selection_file = args.selection or params.get("selection_file")
    if not selection_file:
        logger.error("No selection file specified. Use --selection or set parameters.selection_file in config YAML.")
        return 1

    selection_path = Path(selection_file)
    if not selection_path.is_absolute():
        selection_path = campaign_dir / selection_path

    if not selection_path.exists():
        logger.error(f"Selection Excel not found: {selection_path}")
        return 1

    # --- gnina_scores.csv (auto-discover) ---
    results_base = Path("05_results") / campaign_id
    gnina_scores_csv = results_base / "02c_gnina_scores" / "gnina_scores.csv"

    if not gnina_scores_csv.exists():
        logger.error(f"gnina_scores.csv not found: {gnina_scores_csv}")
        logger.error("Run module 02c (gnina_scores) first.")
        return 1

    # --- Output directory ---
    selection_stem = selection_path.stem  # "druglikeness_v3" from "druglikeness_v3.xlsx"
    output_dir = args.output or str(results_base / "07b_hit_export" / selection_stem)

    # =========================================================================
    # SETUP LOGGING
    # =========================================================================

    if log_level:
        logging.getLogger().setLevel(getattr(logging, log_level.upper(), logging.INFO))

    out_p = Path(output_dir)
    out_p.mkdir(parents=True, exist_ok=True)
    setup_log_file(out_p / "07b_hit_export.log", log_level)

    # =========================================================================
    # EXECUTE
    # =========================================================================

    logger.info("=" * 60)
    logger.info("  MOLECULAR_DOCKING - Module 07b: Hit Export")
    logger.info("=" * 60)
    logger.info(f"Campaign:       {campaign_id}")
    logger.info(f"Selection:      {selection_path}")
    logger.info(f"GNINA scores:   {gnina_scores_csv}")
    logger.info(f"Output:         {output_dir}")

    result = run_hit_export(
        selection_xlsx=str(selection_path),
        gnina_scores_csv=str(gnina_scores_csv),
        output_dir=output_dir,
    )

    logger.info("")
    logger.info("=" * 60)
    logger.info(f"  {result['n_exported']}/{result['n_selected']} molecules exported")
    if result["n_missing"] > 0:
        logger.warning(f"  Missing: {result['missing']}")
    logger.info(f"  Summary: {result['summary_csv']}")
    logger.info(f"  SDFs:    {result['sdf_dir']}")
    logger.info("=" * 60)

    return 0 if result["success"] else 1


if __name__ == "__main__":
    sys.exit(main())

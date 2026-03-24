#!/usr/bin/env python3
"""
04a DOCK6 Score Ranking - CLI
================================
Rank molecules by Grid_Score and decompose into vdW + ES.

Input:  05_results/{campaign}/01e_score_collection/dock6_all_poses.csv
Output: 05_results/{campaign}/04_dock6_analysis/04a_score_ranking/

Project: molecular_docking
Module: 04a (DOCK6 analysis)
Version: 1.1 (2026-03-23) — path fix: 05_dock6 → 04_dock6_analysis
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from molecular_docking.m04_dock6_analysis.score_ranking import run_score_ranking

logger = logging.getLogger(__name__)


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def setup_log_file(log_path, log_level="INFO"):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(str(log_path), mode="w", encoding="utf-8")
    fh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    fh.setFormatter(logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s"))
    logging.getLogger().addHandler(fh)


def main():
    parser = argparse.ArgumentParser(
        description="04a DOCK6 Score Ranking — Grid_Score ranking + vdW/ES decomposition",
    )
    parser.add_argument("--config", "-c", type=str, help="Module config YAML")
    parser.add_argument("--campaign", type=str, help="Campaign config YAML")
    parser.add_argument("--all-poses", type=str, default=None,
                        help="Direct path to dock6_all_poses.csv")
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--top-n", type=int, default=None)
    parser.add_argument("--log-level", type=str, default=None)

    args = parser.parse_args()

    # Defaults
    all_poses_csv = None
    output_dir = None
    campaign_id = "direct"
    score_key = "Grid_Score"
    top_n = 0
    log_level = "INFO"

    # --- Campaign config ---
    if args.campaign:
        cc = load_yaml(args.campaign)
        campaign_id = cc.get("campaign_id", Path(args.campaign).parent.name)
        all_poses_csv = str(
            Path("05_results") / campaign_id / "01e_score_collection" / "dock6_all_poses.csv"
        )
        output_dir = str(
            Path("05_results") / campaign_id / "04_dock6_analysis" / "04a_score_ranking"
        )

    # --- Module config ---
    if args.config:
        mc = load_yaml(args.config)
        params = mc.get("parameters", {})
        score_key = params.get("score_key", score_key)
        top_n = params.get("top_n", top_n)
        log_level = params.get("log_level", log_level)

    # --- CLI overrides ---
    if args.all_poses:
        all_poses_csv = args.all_poses
    if args.output:
        output_dir = args.output
    if args.top_n is not None:
        top_n = args.top_n
    if args.log_level:
        log_level = args.log_level

    if not all_poses_csv or not Path(all_poses_csv).exists():
        logger.error(f"All poses CSV not found: {all_poses_csv}")
        logger.error("Run module 01e first.")
        return 1

    if not output_dir:
        output_dir = "05_results/04_dock6_analysis/04a_score_ranking"

    logging.getLogger().setLevel(getattr(logging, log_level.upper(), logging.INFO))
    setup_log_file(Path(output_dir) / "04a_score_ranking.log", log_level)

    logger.info("=" * 60)
    logger.info("  MOLECULAR_DOCKING — Module 04a: DOCK6 Score Ranking")
    logger.info("=" * 60)
    logger.info(f"  Campaign: {campaign_id}")

    result = run_score_ranking(
        all_poses_csv=all_poses_csv,
        output_dir=output_dir,
        score_key=score_key,
        top_n=top_n,
    )

    if not result.get("success"):
        logger.error(f"Error: {result.get('error')}")
        return 1

    logger.info(f"\nNext: python 02_scripts/04b_footprint_analysis.py "
                f"--config 03_configs/04b_footprint_analysis.yaml "
                f"--campaign {args.campaign or '<campaign_config.yaml>'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

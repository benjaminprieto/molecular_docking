#!/usr/bin/env python3
"""
07a Cross-Engine Comparison - CLI
=====================================
Compare GNINA and DOCK6 results to identify consensus hits.

Usage:
    python 02_scripts/07a_cross_engine_comparison.py \
        --config 03_configs/07a_cross_engine_comparison.yaml \
        --campaign 04_data/campaigns/UDX_namiki_pH63/campaign_config.yaml

Project: molecular_docking
Module: 07a
Version: 1.0
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)-8s | %(message)s')
sys.path.insert(0, str(Path(__file__).parent.parent / '01_src'))

from molecular_docking.m07_cross_engine.cross_engine_comparison import run_cross_engine_comparison

logger = logging.getLogger(__name__)


def load_yaml(path):
    with open(path) as f:
        return yaml.safe_load(f)


def setup_log_file(log_path: Path, log_level: str = "INFO"):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(str(log_path), encoding="utf-8")
    fh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    fh.setFormatter(logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s"))
    logging.getLogger().addHandler(fh)


def main():
    parser = argparse.ArgumentParser(description='07a Cross-Engine Comparison')
    parser.add_argument('--config', '-c', type=str, help='Module YAML config')
    parser.add_argument('--campaign', type=str, required=True, help='Campaign config YAML')
    parser.add_argument('--output', '-o', type=str, default=None)
    parser.add_argument('--top-n', type=int, default=None, help='Top N for consensus (default: 20)')
    parser.add_argument('--log-level', type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = parser.parse_args()

    cc = load_yaml(args.campaign)
    campaign_dir = Path(args.campaign).parent
    campaign_id = cc.get('campaign_id', campaign_dir.name)

    # Module config
    params = {}
    if args.config:
        mc = load_yaml(args.config)
        params = mc.get('parameters', {})

    top_n = args.top_n or params.get('top_n', 20)
    log_level = args.log_level or params.get('log_level', 'INFO')

    # Output directory
    results_base = Path('05_results') / campaign_id
    output_dir = args.output or str(results_base / '07_cross_engine')

    out_p = Path(output_dir)
    out_p.mkdir(parents=True, exist_ok=True)
    setup_log_file(out_p / "07a_cross_engine.log", log_level)

    # Auto-detect input files
    gnina_scores_csv = params.get('gnina_scores_csv')
    if not gnina_scores_csv:
        auto = results_base / '02c_gnina_scores' / 'gnina_scores.csv'
        if auto.exists():
            gnina_scores_csv = str(auto)

    dock6_scores_csv = params.get('dock6_scores_csv')
    if not dock6_scores_csv:
        auto = results_base / '01e_score_collection' / 'dock6_scores.csv'
        if auto.exists():
            dock6_scores_csv = str(auto)

    gnina_composite_csv = params.get('gnina_composite_csv')
    if not gnina_composite_csv:
        auto = results_base / '05_gnina_analysis' / '05g_campaign_report' / 'composite_ranking.csv'
        if auto.exists():
            gnina_composite_csv = str(auto)

    dock6_all_poses_csv = params.get('dock6_all_poses_csv')
    if not dock6_all_poses_csv:
        auto = results_base / '01e_score_collection' / 'dock6_all_poses.csv'
        if auto.exists():
            dock6_all_poses_csv = str(auto)

    dock6_composite_csv = params.get('dock6_composite_csv')
    if not dock6_composite_csv:
        auto = results_base / '04_dock6_analysis' / '04e_campaign_report' / 'composite_ranking.csv'
        if auto.exists():
            dock6_composite_csv = str(auto)

    # Validate mandatory inputs
    if not gnina_scores_csv or not Path(gnina_scores_csv).exists():
        logger.error(f"GNINA scores not found. Expected: {results_base / '02c_gnina_scores' / 'gnina_scores.csv'}")
        logger.error("Run GNINA pipeline first (02a-02c)")
        sys.exit(1)

    if not dock6_scores_csv or not Path(dock6_scores_csv).exists():
        logger.error(f"DOCK6 scores not found. Expected: {results_base / '01e_score_collection' / 'dock6_scores.csv'}")
        logger.error("Run DOCK6 pipeline first (01a-01e)")
        sys.exit(1)

    logger.info("=" * 60)
    logger.info("M07a: CROSS-ENGINE COMPARISON v1.0")
    logger.info("=" * 60)
    logger.info(f"Campaign:         {campaign_id}")
    logger.info(f"GNINA scores:     {gnina_scores_csv}")
    logger.info(f"GNINA composite:  {gnina_composite_csv or 'not found'}")
    logger.info(f"DOCK6 scores:     {dock6_scores_csv}")
    logger.info(f"DOCK6 all poses:  {dock6_all_poses_csv or 'not found'}")
    logger.info(f"DOCK6 composite:  {dock6_composite_csv or 'not found'}")
    logger.info(f"Top N:            {top_n}")
    logger.info(f"Output:           {output_dir}")

    result = run_cross_engine_comparison(
        gnina_scores_csv=gnina_scores_csv,
        dock6_scores_csv=dock6_scores_csv,
        output_dir=output_dir,
        campaign_id=campaign_id,
        gnina_composite_csv=gnina_composite_csv,
        dock6_all_poses_csv=dock6_all_poses_csv,
        dock6_composite_csv=dock6_composite_csv,
        top_n=top_n,
    )

    if not result.get('success'):
        logger.error(f"Failed: {result.get('error', 'unknown')}")
        sys.exit(1)

    logger.info(f"\nOpen in browser: {result['report_path']}")


if __name__ == '__main__':
    main()

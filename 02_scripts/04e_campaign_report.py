#!/usr/bin/env python3
"""
04e DOCK6 Campaign Report - CLI
===================================
Generates consolidated HTML report for DOCK6 analysis with composite ranking.

Usage:
    python 02_scripts/04e_campaign_report.py -c 03_configs/04e_campaign_report.yaml \
        --campaign 04_data/campaigns/UDX_namiki_pH63/campaign_config.yaml

Project: molecular_docking
Module: 04e (DOCK6 analysis)
Version: 3.0 — multiplicative composite (Grid_Score × weighted_consistency)
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)-8s | %(message)s')
sys.path.insert(0, str(Path(__file__).parent.parent / '01_src'))

from molecular_docking.m04_dock6_analysis.campaign_report import run_campaign_report

logger = logging.getLogger(__name__)


def load_yaml(path):
    with open(path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


def setup_log_file(log_path: Path, log_level: str = "INFO"):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(str(log_path), encoding="utf-8")
    fh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    fh.setFormatter(logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s"))
    logging.getLogger().addHandler(fh)


def main():
    parser = argparse.ArgumentParser(description='04e DOCK6 Campaign Report')
    parser.add_argument('--config', '-c', type=str, help='Module YAML config')
    parser.add_argument('--campaign', type=str, required=True, help='Campaign config YAML')
    parser.add_argument('--output', '-o', type=str, default=None)
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

    # Results directory
    results_base = Path('05_results') / campaign_id / '04_dock6_analysis'
    output_subdir = params.get('output_subdir', '04e_campaign_report')
    output_dir = args.output or str(results_base / output_subdir)

    log_level = args.log_level or params.get("log_level", "INFO")
    out_p = Path(output_dir)
    out_p.mkdir(parents=True, exist_ok=True)
    setup_log_file(out_p / "04e_campaign_report.log", log_level)

    # Drug-likeness source (try gnina_scores.csv, then unique_molecules.csv)
    gnina_scores_csv = params.get("gnina_scores_csv")
    if not gnina_scores_csv:
        auto_gnina = Path('05_results') / campaign_id / '02c_gnina_scores' / 'gnina_scores.csv'
        if auto_gnina.exists():
            gnina_scores_csv = str(auto_gnina)

    molecules_csv = None
    auto_mol = Path('05_results') / campaign_id / '00a_molecule_parser' / 'unique_molecules.csv'
    if auto_mol.exists():
        molecules_csv = str(auto_mol)

    # All poses CSV from 01e (for mean Grid_Score composite)
    all_poses_csv = None
    auto_poses = Path('05_results') / campaign_id / '01e_score_collection' / 'dock6_all_poses.csv'
    if auto_poses.exists():
        all_poses_csv = str(auto_poses)

    logger.info("=" * 60)
    logger.info("M04e: DOCK6 CAMPAIGN REPORT v4.0")
    logger.info("=" * 60)
    logger.info(f"Campaign:       {campaign_id}")
    logger.info(f"Results base:   {results_base}")
    logger.info(f"Output:         {output_dir}")
    logger.info(f"All poses:      {all_poses_csv or 'not found'}")
    logger.info(f"Drug-likeness:  {gnina_scores_csv or molecules_csv or 'none'}")

    result = run_campaign_report(
        results_base=str(results_base),
        output_dir=output_dir,
        campaign_id=campaign_id,
        gnina_scores_csv=gnina_scores_csv,
        molecules_csv=molecules_csv,
        all_poses_csv=all_poses_csv,
    )

    if not result.get('success'):
        logger.error(f"Failed: {result.get('error')}")
        return 1

    logger.info(f"\nOpen in browser: {result['report_path']}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
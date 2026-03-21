#!/usr/bin/env python3
"""
02c GNINA Score Collector - CLI (Module 02c)
=============================================
Collects, ranks, and exports GNINA scores to CSV + Excel.

Hardcoded upstream paths:
    - 05_results/{campaign_id}/02b_gnina_run/gnina_results.csv
    - 05_results/{campaign_id}/00a_molecule_parser/unique_molecules.csv
    - Output: 05_results/{campaign_id}/02c_gnina_scores/

Project: molecular_docking
Module: 02c (GNINA engine)
Version: 1.0
"""

import argparse
import yaml
import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)-8s | %(message)s')

sys.path.insert(0, str(Path(__file__).parent.parent / '01_src'))

from molecular_docking.m02_gnina.gnina_score_collector import run_gnina_score_collection

logger = logging.getLogger(__name__)


def load_yaml(path):
    with open(path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


def setup_log_file(log_path, log_level="INFO"):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(str(log_path), mode='w', encoding='utf-8')
    fh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    fh.setFormatter(logging.Formatter('%(asctime)s | %(levelname)-8s | %(name)s | %(message)s'))
    logging.getLogger().addHandler(fh)


def main():
    parser = argparse.ArgumentParser(
        description='02c GNINA Score Collector — rank & export GNINA scores (Module 02c)',
    )
    parser.add_argument('--config', '-c', type=str, help='Module YAML config')
    parser.add_argument('--campaign', type=str, help='Campaign config YAML')
    parser.add_argument('--results', type=str, default=None, help='GNINA results CSV (direct mode)')
    parser.add_argument('--molecules', type=str, default=None)
    parser.add_argument('--output', '-o', type=str, default=None)
    parser.add_argument('--top-n', type=int, default=None)
    parser.add_argument('--rank-by', type=str, default=None,
                        choices=['cnn_affinity', 'cnn_score', 'vina_affinity', 'consensus'])
    parser.add_argument('--no-excel', action='store_true')
    parser.add_argument('--log-level', type=str, default=None)
    parser.add_argument('--verbose', '-v', action='store_true')

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Defaults
    gnina_results_csv = None
    molecules_csv = None
    output_dir = None
    campaign_id = "direct"
    top_n = 10
    rank_by = "cnn_affinity"
    export_excel = True
    keep_all_poses = False
    name_column = "Name"
    log_level = "INFO"

    if args.campaign:
        cc = load_yaml(args.campaign)
        campaign_id = cc.get("campaign_id", Path(args.campaign).parent.name)
        gnina_results_csv = str(Path("05_results") / campaign_id / "02b_gnina_run" / "gnina_results.csv")
        molecules_csv = str(Path("05_results") / campaign_id / "00a_molecule_parser" / "unique_molecules.csv")
        output_dir = str(Path("05_results") / campaign_id / "02c_gnina_scores")
    elif args.results:
        gnina_results_csv = args.results
        output_dir = args.output or "05_results/02c_gnina_scores"
    else:
        parser.error("Requires --campaign or --results gnina_results.csv")

    if args.config:
        mc = load_yaml(args.config)
        params = mc.get("parameters", {})
        top_n = params.get("top_n", top_n)
        rank_by = params.get("rank_by", rank_by)
        export_excel = params.get("export_excel", export_excel)
        keep_all_poses = params.get("keep_all_poses", False)
        name_column = params.get("name_column", name_column)
        log_level = params.get("log_level", log_level)
        subdir = mc.get("outputs", {}).get("subdir", "02c_gnina_scores")
        if args.campaign and not args.output:
            output_dir = str(Path("05_results") / campaign_id / subdir)

    if args.output: output_dir = args.output
    if args.molecules: molecules_csv = args.molecules
    if args.top_n: top_n = args.top_n
    if args.rank_by: rank_by = args.rank_by
    if args.no_excel: export_excel = False
    if args.log_level: log_level = args.log_level

    if not gnina_results_csv or not Path(gnina_results_csv).exists():
        logger.error(f"GNINA results not found: {gnina_results_csv}. Run 02b first.")
        return 1

    log_path = Path(output_dir) / '02c_gnina_scores.log'
    setup_log_file(log_path, log_level)

    try:
        logger.info("=" * 60)
        logger.info("M02c: GNINA SCORE COLLECTION (v1.0)")
        logger.info("=" * 60)
        logger.info(f"Campaign:         {campaign_id}")
        logger.info(f"GNINA results:    {gnina_results_csv}")
        logger.info(f"Rank by:          {rank_by}")
        logger.info(f"Top N:            {top_n}")
        logger.info("-" * 60)

        result = run_gnina_score_collection(
            gnina_results_csv=gnina_results_csv,
            output_dir=output_dir,
            molecules_csv=molecules_csv,
            name_column=name_column,
            top_n=top_n,
            rank_by=rank_by,
            export_excel=export_excel,
            keep_all_poses=keep_all_poses,
        )

        if not result.get("success"):
            logger.error("Score collection failed")
            return 1

        logger.info("")
        logger.info("OUTPUT FILES")
        logger.info(f"  CSV:     {result['scores_csv']}")
        if result.get('scores_xlsx'):
            logger.info(f"  Excel:   {result['scores_xlsx']}")
        logger.info(f"  Summary: {result['summary_txt']}")
        logger.info("=" * 60)
        logger.info("M02c COMPLETE — GNINA pipeline finished")
        logger.info("=" * 60)
        return 0

    except Exception as e:
        logger.error(f"Failed: {e}")
        if args.verbose:
            import traceback; traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
#!/usr/bin/env python3
"""
02c Vina Score Collector - CLI (Module 02c)
=============================================
Collects, ranks, and exports Vina docking scores to CSV + Excel.

Reads campaign_config.yaml to determine:
    - campaign_id       → path construction

Hardcoded upstream paths:
    - 05_results/{campaign_id}/02b_vina_run/vina_results.csv
    - 05_results/{campaign_id}/00a_molecule_parser/unique_molecules.csv
    - Output: 05_results/{campaign_id}/02c_vina_scores/

Preconditions:
    - Module 02b must have been run (vina_results.csv)
    - Module 00a should have been run (optional metadata enrichment)

Usage:
    python 02_scripts/02c_vina_score_collector.py --config 03_configs/02c_vina_score_collector.yaml --campaign 04_data/campaigns/example_campaign/campaign_config.yaml

Project: molecular_docking
Module: 02c
Version: 1.0
"""

import argparse
import yaml
import logging
import sys
from pathlib import Path
from datetime import datetime

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)-8s | %(message)s'
)

sys.path.insert(0, str(Path(__file__).parent.parent / '01_src'))

from molecular_docking.m02_vina.vina_score_collector import run_vina_score_collection

logger = logging.getLogger(__name__)


def load_yaml(path: str) -> dict:
    with open(path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


def setup_log_file(log_path: Path, log_level: str = "INFO"):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(str(log_path), mode='w', encoding='utf-8')
    fh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    fh.setFormatter(logging.Formatter(
        '%(asctime)s | %(levelname)-8s | %(name)s | %(message)s'
    ))
    logging.getLogger().addHandler(fh)
    logger.info(f"Log file: {log_path}")


def main():
    parser = argparse.ArgumentParser(
        description='02c Vina Score Collector — rank & export Vina scores (Module 02c)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Standard campaign mode
  python 02_scripts/02c_vina_score_collector.py --config 03_configs/02c_vina_score_collector.yaml --campaign 04_data/campaigns/example_campaign/campaign_config.yaml

  # Direct mode
  python 02_scripts/02c_vina_score_collector.py --results vina_results.csv --output results/02c_vina_scores/

  # With metadata enrichment
  python 02_scripts/02c_vina_score_collector.py --results vina_results.csv --molecules unique_molecules.csv --output results/02c_vina_scores/

  # No Excel
  python 02_scripts/02c_vina_score_collector.py --config 03_configs/02c_vina_score_collector.yaml --campaign 04_data/campaigns/example_campaign/campaign_config.yaml --no-excel
        '''
    )

    parser.add_argument('--config', '-c', type=str, help='Module YAML config')
    parser.add_argument('--campaign', type=str, help='Campaign config YAML')

    # Direct mode
    parser.add_argument('--results', type=str, default=None,
                        help='Vina results CSV from 02b (direct mode)')
    parser.add_argument('--molecules', type=str, default=None,
                        help='Molecules CSV from 00a (optional)')

    # Overrides
    parser.add_argument('--output', '-o', type=str, default=None)
    parser.add_argument('--top-n', type=int, default=None)
    parser.add_argument('--no-excel', action='store_true', help='Disable Excel export')
    parser.add_argument('--log-level', type=str, default=None)
    parser.add_argument('--verbose', '-v', action='store_true')

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.log_level:
        logging.getLogger().setLevel(getattr(logging, args.log_level.upper(), logging.INFO))

    # =========================================================================
    # DEFAULTS
    # =========================================================================
    vina_results_csv = None
    molecules_csv = None
    output_dir = None
    campaign_id = "direct"
    top_n = 10
    export_excel = True
    name_column = "Name"
    log_level = "INFO"

    # =========================================================================
    # RESOLVE PARAMETERS
    # =========================================================================
    if args.campaign:
        cc = load_yaml(args.campaign)
        campaign_dir = Path(args.campaign).parent
        campaign_id = cc.get("campaign_id", campaign_dir.name)

        # Hardcoded upstream paths
        vina_results_csv = str(
            Path("05_results") / campaign_id / "02b_vina_run" / "vina_results.csv"
        )
        molecules_csv = str(
            Path("05_results") / campaign_id / "00a_molecule_parser" / "unique_molecules.csv"
        )
        output_dir = str(Path("05_results") / campaign_id / "02c_vina_scores")

    elif args.results:
        vina_results_csv = args.results
        output_dir = args.output or "05_results/02c_vina_scores"

    else:
        parser.error("Requires --campaign or --results vina_results.csv")

    # Module config
    if args.config:
        mc = load_yaml(args.config)
        params = mc.get("parameters", {})
        top_n = params.get("top_n", top_n)
        export_excel = params.get("export_excel", export_excel)
        name_column = params.get("name_column", name_column)
        log_level = params.get("log_level", log_level)

        subdir = mc.get("outputs", {}).get("subdir", "02c_vina_scores")
        if args.campaign and not args.output:
            output_dir = str(Path("05_results") / campaign_id / subdir)

    # CLI overrides
    if args.output:
        output_dir = args.output
    if args.molecules:
        molecules_csv = args.molecules
    if args.top_n is not None:
        top_n = args.top_n
    if args.no_excel:
        export_excel = False
    if args.log_level:
        log_level = args.log_level

    # Validate
    if not vina_results_csv or not Path(vina_results_csv).exists():
        logger.error(f"Vina results not found: {vina_results_csv}")
        logger.error("Run module 02b first.")
        return 1

    # =========================================================================
    # SETUP LOG FILE
    # =========================================================================
    log_path = Path(output_dir) / '02c_vina_scores.log'
    setup_log_file(log_path, log_level)

    # =========================================================================
    # RUN
    # =========================================================================
    try:
        logger.info("=" * 60)
        logger.info("M02c: VINA SCORE COLLECTION (v1.0)")
        logger.info("=" * 60)
        logger.info(f"Timestamp:        {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info(f"Campaign:         {campaign_id}")
        logger.info(f"Vina results:     {vina_results_csv}")
        logger.info(f"Molecules CSV:    {molecules_csv or 'None'}")
        logger.info(f"Output dir:       {output_dir}")
        logger.info(f"Top N:            {top_n}")
        logger.info(f"Export Excel:     {export_excel}")
        logger.info(f"Log file:         {log_path}")
        logger.info("-" * 60)

        result = run_vina_score_collection(
            vina_results_csv=vina_results_csv,
            output_dir=output_dir,
            molecules_csv=molecules_csv,
            name_column=name_column,
            top_n=top_n,
            export_excel=export_excel,
        )

        if not result.get("success"):
            logger.error(f"Score collection failed")
            return 1

        # Summary
        logger.info("")
        logger.info("-" * 60)
        logger.info("OUTPUT FILES")
        logger.info("-" * 60)
        logger.info(f"  CSV:     {result['scores_csv']}")
        if result.get('scores_xlsx'):
            logger.info(f"  Excel:   {result['scores_xlsx']}")
        logger.info(f"  JSON:    {result['scores_json']}")
        logger.info(f"  Summary: {result['summary_txt']}")
        logger.info(f"  Log:     {log_path}")
        logger.info("=" * 60)
        logger.info("M02c COMPLETE — Vina pipeline finished")
        logger.info("=" * 60)

        return 0

    except Exception as e:
        logger.error(f"Failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
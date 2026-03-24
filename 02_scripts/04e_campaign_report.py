#!/usr/bin/env python3
"""
04e DOCK6 Campaign Report - CLI
===================================
HTML report with composite ranking from 04a-04d.

Input:  05_results/{campaign}/04_dock6_analysis/04a-04d outputs
Output: 05_results/{campaign}/04_dock6_analysis/04e_campaign_report/

Project: molecular_docking
Module: 04e (DOCK6 analysis)
Version: 1.1 (2026-03-23) — path fix: 05_dock6 → 04_dock6_analysis
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from molecular_docking.m04_dock6_analysis.campaign_report import run_campaign_report

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
        description="04e DOCK6 Campaign Report — HTML report + composite ranking",
    )
    parser.add_argument("--config", "-c", type=str, help="Module config YAML")
    parser.add_argument("--campaign", type=str, help="Campaign config YAML")
    parser.add_argument("--results-base", type=str, default=None,
                        help="Base dir with 04a-04d outputs (04_dock6_analysis/)")
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--log-level", type=str, default=None)

    args = parser.parse_args()

    results_base = None
    output_dir = None
    campaign_id = "direct"
    composite_weights = None
    log_level = "INFO"

    if args.campaign:
        cc = load_yaml(args.campaign)
        campaign_id = cc.get("campaign_id", Path(args.campaign).parent.name)
        base = Path("05_results") / campaign_id / "04_dock6_analysis"
        results_base = str(base)
        output_dir = str(base / "04e_campaign_report")

    if args.config:
        mc = load_yaml(args.config)
        params = mc.get("parameters", {})
        log_level = params.get("log_level", log_level)
        if "composite_weights" in params:
            composite_weights = params["composite_weights"]

    if args.results_base:
        results_base = args.results_base
    if args.output:
        output_dir = args.output
    if args.log_level:
        log_level = args.log_level

    if not results_base or not Path(results_base).exists():
        logger.error(f"Results base not found: {results_base}")
        logger.error("Run modules 04a-04d first.")
        return 1

    if not output_dir:
        output_dir = str(Path(results_base) / "04e_campaign_report")

    logging.getLogger().setLevel(getattr(logging, log_level.upper(), logging.INFO))
    setup_log_file(Path(output_dir) / "04e_campaign_report.log", log_level)

    logger.info("=" * 60)
    logger.info("  MOLECULAR_DOCKING — Module 04e: DOCK6 Campaign Report")
    logger.info("=" * 60)
    logger.info(f"  Campaign: {campaign_id}")

    result = run_campaign_report(
        results_base=results_base,
        output_dir=output_dir,
        campaign_id=campaign_id,
        composite_weights=composite_weights,
    )

    if not result.get("success"):
        logger.error(f"Error: {result.get('error')}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

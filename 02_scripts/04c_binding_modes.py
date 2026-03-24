#!/usr/bin/env python3
"""
04c DOCK6 Binding Modes - CLI
=================================
Characterize DOCK6 binding modes per molecule.

Input:  05_results/{campaign}/01c_dock6_run/
Output: 05_results/{campaign}/04_dock6_analysis/04c_binding_modes/

Project: molecular_docking
Module: 04c (DOCK6 analysis)
Version: 1.1 (2026-03-23) — path fix: 05_dock6 → 04_dock6_analysis
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from molecular_docking.m04_dock6_analysis.binding_modes import run_binding_modes

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
        description="04c DOCK6 Binding Modes — mode characterization + RMSD",
    )
    parser.add_argument("--config", "-c", type=str, help="Module config YAML")
    parser.add_argument("--campaign", type=str, help="Campaign config YAML")
    parser.add_argument("--docking-dir", type=str, default=None)
    parser.add_argument("--footprint-csv", type=str, default=None,
                        help="footprint_per_molecule.csv from 04b")
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--log-level", type=str, default=None)

    args = parser.parse_args()

    docking_dir = None
    output_dir = None
    footprint_csv = None
    campaign_id = "direct"
    log_level = "INFO"

    if args.campaign:
        cc = load_yaml(args.campaign)
        campaign_id = cc.get("campaign_id", Path(args.campaign).parent.name)
        base = Path("05_results") / campaign_id
        docking_dir = str(base / "01c_dock6_run")
        output_dir = str(base / "04_dock6_analysis" / "04c_binding_modes")
        footprint_csv = str(base / "04_dock6_analysis" / "04b_footprint_analysis" / "footprint_per_molecule.csv")

    if args.config:
        mc = load_yaml(args.config)
        log_level = mc.get("parameters", {}).get("log_level", log_level)

    if args.docking_dir:
        docking_dir = args.docking_dir
    if args.footprint_csv:
        footprint_csv = args.footprint_csv
    if args.output:
        output_dir = args.output
    if args.log_level:
        log_level = args.log_level

    if not docking_dir or not Path(docking_dir).exists():
        logger.error(f"Docking dir not found: {docking_dir}")
        return 1

    if not output_dir:
        output_dir = "05_results/04_dock6_analysis/04c_binding_modes"

    if log_level:
        logging.getLogger().setLevel(getattr(logging, log_level.upper(), logging.INFO))
    setup_log_file(Path(output_dir) / "04c_binding_modes.log", log_level)

    logger.info("=" * 60)
    logger.info("  MOLECULAR_DOCKING — Module 04c: DOCK6 Binding Modes")
    logger.info("=" * 60)
    logger.info(f"  Campaign: {campaign_id}")

    # Only pass footprint if it exists
    fp_csv = footprint_csv if footprint_csv and Path(footprint_csv).exists() else None

    result = run_binding_modes(
        docking_dir=docking_dir,
        output_dir=output_dir,
        footprint_csv=fp_csv,
    )

    if not result.get("success"):
        logger.error(f"Error: {result.get('error')}")
        return 1

    logger.info(f"\nNext: python 02_scripts/04d_contact_mapping.py "
                f"--config 03_configs/04d_contact_mapping.yaml "
                f"--campaign {args.campaign or '<campaign_config.yaml>'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

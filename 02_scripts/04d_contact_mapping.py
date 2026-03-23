#!/usr/bin/env python3
"""
04d DOCK6 Contact Mapping - CLI
===================================
Distance-based contacts cross-referenced with footprint energies.

Input:  01c_dock6_run/, 00b receptor, 04b footprint (optional)
Output: 05_results/{campaign}/05_dock6/04d_contact_mapping/

Project: molecular_docking
Module: 04d (DOCK6 analysis)
Version: 1.0 (2026-03-22)
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from molecular_docking.m04_dock6_analysis.contact_mapping import run_contact_mapping

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
        description="04d DOCK6 Contact Mapping — contacts + footprint cross-reference",
    )
    parser.add_argument("--config", "-c", type=str, help="Module config YAML")
    parser.add_argument("--campaign", type=str, help="Campaign config YAML")
    parser.add_argument("--docking-dir", type=str, default=None)
    parser.add_argument("--receptor", type=str, default=None,
                        help="Receptor PDB or mol2")
    parser.add_argument("--footprint-csv", type=str, default=None)
    parser.add_argument("--all-poses", type=str, default=None)
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--contact-cutoff", type=float, default=None)
    parser.add_argument("--log-level", type=str, default=None)

    args = parser.parse_args()

    docking_dir = None
    receptor_path = None
    output_dir = None
    footprint_csv = None
    all_poses_csv = None
    campaign_id = "direct"
    contact_cutoff = 4.5
    log_level = "INFO"

    if args.campaign:
        cc = load_yaml(args.campaign)
        campaign_id = cc.get("campaign_id", Path(args.campaign).parent.name)
        base = Path("05_results") / campaign_id
        docking_dir = str(base / "01c_dock6_run")
        output_dir = str(base / "05_dock6" / "04d_contact_mapping")
        all_poses_csv = str(base / "01d_score_collection" / "dock6_all_poses.csv")
        footprint_csv = str(base / "05_dock6" / "04b_footprint_analysis" / "footprint_per_molecule.csv")

        # Receptor: try PDB first (from 00b), then mol2
        rec_pdb = base / "00b_receptor_preparation" / "rec_noH.pdb"
        rec_mol2 = base / "00b_receptor_preparation" / "rec_charged.mol2"
        if rec_pdb.exists():
            receptor_path = str(rec_pdb)
        elif rec_mol2.exists():
            receptor_path = str(rec_mol2)

    if args.config:
        mc = load_yaml(args.config)
        params = mc.get("parameters", {})
        contact_cutoff = params.get("contact_cutoff", contact_cutoff)
        log_level = params.get("log_level", log_level)

    if args.docking_dir:
        docking_dir = args.docking_dir
    if args.receptor:
        receptor_path = args.receptor
    if args.footprint_csv:
        footprint_csv = args.footprint_csv
    if args.all_poses:
        all_poses_csv = args.all_poses
    if args.output:
        output_dir = args.output
    if args.contact_cutoff is not None:
        contact_cutoff = args.contact_cutoff
    if args.log_level:
        log_level = args.log_level

    if not docking_dir or not Path(docking_dir).exists():
        logger.error(f"Docking dir not found: {docking_dir}")
        return 1
    if not receptor_path or not Path(receptor_path).exists():
        logger.error(f"Receptor not found: {receptor_path}")
        return 1
    if not output_dir:
        output_dir = "05_results/05_dock6/04d_contact_mapping"

    if log_level:
        logging.getLogger().setLevel(getattr(logging, log_level.upper(), logging.INFO))
    setup_log_file(Path(output_dir) / "04d_contact_mapping.log", log_level)

    logger.info("=" * 60)
    logger.info("  MOLECULAR_DOCKING — Module 04d: DOCK6 Contact Mapping")
    logger.info("=" * 60)
    logger.info(f"  Campaign: {campaign_id}")

    fp_csv = footprint_csv if footprint_csv and Path(footprint_csv).exists() else None
    ap_csv = all_poses_csv if all_poses_csv and Path(all_poses_csv).exists() else None

    result = run_contact_mapping(
        docking_dir=docking_dir,
        receptor_path=receptor_path,
        output_dir=output_dir,
        footprint_csv=fp_csv,
        contact_cutoff=contact_cutoff,
        all_poses_csv=ap_csv,
    )

    if not result.get("success"):
        logger.error(f"Error: {result.get('error')}")
        return 1

    logger.info(f"\nNext: python 02_scripts/04e_campaign_report.py "
                f"--config 03_configs/04e_campaign_report.yaml "
                f"--campaign {args.campaign or '<campaign_config.yaml>'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

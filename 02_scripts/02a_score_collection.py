#!/usr/bin/env python3
"""
02a Score Collection - CLI
============================
Collects DOCK6 docking scores, ranks molecules, and produces Excel output.

Scans 01b output for scored mol2 files, extracts Grid_Score per pose,
selects best pose per molecule, and generates ranked Excel + CSV.

Reads campaign_config.yaml:
    - campaign_id

Reads module YAML (03_configs/02a_score_collection.yaml):
    - score_key, max_molecules, extract_best_pose_mol2, compute_properties

Usage (PyCharm Pro terminal — single line):
    python 02_scripts/02a_score_collection.py --config 03_configs/02a_score_collection.yaml --campaign 04_data/campaigns/example_pH72/campaign_config.yaml

PyCharm Run Configuration:
    Script:     02_scripts/02a_score_collection.py
    Parameters: --config 03_configs/02a_score_collection.yaml --campaign 04_data/campaigns/example_campaign/campaign_config.yaml
    Working dir: (project root)

Project: molecular_docking
Module: 02a
Version: 1.0
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-8s | %(message)s",
)

sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from molecular_docking.m02_collection.score_collector import run_score_collection

logger = logging.getLogger(__name__)


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def setup_log_file(log_path: Path, log_level: str = "INFO"):
    """Add file handler to root logger."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(str(log_path), encoding="utf-8")
    fh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    fh.setFormatter(logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s"))
    logging.getLogger().addHandler(fh)


def main():
    parser = argparse.ArgumentParser(
        description="Collect DOCK6 scores and produce ranked Excel",
    )
    parser.add_argument("--config", "-c", type=str, required=True,
                        help="Module config YAML")
    parser.add_argument("--campaign", type=str, required=True,
                        help="Campaign config YAML")
    parser.add_argument("--output", "-o", type=str, default=None,
                        help="Output directory")
    parser.add_argument("--docking-dir", type=str, default=None,
                        help="Override docking directory (default: 01b output)")
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    args = parser.parse_args()

    # --- Load configs ---
    cc = load_yaml(args.campaign)
    campaign_id = cc.get("campaign_id", Path(args.campaign).parent.name)
    mc = load_yaml(args.config)
    params = mc.get("parameters", {})

    # --- Resolve paths ---
    docking_dir = args.docking_dir or str(
        Path("05_results") / campaign_id / "01b_dock6_run"
    )
    molecules_csv = str(
        Path("05_results") / campaign_id / "00a_molecule_parser" / "unique_molecules.csv"
    )
    output_dir = args.output or str(
        Path("05_results") / campaign_id / "02a_score_collection"
    )

    # --- Params ---
    score_key = params.get("score_key", "Grid_Score")
    max_molecules = params.get("max_molecules", 500)
    extract_best = params.get("extract_best_pose_mol2", True)
    keep_all = params.get("keep_all_poses", False)
    compute_props = params.get("compute_properties", True)
    log_level = args.log_level or params.get("log_level", "INFO")

    # --- Output filenames from campaign config ---
    oc = cc.get("output", {})
    scores_filename = oc.get("scores_filename", "dock6_scores.xlsx")
    mol2_dirname = oc.get("mol2_dirname", "best_poses")

    # --- Setup logging ---
    logging.getLogger().setLevel(getattr(logging, log_level.upper(), logging.INFO))
    log_path = Path(output_dir) / "02a_score_collection.log"
    setup_log_file(log_path, log_level)

    # --- Execute ---
    logger.info("=" * 60)
    logger.info("  MOLECULAR_DOCKING - Module 02a: Score Collection")
    logger.info("=" * 60)
    logger.info(f"Campaign:     {campaign_id}")
    logger.info(f"Docking dir:  {docking_dir}")
    logger.info(f"Score key:    {score_key}")
    logger.info(f"Output:       {output_dir}")

    result = run_score_collection(
        docking_dir=docking_dir,
        output_dir=output_dir,
        molecules_csv=molecules_csv,
        score_key=score_key,
        max_molecules=max_molecules,
        extract_best_pose_mol2=extract_best,
        keep_all_poses=keep_all,
        compute_properties=compute_props,
        scores_filename=scores_filename,
        mol2_dirname=mol2_dirname,
        source_label=cc.get("metadata", {}).get("source"),
    )

    if not result.get("success"):
        logger.error(f"Error: {result.get('error')}")
        return 1

    logger.info("")
    logger.info(f"Next steps:")
    logger.info(f"  - Open {result.get('output_xlsx')} to review scores")
    logger.info(f"  - Best poses in {result.get('mol2_dir')}/")

    return 0


if __name__ == "__main__":
    sys.exit(main())
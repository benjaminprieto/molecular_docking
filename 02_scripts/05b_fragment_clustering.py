#!/usr/bin/env python3
"""
05b Fragment Clustering - CLI
===============================
DBSCAN clustering per rigid fragment per molecule.

Reads:  05_results/{campaign_id}/05a_parse_and_fragment/
Output: 05_results/{campaign_id}/05b_fragment_clustering/

Usage:
    python 02_scripts/05b_fragment_clustering.py \
        --config 03_configs/05b_fragment_clustering.yaml \
        --campaign 04_data/campaigns/UDX_pharmit_pH63/campaign_config.yaml

Project: molecular_docking
Module: 05b
Version: 2.0
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)-8s | %(message)s")

sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from molecular_docking.m05_gnina_analysis.fragment_clustering import run_fragment_clustering

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
        description="05b Fragment Clustering — DBSCAN per fragment per molecule",
    )
    parser.add_argument("--config", "-c", type=str, required=True)
    parser.add_argument("--campaign", type=str, required=True)
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--eps", type=float, default=None,
                        help="DBSCAN eps (Angstrom)")
    parser.add_argument("--min-samples", type=int, default=None)
    parser.add_argument("--name", type=str, default=None)
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    args = parser.parse_args()

    cc = load_yaml(args.campaign)
    mc = load_yaml(args.config)
    params = mc.get("parameters", {})

    campaign_dir = Path(args.campaign).parent
    campaign_id = cc.get("campaign_id", campaign_dir.name)
    results_base = Path("05_results") / campaign_id / "m05_gnina_analysis"

    output_subdir = mc.get("outputs", {}).get("subdir", "05b_fragment_clustering")
    output_dir = Path(args.output) if args.output else results_base / output_subdir

    log_level = args.log_level or params.get("log_level", "INFO")
    setup_log_file(output_dir / "05b_fragment_clustering.log", log_level)

    eps = args.eps or params.get("eps", 3.0)
    min_samples = args.min_samples or params.get("min_samples", 5)
    min_dominant_fraction = params.get("min_dominant_fraction", 0.15)

    parsed_dir = str(results_base / "05a_parse_and_fragment")
    molecule_names = [args.name] if args.name else None

    logger.info(f"Campaign:    {campaign_id}")
    logger.info(f"eps:         {eps} A")
    logger.info(f"min_samples: {min_samples}")
    logger.info(f"Input:       {parsed_dir}")
    logger.info(f"Output:      {output_dir}")

    result = run_fragment_clustering(
        parsed_dir=parsed_dir,
        output_dir=str(output_dir),
        eps=eps,
        min_samples=min_samples,
        min_dominant_fraction=min_dominant_fraction,
        molecule_names=molecule_names,
    )

    if not result.get("success"):
        logger.error(f"Failed: {result.get('error')}")
        sys.exit(1)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
05c Score Decomposition - CLI
===============================
Regress docking scores onto fragment cluster positions.
Identifies per-fragment contributions and sweet spot combinations.

Reads:  05a (parse_and_fragment) + 05b (fragment_clustering)
Output: 05_results/{campaign_id}/05c_score_decomposition/

Usage:
    python 02_scripts/05c_score_decomposition.py \
        --config 03_configs/05c_score_decomposition.yaml \
        --campaign 04_data/campaigns/UDX_pharmit_pH63/campaign_config.yaml

Project: molecular_docking
Module: 05c
Version: 3.0
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)-8s | %(message)s")

sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from molecular_docking.m05_gnina_analysis.score_decomposition import run_score_decomposition

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
        description="05c Score Decomposition — per-fragment score contributions + sweet spots",
    )
    parser.add_argument("--config", "-c", type=str, required=True)
    parser.add_argument("--campaign", type=str, required=True)
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--score-key", type=str, default=None,
                        help="Decompose only this score (e.g. vina_affinity)")
    parser.add_argument("--name", type=str, default=None)
    parser.add_argument("--gnina-scores", type=str, default=None,
                        help="Path to gnina_scores.csv from 02c (auto-detected if omitted)")
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    args = parser.parse_args()

    cc = load_yaml(args.campaign)
    mc = load_yaml(args.config)
    params = mc.get("parameters", {})

    campaign_dir = Path(args.campaign).parent
    campaign_id = cc.get("campaign_id", campaign_dir.name)
    results_base = Path("05_results") / campaign_id / "05_gnina_analysis"

    output_subdir = mc.get("outputs", {}).get("subdir", "05c_score_decomposition")
    output_dir = Path(args.output) if args.output else results_base / output_subdir

    log_level = args.log_level or params.get("log_level", "INFO")
    setup_log_file(output_dir / "05c_score_decomposition.log", log_level)

    n_sweet_spots = params.get("n_sweet_spots", 10)
    score_keys = [args.score_key] if args.score_key else params.get("score_keys", None)

    parsed_dir = str(results_base / "05a_parse_and_fragment")
    cluster_dir = str(results_base / "05b_fragment_clustering")
    molecule_names = [args.name] if args.name else None

    # Resolve gnina_scores.csv from 02c
    gnina_scores_csv = args.gnina_scores or params.get("gnina_scores_csv")
    if not gnina_scores_csv:
        auto_path = Path("05_results") / campaign_id / "02c_gnina_scores" / "gnina_scores.csv"
        if auto_path.exists():
            gnina_scores_csv = str(auto_path)
            logger.info(f"  Auto-detected gnina_scores.csv: {gnina_scores_csv}")

    logger.info(f"Campaign:     {campaign_id}")
    logger.info(f"Score keys:   {score_keys or 'all'}")
    logger.info(f"Sweet spots:  top {n_sweet_spots}")
    logger.info(f"Input (05a):  {parsed_dir}")
    logger.info(f"Input (05b):  {cluster_dir}")
    logger.info(f"Output:       {output_dir}")

    result = run_score_decomposition(
        parsed_dir=parsed_dir,
        cluster_dir=cluster_dir,
        output_dir=str(output_dir),
        score_keys=score_keys,
        n_sweet_spots=n_sweet_spots,
        molecule_names=molecule_names,
        gnina_scores_csv=gnina_scores_csv,
    )

    if not result.get("success"):
        logger.error(f"Failed: {result.get('error')}")
        sys.exit(1)


if __name__ == "__main__":
    main()

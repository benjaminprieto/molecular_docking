#!/usr/bin/env python3
"""
05a Parse & Fragment - CLI
============================
Parse docking poses + identify rigid fragments per molecule.

Reads upstream:
    DOCK6: 05_results/{campaign_id}/01c_dock6_run/{name}/{name}_scored.mol2
    GNINA: 05_results/{campaign_id}/02b_gnina_run/poses/{name}_gnina.sdf

Output: 05_results/{campaign_id}/05a_parse_and_fragment/{name}.npz + {name}.json

Usage:
    python 02_scripts/05a_parse_and_fragment.py \\
        --config 03_configs/05a_parse_and_fragment.yaml \\
        --campaign 04_data/campaigns/UDX_pharmit_pH63/campaign_config.yaml

Project: molecular_docking
Module: 05a
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

from molecular_docking.m05_gnina_analysis.parse_and_fragment import run_parse_and_fragment

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
        description="05a Parse & Fragment — parse poses + rigid fragment identification",
    )
    parser.add_argument("--config", "-c", type=str, required=True)
    parser.add_argument("--campaign", type=str, required=True)
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--engine", type=str, default=None,
                        choices=["dock6", "gnina"])
    parser.add_argument("--min-fragment-size", type=int, default=None)
    parser.add_argument("--name", type=str, default=None,
                        help="Process only this molecule")
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    args = parser.parse_args()

    cc = load_yaml(args.campaign)
    mc = load_yaml(args.config)
    params = mc.get("parameters", {})

    campaign_dir = Path(args.campaign).parent
    campaign_id = cc.get("campaign_id", campaign_dir.name)
    results_base = Path("05_results") / campaign_id

    engine = args.engine or params.get("engine", "gnina")
    engine_base = results_base / "m05_gnina_analysis"
    output_subdir = mc.get("outputs", {}).get("subdir", "05a_parse_and_fragment")
    output_dir = Path(args.output) if args.output else engine_base / output_subdir

    log_level = args.log_level or params.get("log_level", "INFO")
    setup_log_file(output_dir / "05a_parse_and_fragment.log", log_level)

    engine = args.engine or params.get("engine", "gnina")
    min_fragment_size = args.min_fragment_size or params.get("min_fragment_size", 4)

    dock6_dir = str(results_base / "01c_dock6_run")
    gnina_poses_dir = str(results_base / "02b_gnina_run" / "poses")
    molecules_csv = str(results_base / "00a_molecule_parser" / "unique_molecules.csv")
    antechamber_dir = str(results_base / "01a_antechamber" / "mol2")
    ionization_sdf_dir = str(results_base / "00c_ionization_profiling" / "sdf")

    molecule_names = [args.name] if args.name else None

    logger.info(f"Campaign:     {campaign_id}")
    logger.info(f"Engine:       {engine}")
    logger.info(f"min_fragment: {min_fragment_size} heavy atoms")
    logger.info(f"Output:       {output_dir}")

    result = run_parse_and_fragment(
        output_dir=str(output_dir),
        engine=engine,
        dock6_dir=dock6_dir,
        gnina_poses_dir=gnina_poses_dir,
        molecules_csv=molecules_csv,
        antechamber_dir=antechamber_dir,
        ionization_sdf_dir=ionization_sdf_dir,
        min_fragment_size=min_fragment_size,
        molecule_names=molecule_names,
    )

    if not result.get("success"):
        logger.error(f"Failed: {result.get('error')}")
        sys.exit(1)


if __name__ == "__main__":
    main()

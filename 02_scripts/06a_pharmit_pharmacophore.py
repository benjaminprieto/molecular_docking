#!/usr/bin/env python3
"""
06a Pharmit Pharmacophore Generator - CLI (v4.0)
===================================================
PLIP IS the pharmacophore. No RDKit intermediate.

Usage:
    python 02_scripts/06a_pharmit_pharmacophore.py \
        --plip 05_results/.../03a_plip_analysis/interactions.json \
        --dock6 05_results/.../04_dock6_analysis/04b_footprint_analysis/residue_consensus.csv

    # With campaign (auto-resolves):
    python 02_scripts/06a_pharmit_pharmacophore.py \
        --campaign 04_data/campaigns/UDX_reference_pH63/campaign_config.yaml \
        --config 03_configs/06a_pharmit_pharmacophore.yaml

Project: molecular_docking
Module: 06a
Version: 4.0
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)-8s | %(message)s")

sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from molecular_docking.m06_pharmit.pharmit_pharmacophore import (
    generate_pharmit_pharmacophore,
)

logger = logging.getLogger(__name__)


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def setup_log_file(log_path: Path, log_level: str = "INFO"):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(str(log_path), mode="w", encoding="utf-8")
    fh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    fh.setFormatter(logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s"))
    logging.getLogger().addHandler(fh)


def main():
    parser = argparse.ArgumentParser(
        description="06a Pharmit Pharmacophore — PLIP → sub-pocket Pharmit queries",
    )

    parser.add_argument("--plip", type=str, default=None,
                        help="PLIP interactions JSON (from 03a)")
    parser.add_argument("--dock6", type=str, default=None,
                        help="DOCK6 residue_consensus.csv (from 04b)")
    parser.add_argument("--ligand", "-l", type=str, default=None,
                        help="Ligand SDF/mol2 (optional, for GNINA rescore)")
    parser.add_argument("--receptor", "-r", type=str, default=None,
                        help="Receptor PDB (optional, for GNINA rescore)")

    parser.add_argument("--campaign", type=str, default=None)
    parser.add_argument("--config", "-c", type=str, default=None)

    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--name", type=str, default="pharmacophore")
    parser.add_argument("--n-required", type=int, default=None)
    parser.add_argument("--n-optional", type=int, default=None)
    parser.add_argument("--radius", type=float, default=None)
    parser.add_argument("--energy-cutoff", type=float, default=None)
    parser.add_argument("--group-tolerance", type=float, default=None)
    parser.add_argument("--strategies", type=str, default=None,
                        help="Comma-separated: xylose,uracil,combined,analogues,druglike")
    parser.add_argument("--include-hydrophobic", action="store_true")
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    args = parser.parse_args()

    # Defaults
    n_required = 3
    n_optional = 2
    radius = 1.0
    energy_cutoff = 0.0
    group_tolerance = 1.0
    include_hydrophobic = False
    output_name = args.name
    strategies = None
    log_level = "INFO"

    if args.config:
        mc = load_yaml(args.config)
        params = mc.get("parameters", {})
        n_required = params.get("n_required", n_required)
        n_optional = params.get("n_optional", n_optional)
        radius = params.get("radius", radius)
        energy_cutoff = params.get("energy_cutoff", energy_cutoff)
        group_tolerance = params.get("group_tolerance", group_tolerance)
        include_hydrophobic = params.get("include_hydrophobic", include_hydrophobic)
        output_name = params.get("output_name", output_name)
        strategies = params.get("strategies", strategies)
        log_level = params.get("log_level", log_level)

    # CLI overrides
    if args.n_required is not None: n_required = args.n_required
    if args.n_optional is not None: n_optional = args.n_optional
    if args.radius is not None: radius = args.radius
    if args.energy_cutoff is not None: energy_cutoff = args.energy_cutoff
    if args.group_tolerance is not None: group_tolerance = args.group_tolerance
    if args.include_hydrophobic: include_hydrophobic = True
    if args.strategies: strategies = [s.strip() for s in args.strategies.split(",")]
    if args.log_level: log_level = args.log_level

    logging.getLogger().setLevel(getattr(logging, log_level))

    # Resolve paths
    plip_json = args.plip
    dock6_csv = args.dock6
    ligand_path = args.ligand
    receptor_path = args.receptor
    output_dir = args.output

    if args.campaign:
        cc = load_yaml(args.campaign)
        campaign_dir = Path(args.campaign).parent
        campaign_id = cc.get("campaign_id", campaign_dir.name)
        results_base = Path("05_results") / campaign_id

        # PLIP JSON from 03a
        if plip_json is None:
            candidate = results_base / "03a_plip_analysis" / "interactions.json"
            if candidate.exists():
                plip_json = str(candidate)
                logger.info(f"  Auto-detected PLIP: {candidate}")

        # DOCK6 footprint from 04b
        if dock6_csv is None:
            for subdir in ["04_dock6_analysis/04b_footprint_analysis",
                           "04b_footprint_analysis"]:
                candidate = results_base / subdir / "residue_consensus.csv"
                if candidate.exists():
                    dock6_csv = str(candidate)
                    logger.info(f"  Auto-detected DOCK6: {candidate}")
                    break

        # Ligand for GNINA
        if ligand_path is None:
            ref = cc.get("grids", {}).get("binding_site", {}).get("reference_mol2")
            if ref:
                candidate = campaign_dir / ref
                if candidate.exists():
                    ligand_path = str(candidate)

        # Receptor for GNINA
        if receptor_path is None:
            for name in ["receptor_clean.pdb", "rec_noH.pdb"]:
                candidate = results_base / "00b_receptor_preparation" / name
                if candidate.exists():
                    receptor_path = str(candidate)
                    break

        # Output
        if output_dir is None:
            output_dir = str(results_base / "06a_pharmit")

    # Validate
    if plip_json is None or not Path(plip_json).exists():
        logger.error(f"PLIP JSON not found: {plip_json}")
        logger.error("Run 03a first, or provide --plip")
        sys.exit(1)

    if output_dir is None:
        output_dir = "05_results/06a_pharmit"

    # Log file
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    setup_log_file(out_path / "06a_pharmit.log", log_level)

    # Run
    result = generate_pharmit_pharmacophore(
        plip_json_path=plip_json,
        output_dir=output_dir,
        output_name=output_name,
        dock6_csv_path=dock6_csv,
        ligand_path=ligand_path,
        receptor_path=receptor_path,
        include_hydrophobic=include_hydrophobic,
        energy_cutoff=energy_cutoff,
        n_required=n_required,
        n_optional=n_optional,
        radius=radius,
        group_tolerance=group_tolerance,
        strategies=strategies,
    )

    if result["success"]:
        logger.info(f"\nUpload to pharmit.csb.pitt.edu:")
        for strat, path in result["queries"].items():
            logger.info(f"  {strat}: {path}")
        sys.exit(0)
    else:
        logger.error(f"Failed: {result.get('error')}")
        sys.exit(1)


if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""
01a Antechamber Preparation - CLI (DOCK6-specific)
=====================================================
Converts protonated SDF files from 00c into DOCK6-ready mol2 files via antechamber.

This module is DOCK6-specific — Vina does not need mol2 files.

Input: Individual SDF files from 00c (structures/pH{value}/*.sdf)
Output: mol2/*.mol2 (AM1-BCC charges, GAFF2 atom types)

Reads campaign_config.yaml:
    - docking_ph          -> selects pH folder from 00c

Reads module YAML (03_configs/01a_antechamber_preparation.yaml):
    - charge_method, atom_type, timeout, obabel_fallback

Hardcoded upstream paths:
    - 05_results/{campaign_id}/00c_ionization_profiling/.../structures/pH{xx}/
    - Output: 05_results/{campaign_id}/01a_antechamber/

Usage:
    python 02_scripts/01a_antechamber_preparation.py --config 03_configs/01a_antechamber_preparation.yaml --campaign 04_data/campaigns/example_campaign/campaign_config.yaml

Project: molecular_docking
Module: 01a (DOCK6 engine)
Version: 1.2 — renumbered from 00d (2026-03-16)
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

from molecular_docking.m01_docking.antechamber_preparation import run_antechamber_preparation

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


def _find_latest_run(ionization_dir: Path) -> Path:
    """Find the latest run_id directory in 00c output."""
    run_dirs = [
        d for d in ionization_dir.iterdir()
        if d.is_dir() and d.name not in ("_work",) and not d.name.startswith(".")
    ]
    if not run_dirs:
        return ionization_dir
    run_dirs.sort(key=lambda d: d.name, reverse=True)
    return run_dirs[0]


def _find_ph_folder(structures_dir: Path, docking_ph: float) -> Path:
    """Find the pH folder matching docking_ph."""
    ph_label = f"pH{int(round(docking_ph * 10))}"
    ph_dir = structures_dir / ph_label

    if ph_dir.exists():
        return ph_dir

    available = [
        d.name for d in structures_dir.iterdir()
        if d.is_dir() and d.name.startswith("pH")
    ]

    if available:
        logger.warning(f"pH folder '{ph_label}' not found. "
                       f"Available: {sorted(available)}")
        if len(available) == 1:
            logger.info(f"  Using only available folder: {available[0]}")
            return structures_dir / available[0]

    return ph_dir


def main():
    parser = argparse.ArgumentParser(
        description="01a Antechamber — Convert protonated SDFs to DOCK6-ready mol2 (DOCK6 engine)",
    )
    parser.add_argument("--config", "-c", type=str, default=None,
                        help="Module config YAML")
    parser.add_argument("--campaign", type=str, default=None,
                        help="Campaign config YAML")

    parser.add_argument("--sdf-dir", type=str, default=None,
                        help="Path to SDF directory (overrides auto-detection)")
    parser.add_argument("--output", "-o", type=str, default=None,
                        help="Output directory")

    parser.add_argument("--charge-method", type=str, default=None,
                        choices=["bcc", "gas"])
    parser.add_argument("--atom-type", type=str, default=None,
                        choices=["gaff2", "gaff"])
    parser.add_argument("--timeout", type=int, default=None,
                        help="Timeout per molecule (seconds)")
    parser.add_argument("--no-fallback", action="store_true",
                        help="Disable obabel fallback")
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    args = parser.parse_args()

    sdf_dir = None
    output_dir = None
    campaign_id = "direct"
    docking_ph = 7.2
    charge_method = "bcc"
    atom_type = "gaff2"
    obabel_fallback = True
    timeout_per_molecule = 300
    log_level = "INFO"

    if args.campaign:
        cc = load_yaml(args.campaign)
        campaign_dir = Path(args.campaign).parent
        campaign_id = cc.get("campaign_id", campaign_dir.name)
        docking_ph = cc.get("docking_ph", docking_ph)

        ionization_dir = Path("05_results") / campaign_id / "00c_ionization_profiling"
        if ionization_dir.exists():
            latest_run = _find_latest_run(ionization_dir)
            structures = latest_run / "structures"
            if structures.exists():
                sdf_dir = str(_find_ph_folder(structures, docking_ph))

        # --- CHANGED: 00d_antechamber → 01a_antechamber ---
        output_dir = str(Path("05_results") / campaign_id / "01a_antechamber")

    if args.config:
        mc = load_yaml(args.config)
        params = mc.get("parameters", {})
        charge_method = params.get("charge_method", charge_method)
        atom_type = params.get("atom_type", atom_type)
        obabel_fallback = params.get("obabel_fallback", obabel_fallback)
        timeout_per_molecule = params.get("timeout_per_molecule", timeout_per_molecule)
        log_level = params.get("log_level", log_level)

    if args.sdf_dir:
        sdf_dir = args.sdf_dir
    if args.output:
        output_dir = args.output
    if args.charge_method:
        charge_method = args.charge_method
    if args.atom_type:
        atom_type = args.atom_type
    if args.timeout:
        timeout_per_molecule = args.timeout
    if args.no_fallback:
        obabel_fallback = False
    if args.log_level:
        log_level = args.log_level

    if not sdf_dir:
        parser.error("Provide --campaign or --sdf-dir. "
                     "Run 00c first to generate protonated SDF files.")
    if not output_dir:
        output_dir = "05_results/01a_antechamber"

    if not Path(sdf_dir).exists():
        logger.error(f"SDF directory not found: {sdf_dir}")
        logger.error("Run module 00c first.")
        return 1

    n_sdf = len(list(Path(sdf_dir).glob("*.sdf")))
    if n_sdf == 0:
        logger.error(f"No SDF files found in: {sdf_dir}")
        return 1

    if log_level:
        logging.getLogger().setLevel(getattr(logging, log_level.upper(), logging.INFO))

    # --- CHANGED: log file name ---
    log_path = Path(output_dir) / "01a_antechamber.log"
    setup_log_file(log_path, log_level)

    logger.info("=" * 60)
    logger.info("  MOLECULAR_DOCKING - Module 01a: Antechamber (DOCK6)")
    logger.info("=" * 60)
    logger.info(f"Campaign:     {campaign_id}")
    logger.info(f"SDF dir:      {sdf_dir} ({n_sdf} files)")
    logger.info(f"Output:       {output_dir}")
    logger.info(f"Charges:      {charge_method} ({atom_type})")
    logger.info(f"Timeout:      {timeout_per_molecule}s per molecule")

    result = run_antechamber_preparation(
        sdf_dir=sdf_dir,
        output_dir=output_dir,
        charge_method=charge_method,
        atom_type=atom_type,
        obabel_fallback=obabel_fallback,
        timeout_per_molecule=timeout_per_molecule,
    )

    if not result.get("success"):
        logger.error(f"Error: {result.get('error')}")
        return 1

    logger.info("")
    # --- CHANGED: next step references ---
    logger.info(f"Next: python 02_scripts/01b_grid_generation.py "
                f"--config 03_configs/01b_grid_generation.yaml "
                f"--campaign {args.campaign or '<campaign_config.yaml>'}")

    return 0 if result["n_failed"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

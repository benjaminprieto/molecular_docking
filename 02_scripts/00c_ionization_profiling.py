#!/usr/bin/env python3
"""
00c Ionization Profiling - CLI
================================
Protona moléculas al pH de docking usando ionprofile.
Genera SDF protonados con 3D + Excel/CSV con perfiles de carga.

El usuario revisa los outputs antes de pasar a 00d (antechamber).

Reads campaign_config.yaml:
    - docking_ph                         -> pH o lista de pH
    - molecules.protonation.enabled      -> true/false
    - molecules.protonation.engine       -> openbabel | dimorphite | qupkake

Reads module YAML (03_configs/00c_ionization_profiling.yaml):
    - ph_range, ph_step, precision
    - output_formats

Input: 05_results/{campaign_id}/00a_molecule_parser/unique_molecules.csv
Output: 05_results/{campaign_id}/00c_ionization_profiling/{run_id}/
    - ionization_profiling.csv
    - ionization_profiling.xlsx
    - structures/pH72/*.sdf  (individuales, 3D, H explícitos)
    - ionprofile.log

Usage:
    python 02_scripts/00c_ionization_profiling.py --config 03_configs/00c_ionization_profiling.yaml --campaign 04_data/campaigns/example_campaign/campaign_config.yaml

Project: molecular_docking
Module: 00c
Version: 2.0
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

from molecular_docking.m00_preparation.ionization_profiling import run_ionization_profiling

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
        description="Protonate molecules at docking pH using ionprofile",
    )
    parser.add_argument("--config", "-c", type=str, default=None,
                        help="Module config YAML")
    parser.add_argument("--campaign", type=str, default=None,
                        help="Campaign config YAML")

    # Direct mode
    parser.add_argument("--input-csv", type=str, default=None,
                        help="Path to unique_molecules.csv")
    parser.add_argument("--output", "-o", type=str, default=None,
                        help="Output directory")

    # Overrides
    parser.add_argument("--ph", type=float, nargs="+", default=None,
                        help="Docking pH(s), e.g., --ph 7.2 or --ph 7.2 6.3")
    parser.add_argument("--engine", type=str, default=None,
                        choices=["openbabel", "dimorphite", "qupkake"])
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    args = parser.parse_args()

    # =========================================================================
    # RESOLVE PARAMETERS
    # =========================================================================

    molecules_csv = None
    output_dir = None
    campaign_id = "direct"
    docking_ph = 7.2
    ph_range = None
    ph_step = 0.1
    engine = "openbabel"
    precision = 0.5
    output_formats = ["csv", "excel", "sdf"]
    log_level = "INFO"

    # --- Campaign config ---
    if args.campaign:
        cc = load_yaml(args.campaign)
        campaign_dir = Path(args.campaign).parent
        campaign_id = cc.get("campaign_id", campaign_dir.name)
        docking_ph = cc.get("docking_ph", docking_ph)

        mc = cc.get("molecules", {})
        pc = mc.get("protonation", {})
        engine = pc.get("engine", engine)

        molecules_csv = str(
            Path("05_results") / campaign_id / "00a_molecule_parser" / "unique_molecules.csv"
        )
        output_dir = str(Path("05_results") / campaign_id / "00c_ionization_profiling")

    # --- Module config ---
    if args.config:
        mc = load_yaml(args.config)
        params = mc.get("parameters", {})
        ph_range = params.get("ph_range", ph_range)
        ph_step = params.get("ph_step", ph_step)
        precision = params.get("precision", precision)
        engine = params.get("engine", engine)
        output_formats = params.get("output_formats", output_formats)
        log_level = params.get("log_level", log_level)

    # --- CLI overrides ---
    if args.input_csv:
        molecules_csv = args.input_csv
    if args.output:
        output_dir = args.output
    if args.ph:
        if len(args.ph) == 1:
            docking_ph = args.ph[0]
        else:
            docking_ph = args.ph
    if args.engine:
        engine = args.engine
    if args.log_level:
        log_level = args.log_level

    if not molecules_csv:
        parser.error("Provide --campaign or --input-csv")
    if not output_dir:
        output_dir = "05_results/00c_ionization_profiling"

    if not Path(molecules_csv).exists():
        logger.error(f"Input CSV not found: {molecules_csv}")
        logger.error("Run module 00a first.")
        return 1

    # =========================================================================
    # SETUP LOGGING
    # =========================================================================

    if log_level:
        logging.getLogger().setLevel(getattr(logging, log_level.upper(), logging.INFO))

    log_path = Path(output_dir) / "00c_ionization_profiling.log"
    setup_log_file(log_path, log_level)

    # =========================================================================
    # EXECUTE
    # =========================================================================

    logger.info("=" * 60)
    logger.info("  MOLECULAR_DOCKING - Module 00c: Ionization Profiling")
    logger.info("=" * 60)
    logger.info(f"Campaign:     {campaign_id}")
    logger.info(f"Input:        {molecules_csv}")
    logger.info(f"Output:       {output_dir}")
    logger.info(f"pH:           {docking_ph}")
    logger.info(f"Engine:       {engine}")

    result = run_ionization_profiling(
        molecules_csv=molecules_csv,
        output_dir=output_dir,
        docking_ph=docking_ph,
        ph_range=ph_range,
        ph_step=ph_step,
        engine=engine,
        precision=precision,
        output_formats=output_formats,
    )

    if not result.get("success"):
        logger.error(f"Error: {result.get('error')}")
        return 1

    logger.info("")
    logger.info(f"Next: Review outputs in {result['output_dir']}")
    logger.info(f"  - ionization_profiling.csv/xlsx (charge profiles)")
    logger.info(f"  - structures/pH*/  (protonated SDF files)")
    logger.info(f"Then run 00d:")
    logger.info(f"  python 02_scripts/00d_antechamber.py "
                f"--config 03_configs/00d_antechamber.yaml "
                f"--campaign {args.campaign or '<campaign_config.yaml>'}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

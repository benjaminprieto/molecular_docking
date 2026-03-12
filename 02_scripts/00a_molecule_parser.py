#!/usr/bin/env python3
"""
00a Molecule Parser - CLI
===========================
Parsea moleculas de cualquier formato -> tabla normalizada + SDF limpio.

Formatos soportados: SDF, CSV/XLSX, mol2, SMILES txt, directorio mixto.

Reads campaign_config.yaml:
    - molecules.input_file       -> archivo o directorio de entrada
    - molecules.name_source      -> como extraer nombre del SDF
    - molecules.name_property    -> propiedad SDF para nombre
    - molecules.name_column      -> columna CSV/XLSX para nombre
    - molecules.smiles_column    -> columna CSV/XLSX para SMILES

Reads module YAML (03_configs/00a_molecule_parser.yaml):
    - conformer_strategy         -> best_rmsd | first | all
    - compute_properties         -> true/false
    - detect_smiles_duplicates   -> true/false

Output: 05_results/{campaign_id}/00a_molecule_parser/
    - unique_molecules.csv
    - unique_molecules.sdf
    - parse_summary.txt
    - 00a_molecule_parser.log

Usage:
    python 02_scripts/00a_molecule_parser.py \\
        --config 03_configs/00a_molecule_parser.yaml \\
        --campaign 04_data/campaigns/phermit_groove/campaign_config.yaml

    # Direct mode (no campaign):
    python 02_scripts/00a_molecule_parser.py \\
        --input molecules.sdf \\
        --output results/00a/

Project: molecular_docking
Module: 00a
Version: 1.0
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path
from datetime import datetime

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-8s | %(message)s",
)

sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from molecular_docking.m00_preparation.molecule_parser import run_molecule_parser

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
        description="Parse molecules from any format -> normalized table + SDF",
    )
    # Modes
    parser.add_argument("--config", "-c", type=str, default=None,
                        help="Module config YAML (03_configs/00a_molecule_parser.yaml)")
    parser.add_argument("--campaign", type=str, default=None,
                        help="Campaign config YAML")

    # Direct mode overrides
    parser.add_argument("--input", "-i", type=str, default=None,
                        help="Input file or directory (overrides campaign)")
    parser.add_argument("--output", "-o", type=str, default=None,
                        help="Output directory (overrides campaign)")

    # Parameter overrides
    parser.add_argument("--strategy", type=str, default=None,
                        choices=["best_rmsd", "first", "all"],
                        help="Conformer resolution strategy")
    parser.add_argument("--no-properties", action="store_true",
                        help="Skip molecular property calculation")
    parser.add_argument("--no-duplicates", action="store_true",
                        help="Skip SMILES duplicate detection")
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    args = parser.parse_args()

    # =========================================================================
    # RESOLVE PARAMETERS
    # =========================================================================

    # Defaults
    input_file = None
    output_dir = None
    campaign_id = "direct"
    name_source = "header_first_word"
    name_property = "_Name"
    name_column = None
    smiles_column = None
    conformer_strategy = "best_rmsd"
    compute_properties = True
    detect_duplicates = True
    log_level = "INFO"

    # --- Campaign config ---
    if args.campaign:
        cc = load_yaml(args.campaign)
        campaign_dir = Path(args.campaign).parent
        campaign_id = cc.get("campaign_id", campaign_dir.name)

        mc = cc.get("molecules", {})
        input_rel = mc.get("input_file")
        if input_rel:
            input_file = str(campaign_dir / input_rel)

        name_source = mc.get("name_source", name_source)
        name_property = mc.get("name_property", name_property)
        name_column = mc.get("name_column", name_column)
        smiles_column = mc.get("smiles_column", smiles_column)

        output_dir = str(Path("05_results") / campaign_id / "00a_molecule_parser")

    # --- Module config ---
    if args.config:
        mc = load_yaml(args.config)
        params = mc.get("parameters", {})
        conformer_strategy = params.get("conformer_strategy", conformer_strategy)
        compute_properties = params.get("compute_properties", compute_properties)
        detect_duplicates = params.get("detect_smiles_duplicates", detect_duplicates)
        log_level = params.get("log_level", log_level)

    # --- CLI overrides (highest priority) ---
    if args.input:
        input_file = args.input
    if args.output:
        output_dir = args.output
    if args.strategy:
        conformer_strategy = args.strategy
    if args.no_properties:
        compute_properties = False
    if args.no_duplicates:
        detect_duplicates = False
    if args.log_level:
        log_level = args.log_level

    # --- Validate ---
    if not input_file:
        parser.error("Provide --campaign or --input")
    if not output_dir:
        output_dir = "05_results/00a_molecule_parser"

    # =========================================================================
    # SETUP LOGGING
    # =========================================================================

    if log_level:
        logging.getLogger().setLevel(getattr(logging, log_level.upper(), logging.INFO))

    log_path = Path(output_dir) / "00a_molecule_parser.log"
    setup_log_file(log_path, log_level)

    # =========================================================================
    # EXECUTE
    # =========================================================================

    logger.info("=" * 60)
    logger.info("  MOLECULAR_DOCKING - Module 00a: Molecule Parser")
    logger.info("=" * 60)
    logger.info(f"Campaign:     {campaign_id}")
    logger.info(f"Input:        {input_file}")
    logger.info(f"Output:       {output_dir}")
    logger.info(f"Strategy:     {conformer_strategy}")
    logger.info(f"Properties:   {compute_properties}")
    logger.info(f"Duplicates:   {detect_duplicates}")

    result = run_molecule_parser(
        input_file=input_file,
        output_dir=output_dir,
        name_source=name_source,
        name_property=name_property,
        name_column=name_column,
        smiles_column=smiles_column,
        conformer_strategy=conformer_strategy,
        compute_properties=compute_properties,
        detect_smiles_duplicates=detect_duplicates,
    )

    if result.get("error"):
        logger.error(f"Error: {result['error']}")
        return 1

    logger.info("")
    logger.info(f"{'=' * 60}")
    logger.info(f"  {result['n_unique']} unique molecules saved")
    logger.info(f"  CSV: {result.get('output_csv')}")
    logger.info(f"{'=' * 60}")
    logger.info(f"Next: python 02_scripts/00c_ligand_preparation.py "
                f"--config 03_configs/00c_ligand_preparation.yaml "
                f"--campaign {args.campaign or '<campaign_config.yaml>'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
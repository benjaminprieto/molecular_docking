#!/usr/bin/env python3
"""
02a GNINA Preparation - CLI (Module 02a)
==========================================
Prepares inputs for GNINA docking: binding box + gnina_inputs.json.

GNINA accepts PDB + SDF/mol2 directly — no format conversion needed.
Much simpler than Vina (no PDBQT/Meeko dependency).

Hardcoded upstream paths:
    - 05_results/{campaign_id}/00b_receptor_preparation/rec_noH.pdb
    - 05_results/{campaign_id}/00c_ionization_profiling/{run_id}/structures/pH{xx}/
    - 05_results/{campaign_id}/00a_molecule_parser/unique_molecules.csv
    - Output: 05_results/{campaign_id}/02a_gnina_preparation/

Project: molecular_docking
Module: 02a (GNINA engine)
Version: 1.0
"""

import argparse
import yaml
import logging
import sys
from pathlib import Path
from datetime import datetime

logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)-8s | %(message)s')

sys.path.insert(0, str(Path(__file__).parent.parent / '01_src'))

from molecular_docking.m02_gnina.gnina_preparation import (
    run_gnina_preparation,
    find_reference_ligand,
)
from molecular_docking.m02_gnina.gnina_runner import check_gnina_available, find_gnina

logger = logging.getLogger(__name__)


def load_yaml(path):
    with open(path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


def setup_log_file(log_path, log_level="INFO"):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(str(log_path), mode='w', encoding='utf-8')
    fh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    fh.setFormatter(logging.Formatter('%(asctime)s | %(levelname)-8s | %(name)s | %(message)s'))
    logging.getLogger().addHandler(fh)
    logger.info(f"Log file: {log_path}")


def _find_latest_run(ionization_dir):
    run_dirs = [d for d in ionization_dir.iterdir()
                if d.is_dir() and d.name not in ("_work",) and not d.name.startswith(".")]
    if not run_dirs:
        return ionization_dir
    run_dirs.sort(key=lambda d: d.name, reverse=True)
    return run_dirs[0]


def _find_ph_folder(structures_dir, docking_ph):
    ph_label = f"pH{int(round(docking_ph * 10))}"
    ph_dir = structures_dir / ph_label
    if ph_dir.exists():
        return ph_dir
    available = [d.name for d in structures_dir.iterdir()
                 if d.is_dir() and d.name.startswith("pH")]
    if available and len(available) == 1:
        logger.info(f"  Using only available folder: {available[0]}")
        return structures_dir / available[0]
    return ph_dir


def main():
    parser = argparse.ArgumentParser(
        description='02a GNINA Preparation — binding box + inputs JSON (Module 02a)',
    )
    parser.add_argument('--config', '-c', type=str, help='Module YAML config')
    parser.add_argument('--campaign', type=str, help='Campaign config YAML')
    parser.add_argument('--receptor', type=str, default=None)
    parser.add_argument('--sdf-dir', type=str, default=None)
    parser.add_argument('--molecules', type=str, default=None)
    parser.add_argument('--reference', type=str, default=None)
    parser.add_argument('--output', '-o', type=str, default=None)
    parser.add_argument('--padding', type=float, default=None)
    parser.add_argument('--log-level', type=str, default=None)
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('--check', action='store_true', help='Check GNINA installation')

    args = parser.parse_args()

    if args.check:
        ok, msg = check_gnina_available()
        print(f"{'✓' if ok else '✗'} {msg}")
        return 0 if ok else 1

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Defaults
    receptor_pdb = None
    sdf_dir = None
    molecules_csv = None
    reference_ligand = None
    output_dir = None
    campaign_id = "direct"
    binding_box_padding = 6.0
    use_autobox = False
    autobox_ligand = None
    clean_protein = True
    include_reference = True
    ligand_format = "sdf"
    name_column = "Name"
    smiles_column = "SMILES_mol2"
    log_level = "INFO"

    if args.campaign:
        cc = load_yaml(args.campaign)
        campaign_dir = Path(args.campaign).parent
        campaign_id = cc.get("campaign_id", campaign_dir.name)

        receptor_pdb = str(Path("05_results") / campaign_id / "00b_receptor_preparation" / "rec_noH.pdb")

        # SDF dir from 00c (timestamped output)
        docking_ph = cc.get("docking_ph", 7.2)
        ionization_dir = Path("05_results") / campaign_id / "00c_ionization_profiling"
        sdf_dir = None
        if ionization_dir.exists():
            latest_run = _find_latest_run(ionization_dir)
            structures = latest_run / "structures"
            if structures.exists():
                sdf_dir = str(_find_ph_folder(structures, docking_ph))
        if not sdf_dir:
            sdf_dir = str(ionization_dir / "sdf")

        molecules_csv = str(Path("05_results") / campaign_id / "00a_molecule_parser" / "unique_molecules.csv")
        reference_ligand = find_reference_ligand(cc, campaign_dir)
        output_dir = str(Path("05_results") / campaign_id / "02a_gnina_preparation")

    if args.config:
        mc = load_yaml(args.config)
        params = mc.get("parameters", {})
        binding_box_padding = params.get("binding_box_padding", binding_box_padding)
        use_autobox = params.get("use_autobox", use_autobox)
        autobox_ligand = params.get("autobox_ligand", autobox_ligand)
        clean_protein = params.get("clean_protein", clean_protein)
        include_reference = params.get("include_reference", include_reference)
        ligand_format = params.get("ligand_format", ligand_format)
        name_column = params.get("name_column", name_column)
        log_level = params.get("log_level", log_level)

    # CLI overrides
    if args.receptor: receptor_pdb = args.receptor
    if args.sdf_dir: sdf_dir = args.sdf_dir
    if args.molecules: molecules_csv = args.molecules
    if args.reference: reference_ligand = args.reference
    if args.output: output_dir = args.output
    if args.padding: binding_box_padding = args.padding
    if args.log_level: log_level = args.log_level

    if not receptor_pdb or not sdf_dir:
        parser.error("Provide --campaign or --receptor + --sdf-dir")
    if not output_dir:
        output_dir = "05_results/02a_gnina_preparation"

    if not Path(receptor_pdb).exists():
        logger.error(f"Receptor not found: {receptor_pdb}. Run 00b first.")
        return 1
    if not Path(sdf_dir).exists():
        logger.error(f"SDF directory not found: {sdf_dir}. Run 00c first.")
        return 1

    sdf_files = sorted(Path(sdf_dir).glob("*.sdf"))
    molecule_names = [f.stem for f in sdf_files]
    if not molecule_names:
        logger.error(f"No SDF files in {sdf_dir}")
        return 1

    log_path = Path(output_dir) / '02a_gnina_preparation.log'
    setup_log_file(log_path, log_level)

    try:
        logger.info("=" * 60)
        logger.info("M02a: GNINA PREPARATION (v1.0)")
        logger.info("=" * 60)
        logger.info(f"Campaign:         {campaign_id}")
        logger.info(f"Receptor PDB:     {receptor_pdb}")
        logger.info(f"SDF directory:    {sdf_dir}")
        logger.info(f"Molecules:        {len(molecule_names)}")
        logger.info(f"Reference ligand: {reference_ligand or 'None'}")
        logger.info(f"Box padding:      {binding_box_padding}Å")
        logger.info(f"Autobox:          {use_autobox}")
        logger.info("-" * 60)

        result = run_gnina_preparation(
            receptor_pdb=receptor_pdb,
            sdf_dir=sdf_dir,
            output_dir=output_dir,
            molecule_names=molecule_names,
            reference_ligand=reference_ligand,
            binding_box_padding=binding_box_padding,
            autobox_ligand=autobox_ligand,
            use_autobox=use_autobox,
            clean_protein=clean_protein,
            include_reference=include_reference,
            ligand_format=ligand_format,
            molecules_csv=molecules_csv,
            name_column=name_column,
            smiles_column=smiles_column,
        )

        if not result.get("success"):
            logger.error(f"Failed: {result.get('error')}")
            return 1

        logger.info("")
        logger.info("OUTPUT FILES")
        logger.info(f"  Protein PDB:    {result['protein_pdb']}")
        logger.info(f"  GNINA inputs:   {result['inputs_json']}")
        logger.info(f"  Log:            {log_path}")
        logger.info("=" * 60)
        logger.info("M02a COMPLETE — Next: 02b_gnina_runner.py")
        logger.info("=" * 60)
        return 0

    except Exception as e:
        logger.error(f"Failed: {e}")
        if args.verbose:
            import traceback; traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())

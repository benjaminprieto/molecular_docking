#!/usr/bin/env python3
"""
02a Vina Preparation - CLI (Module 02a)
=========================================
Prepares receptor PDBQT + ligand PDBQTs + binding box for Vina docking.

Uses Meeko (mk_prepare_ligand.py) for ligand PDBQT conversion.
Receptor via mk_prepare_receptor.py → ADFRsuite → OpenBabel fallback.

Reads campaign_config.yaml for:
    - campaign_id       → path construction
    - receptor.pdb      → receptor source
    - grids.binding_site.reference_mol2 → binding box calculation
    - docking_ph        → ligand protonation (if needed)

Hardcoded upstream paths:
    - 05_results/{campaign_id}/00b_receptor_preparation/rec_noH.pdb
    - 05_results/{campaign_id}/00c_ionization_profiling/sdf/
    - 05_results/{campaign_id}/00a_molecule_parser/unique_molecules.csv
    - Output: 05_results/{campaign_id}/02a_vina_preparation/

Preconditions:
    - Module 00b must have been run (rec_noH.pdb)
    - Module 00c must have been run (protonated SDFs)
    - Meeko installed: pip install meeko

Usage:
    python 02_scripts/02a_vina_preparation.py --config 03_configs/02a_vina_preparation.yaml --campaign 04_data/campaigns/example_campaign/campaign_config.yaml

Project: molecular_docking
Module: 02a
Version: 1.0
"""

import argparse
import yaml
import logging
import sys
from pathlib import Path
from datetime import datetime

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)-8s | %(message)s'
)

sys.path.insert(0, str(Path(__file__).parent.parent / '01_src'))

from molecular_docking.m02_vina.vina_preparation import (
    run_vina_preparation,
    find_reference_ligand,
    check_meeko,
    check_openbabel,
)

logger = logging.getLogger(__name__)


def load_yaml(path: str) -> dict:
    with open(path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


def setup_log_file(log_path: Path, log_level: str = "INFO"):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(str(log_path), mode='w', encoding='utf-8')
    fh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    fh.setFormatter(logging.Formatter(
        '%(asctime)s | %(levelname)-8s | %(name)s | %(message)s'
    ))
    logging.getLogger().addHandler(fh)
    logger.info(f"Log file: {log_path}")


def main():
    parser = argparse.ArgumentParser(
        description='02a Vina Preparation — receptor/ligand PDBQT + binding box (Module 02a)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Standard campaign mode
  python 02_scripts/02a_vina_preparation.py --config 03_configs/02a_vina_preparation.yaml --campaign 04_data/campaigns/example_campaign/campaign_config.yaml

  # Direct mode
  python 02_scripts/02a_vina_preparation.py --receptor rec_noH.pdb --sdf-dir sdf/ --molecules unique_molecules.csv --reference ref.mol2 --output results/02a_vina_preparation/

  # Check dependencies
  python 02_scripts/02a_vina_preparation.py --check
        '''
    )

    parser.add_argument('--config', '-c', type=str, help='Module YAML config')
    parser.add_argument('--campaign', type=str, help='Campaign config YAML')

    # Direct mode
    parser.add_argument('--receptor', type=str, default=None, help='Receptor PDB (direct mode)')
    parser.add_argument('--sdf-dir', type=str, default=None, help='SDF directory from 00c')
    parser.add_argument('--molecules', type=str, default=None, help='Molecules CSV from 00a')
    parser.add_argument('--reference', type=str, default=None, help='Reference ligand for box')

    # Overrides
    parser.add_argument('--output', '-o', type=str, default=None)
    parser.add_argument('--padding', type=float, default=None, help='Binding box padding (Å)')
    parser.add_argument('--log-level', type=str, default=None)
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('--check', action='store_true', help='Check dependencies and exit')

    args = parser.parse_args()

    # =========================================================================
    # CHECK MODE
    # =========================================================================
    if args.check:
        ok_meeko, msg_meeko = check_meeko()
        ok_obabel, msg_obabel = check_openbabel()
        print(f"{'✓' if ok_meeko else '✗'} {msg_meeko}")
        print(f"{'✓' if ok_obabel else '✗'} {msg_obabel}")
        return 0 if ok_meeko else 1

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.log_level:
        logging.getLogger().setLevel(getattr(logging, args.log_level.upper(), logging.INFO))

    # =========================================================================
    # DEFAULTS
    # =========================================================================
    receptor_pdb = None
    sdf_dir = None
    molecules_csv = None
    reference_ligand = None
    output_dir = None
    campaign_id = "direct"
    binding_box_padding = 6.0
    remove_water = True
    remove_hetatm = True
    add_hydrogens = False
    include_reference = True
    name_column = "Name"
    smiles_column = "SMILES_mol2"
    log_level = "INFO"

    # =========================================================================
    # RESOLVE PARAMETERS
    # =========================================================================
    if args.campaign:
        cc = load_yaml(args.campaign)
        campaign_dir = Path(args.campaign).parent
        campaign_id = cc.get("campaign_id", campaign_dir.name)

        # Upstream paths
        receptor_pdb = str(
            Path("05_results") / campaign_id / "00b_receptor_preparation" / "rec_noH.pdb"
        )
        sdf_dir = str(
            Path("05_results") / campaign_id / "00c_ionization_profiling" / "sdf"
        )
        molecules_csv = str(
            Path("05_results") / campaign_id / "00a_molecule_parser" / "unique_molecules.csv"
        )

        # Reference ligand from campaign config
        reference_ligand = find_reference_ligand(cc, campaign_dir)
        if reference_ligand:
            logger.info(f"Reference ligand: {reference_ligand}")

        output_dir = str(Path("05_results") / campaign_id / "02a_vina_preparation")

    # Module config
    if args.config:
        mc = load_yaml(args.config)
        params = mc.get("parameters", {})
        binding_box_padding = params.get("binding_box_padding", binding_box_padding)
        remove_water = params.get("remove_water", remove_water)
        remove_hetatm = params.get("remove_hetatm", remove_hetatm)
        add_hydrogens = params.get("add_hydrogens", add_hydrogens)
        include_reference = params.get("include_reference", include_reference)
        name_column = params.get("name_column", name_column)
        smiles_column = params.get("smiles_column", smiles_column)
        log_level = params.get("log_level", log_level)

        subdir = mc.get("outputs", {}).get("subdir", "02a_vina_preparation")
        if args.campaign and not args.output:
            output_dir = str(Path("05_results") / campaign_id / subdir)

    # CLI overrides
    if args.receptor:
        receptor_pdb = args.receptor
    if args.sdf_dir:
        sdf_dir = args.sdf_dir
    if args.molecules:
        molecules_csv = args.molecules
    if args.reference:
        reference_ligand = args.reference
    if args.output:
        output_dir = args.output
    if args.padding:
        binding_box_padding = args.padding
    if args.log_level:
        log_level = args.log_level

    # Validate
    if not receptor_pdb or not sdf_dir:
        parser.error("Provide --campaign or --receptor + --sdf-dir")
    if not output_dir:
        output_dir = "05_results/02a_vina_preparation"

    if not Path(receptor_pdb).exists():
        logger.error(f"Receptor not found: {receptor_pdb}")
        logger.error("Run module 00b first.")
        return 1

    if not Path(sdf_dir).exists():
        logger.error(f"SDF directory not found: {sdf_dir}")
        logger.error("Run module 00c first.")
        return 1

    # Resolve molecule names from SDF directory
    sdf_files = sorted(Path(sdf_dir).glob("*.sdf"))
    molecule_names = [f.stem for f in sdf_files]
    if not molecule_names:
        logger.error(f"No SDF files found in {sdf_dir}")
        return 1

    # =========================================================================
    # SETUP LOG FILE
    # =========================================================================
    log_path = Path(output_dir) / '02a_vina_preparation.log'
    setup_log_file(log_path, log_level)

    # =========================================================================
    # RUN
    # =========================================================================
    try:
        logger.info("=" * 60)
        logger.info("M02a: VINA PREPARATION (v1.0)")
        logger.info("=" * 60)
        logger.info(f"Timestamp:        {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info(f"Campaign:         {campaign_id}")
        logger.info(f"Receptor PDB:     {receptor_pdb}")
        logger.info(f"SDF directory:    {sdf_dir}")
        logger.info(f"Molecules:        {len(molecule_names)}")
        logger.info(f"Reference ligand: {reference_ligand or 'None'}")
        logger.info(f"Box padding:      {binding_box_padding}Å")
        logger.info(f"Output:           {output_dir}")
        logger.info("-" * 60)

        result = run_vina_preparation(
            receptor_pdb=receptor_pdb,
            sdf_dir=sdf_dir,
            output_dir=output_dir,
            molecule_names=molecule_names,
            reference_ligand=reference_ligand,
            binding_box_padding=binding_box_padding,
            remove_water=remove_water,
            remove_hetatm=remove_hetatm,
            add_hydrogens=add_hydrogens,
            include_reference=include_reference,
            molecules_csv=molecules_csv,
            name_column=name_column,
            smiles_column=smiles_column,
        )

        if not result.get("success"):
            logger.error(f"Preparation failed: {result.get('error')}")
            return 1

        logger.info("")
        logger.info("-" * 60)
        logger.info("OUTPUT FILES")
        logger.info("-" * 60)
        logger.info(f"  Receptor PDBQT: {result['receptor_pdbqt']}")
        logger.info(f"  Ligands PDBQT:  {result['pdbqt_dir']}")
        logger.info(f"  Vina inputs:    {result['inputs_json']}")
        logger.info(f"  Log:            {log_path}")
        logger.info("=" * 60)
        logger.info("M02a COMPLETE")
        logger.info("Next: 02b_vina_runner.py")
        logger.info("=" * 60)

        return 0

    except Exception as e:
        logger.error(f"Failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
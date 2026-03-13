#!/usr/bin/env python3
"""
01b DOCK6 Run - CLI
=====================
Ejecuta DOCK6 para cada ligando preparado contra los grids.

Reads campaign_config.yaml:
    - grids.grid_dir             -> directorio con grids
    - grids.spheres_file         -> spheres_ligand.sph
    - grids.energy_grid          -> ligand.nrg
    - grids.bump_grid            -> ligand.bmp

Reads module YAML (03_configs/01b_dock6_run.yaml):
    - search_method              -> flex | rigid
    - max_orientations           -> orientaciones a probar
    - num_scored_conformers      -> poses a escribir
    - minimize                   -> simplex minimization
    - simplex_max_iterations     -> iteraciones de minimizacion
    - timeout_per_molecule       -> timeout en segundos

Grid routing (v2.0):
    1. Busca grids en 05_results/{campaign_id}/01a_grid_generation/
    2. Si no existen, busca en campaign_config.grids.grid_dir

Ligand routing (v2.0):
    1. Busca mol2 en 05_results/{campaign_id}/00d_antechamber/mol2/
    2. Si no existe, busca en 00c_ligand_preparation/mol2/

Output: 05_results/{campaign_id}/01b_dock6_run/

Usage (PyCharm Pro terminal — single line):
    python 02_scripts/01b_dock6_run.py --config 03_configs/01b_dock6_run.yaml --campaign 04_data/campaigns/example_campaign/campaign_config.yaml

    # Dock single molecule:
    python 02_scripts/01b_dock6_run.py --config 03_configs/01b_dock6_run.yaml --campaign 04_data/campaigns/example_campaign/campaign_config.yaml --name MolPort-047-542-199

    # Dry run (generate inputs only):
    python 02_scripts/01b_dock6_run.py --config 03_configs/01b_dock6_run.yaml --campaign 04_data/campaigns/example_campaign/campaign_config.yaml --dry-run

PyCharm Run Configuration:
    Script:     02_scripts/01b_dock6_run.py
    Parameters: --config 03_configs/01b_dock6_run.yaml --campaign 04_data/campaigns/example_campaign/campaign_config.yaml
    Working dir: (project root)

Project: molecular_docking
Module: 01b
Version: 2.0 — routing fix + symlinks (2026-03-13)
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

from molecular_docking.m01_docking.dock6_runner import run_dock6_batch, resolve_grid_prefix

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
        description="Run DOCK6 flexible/rigid docking for all prepared ligands",
    )
    # Modes
    parser.add_argument("--config", "-c", type=str, required=True,
                        help="Module config YAML (03_configs/01b_dock6_run.yaml)")
    parser.add_argument("--campaign", type=str, required=True,
                        help="Campaign config YAML")

    # Overrides
    parser.add_argument("--output", "-o", type=str, default=None,
                        help="Output directory")
    parser.add_argument("--name", type=str, default=None,
                        help="Dock only this molecule (by name)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Generate input files only, don't execute dock6")
    parser.add_argument("--timeout", type=int, default=None,
                        help="Timeout per molecule (seconds, overrides YAML)")
    parser.add_argument("--method", type=str, default=None,
                        choices=["flex", "rigid"],
                        help="Search method (overrides YAML)")
    parser.add_argument("--orientations", type=int, default=None,
                        help="Max orientations (overrides YAML)")
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    args = parser.parse_args()

    # =========================================================================
    # RESOLVE PARAMETERS
    # =========================================================================

    cc = load_yaml(args.campaign)
    campaign_dir = Path(args.campaign).parent
    campaign_id = cc.get("campaign_id", campaign_dir.name)

    # --- Resolve grids: check 01a output first, then campaign dir ---
    gc = cc.get("grids", {})
    grid_dir_01a = Path("05_results") / campaign_id / "01a_grid_generation"

    if (grid_dir_01a / gc.get("spheres_file", "spheres_ligand.sph")).exists():
        grid_path = grid_dir_01a
        logger.info(f"Using grids from 01a: {grid_path}")
    else:
        grid_dir = gc.get("grid_dir", "grids/")
        grid_path = Path(grid_dir) if Path(grid_dir).is_absolute() else campaign_dir / grid_dir
        if (grid_path / gc.get("spheres_file", "spheres_ligand.sph")).exists():
            logger.info(f"Using grids from campaign: {grid_path}")
        else:
            logger.warning(f"Grids not found in 01a ({grid_dir_01a}) or campaign ({grid_path})")

    spheres_file = str(grid_path / gc.get("spheres_file", "spheres_ligand.sph"))
    grid_prefix = resolve_grid_prefix(
        str(grid_path), gc.get("energy_grid", "ligand.nrg"),
    )

    # --- Ligand directory: check 00d first, then 00c fallback ---
    ligand_dir_00d = Path("05_results") / campaign_id / "00d_antechamber" / "mol2"
    ligand_dir_00c = Path("05_results") / campaign_id / "00c_ligand_preparation" / "mol2"

    if ligand_dir_00d.exists():
        ligand_dir = str(ligand_dir_00d)
    elif ligand_dir_00c.exists():
        ligand_dir = str(ligand_dir_00c)
    else:
        ligand_dir = str(ligand_dir_00d)  # Will fail with clear error below

    output_dir = args.output or str(
        Path("05_results") / campaign_id / "01b_dock6_run"
    )

    # --- Module config ---
    mc = load_yaml(args.config)
    params = mc.get("parameters", {})

    search_method = params.get("search_method", "flex")
    max_orientations = params.get("max_orientations", 1000)
    num_scored_conformers = params.get("num_scored_conformers", 20)
    minimize = params.get("minimize", True)
    simplex_max_iterations = params.get("simplex_max_iterations", 500)
    timeout_per_molecule = params.get("timeout_per_molecule", 600)
    log_level = params.get("log_level", "INFO")

    # Extra DOCK6 params (passed through to template)
    extra_params = {}
    for key in ["min_anchor_size", "pruning_max_orients", "pruning_clustering_cutoff",
                "pruning_conformer_score_cutoff", "simplex_max_cycles",
                "simplex_score_converge", "simplex_cycle_converge",
                "simplex_trans_step", "simplex_rot_step", "simplex_tors_step",
                "write_orientations"]:
        if key in params:
            extra_params[key] = params[key]

    # --- CLI overrides ---
    if args.method:
        search_method = args.method
    if args.orientations:
        max_orientations = args.orientations
    if args.timeout:
        timeout_per_molecule = args.timeout
    if args.log_level:
        log_level = args.log_level

    molecule_filter = [args.name] if args.name else None

    # =========================================================================
    # VALIDATE
    # =========================================================================

    if not Path(ligand_dir).exists():
        logger.error(f"Ligand directory not found: {ligand_dir}")
        logger.error(f"  Tried: {ligand_dir_00d}")
        logger.error(f"  Tried: {ligand_dir_00c}")
        logger.error("Run module 00d (antechamber) first.")
        return 1

    # =========================================================================
    # SETUP LOGGING
    # =========================================================================

    if log_level:
        logging.getLogger().setLevel(getattr(logging, log_level.upper(), logging.INFO))

    log_path = Path(output_dir) / "01b_dock6_run.log"
    setup_log_file(log_path, log_level)

    # =========================================================================
    # EXECUTE
    # =========================================================================

    logger.info("=" * 60)
    logger.info("  MOLECULAR_DOCKING - Module 01b: DOCK6 Run")
    logger.info("=" * 60)
    logger.info(f"Campaign:      {campaign_id}")
    logger.info(f"Ligands:       {ligand_dir}")
    logger.info(f"Grid prefix:   {grid_prefix}")
    logger.info(f"Spheres:       {spheres_file}")
    logger.info(f"Method:        {search_method}")
    logger.info(f"Orientations:  {max_orientations}")
    logger.info(f"Minimize:      {minimize} (max iter: {simplex_max_iterations})")
    logger.info(f"Poses:         {num_scored_conformers}")
    logger.info(f"Timeout:       {timeout_per_molecule}s per molecule")
    if molecule_filter:
        logger.info(f"Filter:        {molecule_filter}")
    if args.dry_run:
        logger.info("*** DRY RUN ***")

    result = run_dock6_batch(
        ligand_mol2_dir=ligand_dir,
        spheres_file=spheres_file,
        grid_prefix=grid_prefix,
        output_dir=output_dir,
        search_method=search_method,
        max_orientations=max_orientations,
        num_scored_conformers=num_scored_conformers,
        minimize=minimize,
        simplex_max_iterations=simplex_max_iterations,
        timeout_per_molecule=timeout_per_molecule,
        molecule_filter=molecule_filter,
        dry_run=args.dry_run,
        **extra_params,
    )

    if result.get("error"):
        logger.error(f"Error: {result['error']}")
        return 1

    logger.info("")
    logger.info(f"{'=' * 60}")
    if args.dry_run:
        logger.info(f"  DRY RUN: {result['n_dry_run']} input files generated")
    else:
        logger.info(f"  {result['n_ok']}/{result['n_total']} dockings completed "
                     f"({result['total_runtime_sec']:.0f}s)")
    logger.info(f"{'=' * 60}")
    logger.info(f"Next: python 02_scripts/02a_score_collection.py "
                f"--config 03_configs/02a_score_collection.yaml "
                f"--campaign {args.campaign}")

    return 0 if result["n_failed"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
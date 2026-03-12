#!/usr/bin/env python3
"""
01a Grid Generation - CLI
===========================
Genera grids DOCK6. OPCIONAL: skips si grids pre-existentes son validos.

Binding site methods (from campaign_config.yaml):
  A. reference_ligand  -> sphere_selector uses ligand mol2/pdb
  B. residues          -> centroid of specified residues
  C. coordinates       -> explicit (x, y, z) center
  D. (skip)            -> grids already exist in grid_dir

Usage:
    python 02_scripts/01a_grid_generation.py --config 03_configs/01a_grid_generation.yaml --campaign 04_data/campaigns/my_campaign/campaign_config.yaml

    # Force regeneration even if grids exist:
    python 02_scripts/01a_grid_generation.py \\
        --config 03_configs/01a_grid_generation.yaml \\
        --campaign 04_data/campaigns/my_campaign/campaign_config.yaml \\
        --force

Project: molecular_docking
Module: 01a
Version: 1.0
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

from molecular_docking.m01_docking.grid_generation import (
    run_grid_generation,
    validate_existing_grids,
    check_dock6_tools,
)

logger = logging.getLogger(__name__)


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser(description="Generate DOCK6 grids")
    parser.add_argument("--config", type=str, default=None,
                        help="Module config YAML")
    parser.add_argument("--campaign", type=str, default=None,
                        help="Campaign config YAML")
    parser.add_argument("--output", type=str, default=None)
    parser.add_argument("--force", action="store_true",
                        help="Regenerate even if grids exist")
    args = parser.parse_args()

    logger.info("=" * 60)
    logger.info("  MOLECULAR_DOCKING - Module 01a: Grid Generation")
    logger.info("=" * 60)

    if not args.campaign:
        parser.error("--campaign is required")

    # --- Load campaign config ---
    cc = load_yaml(args.campaign)
    campaign_dir = Path(args.campaign).parent
    campaign_id = cc.get("campaign_id", campaign_dir.name)

    # --- Check for pre-existing grids (Option D: skip) ---
    gc = cc.get("grids", {})
    grid_dir = gc.get("grid_dir", "grids/")
    grid_path = Path(grid_dir) if Path(grid_dir).is_absolute() else campaign_dir / grid_dir

    if not args.force and not gc.get("generate", False):
        if validate_existing_grids(
            str(grid_path),
            gc.get("spheres_file", "selected_spheres.sph"),
            gc.get("energy_grid", "grid.nrg"),
            gc.get("bump_grid", "grid.bmp"),
        ):
            logger.info(f"Valid grids found at: {grid_path}")
            logger.info("SKIPPED (use --force to regenerate, or set grids.generate: true)")
            return 0
        elif not gc.get("generate", False):
            logger.error(f"No valid grids at: {grid_path}")
            logger.error("Either provide grids or set grids.generate: true")
            return 1

    # --- Check DOCK6 tools ---
    tools = check_dock6_tools()
    missing = [t for t, ok in tools.items() if not ok]
    if missing:
        logger.error(f"Missing DOCK6 tools: {missing}")
        logger.error("Ensure dms, sphgen, sphere_selector, showbox, grid are in PATH")
        return 1

    # --- Resolve receptor paths ---
    rec_config = cc.get("receptor", {})

    # rec_noH.pdb: prefer trimmed from 00e, fallback to 00b, then campaign
    receptor_prep_dir = Path("05_results") / campaign_id / "00b_receptor_preparation"
    binding_site_dir = Path("05_results") / campaign_id / "00e_binding_site"
    rec_noH_site = binding_site_dir / "rec_noH_site.pdb"

    if rec_noH_site.exists():
        rec_noH = rec_noH_site
        logger.info(f"Using trimmed PDB from 00e: {rec_noH}")
    elif (receptor_prep_dir / "rec_noH.pdb").exists():
        rec_noH = receptor_prep_dir / "rec_noH.pdb"
    else:
        rec_noH = campaign_dir / rec_config.get("pdb", "receptor/receptor.pdb")
        logger.info(f"Using campaign receptor as noH PDB: {rec_noH}")

    rec_mol2 = receptor_prep_dir / "rec_charged.mol2"

    if not rec_mol2.exists():
        # Check prepared_mol2
        prepared = rec_config.get("prepared_mol2")
        if prepared:
            rec_mol2 = Path(prepared) if Path(prepared).is_absolute() else campaign_dir / prepared
        if not rec_mol2.exists():
            logger.error(f"Receptor mol2 not found: {rec_mol2}")
            logger.error("Run 00b (receptor preparation) first, or set receptor.prepared_mol2")
            return 1

    logger.info(f"Receptor (noH):  {rec_noH}")
    logger.info(f"Receptor (mol2): {rec_mol2}")

    # --- Resolve binding site ---
    bs = gc.get("binding_site", {})
    method = bs.get("method", "reference_ligand")
    reference_mol2 = None
    residues = None
    center = None
    site_radius = bs.get("radius", 10.0)
    chain = rec_config.get("chain")

    if method == "reference_ligand":
        ref_path = bs.get("reference_mol2")
        if not ref_path:
            logger.error("binding_site.reference_mol2 not set")
            return 1
        reference_mol2 = str(Path(ref_path) if Path(ref_path).is_absolute()
                             else campaign_dir / ref_path)
        if not Path(reference_mol2).exists():
            logger.error(f"Reference ligand not found: {reference_mol2}")
            return 1

    elif method == "residues":
        residues = bs.get("residues")
        if not residues:
            logger.error("binding_site.residues not set")
            return 1

    elif method == "coordinates":
        center = bs.get("center")
        if not center or len(center) != 3:
            logger.error(f"binding_site.center must be [x, y, z], got: {center}")
            return 1

    else:
        logger.error(f"Unknown binding_site.method: {method}")
        return 1

    logger.info(f"Binding site:    {method}")

    # --- Load module config for algorithmic params ---
    params = {}
    if args.config:
        mc = load_yaml(args.config)
        params = mc.get("parameters", {})

    output_dir = args.output or str(Path("05_results") / campaign_id / "01a_grid_generation")

    # --- Run grid generation ---
    result = run_grid_generation(
        receptor_noH_pdb=str(rec_noH),
        receptor_charged_mol2=str(rec_mol2),
        output_dir=output_dir,
        binding_site_method=method,
        reference_mol2=reference_mol2,
        residues=residues,
        center=center,
        receptor_pdb_for_residues=str(rec_noH),
        chain=chain,
        radius=site_radius,
        probe_radius=params.get("probe_radius", 1.4),
        box_margin=params.get("box_margin", 10.0),
        grid_spacing=params.get("grid_spacing", 0.3),
        energy_cutoff_distance=params.get("energy_cutoff_distance", 9999.0),
        attractive_exponent=params.get("attractive_exponent", 6),
        repulsive_exponent=params.get("repulsive_exponent", 12),
        dielectric_factor=params.get("dielectric_factor", 4),
        bump_overlap=params.get("bump_overlap", 0.75),
    )

    if result.get("success"):
        logger.info("Grid generation successful!")
        return 0
    else:
        logger.error(f"Grid generation failed: {result.get('error')}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
#!/usr/bin/env python3
"""
04b DOCK6 Footprint Analysis - CLI
======================================
Per-residue vdW + ES energy decomposition via DOCK6 footprint scoring.

Two-phase operation:
    Phase 1 (--rescore): Re-score existing poses with footprint_similarity_score_primary
                         This is fast (seconds per molecule), no re-docking.
    Phase 2 (--analyze): Parse per-residue footprint data from re-scored mol2

By default, both phases run sequentially. Use --rescore-only or --analyze-only
to run a single phase.

Input:
    01c_dock6_run/{name}/{name}_scored.mol2
    01d best_poses: reference mol2 (e.g., UDX)
    00b receptor: rec_charged.mol2

Output:
    05_results/{campaign}/05_dock6/04b_footprint_rescore/   (Phase 1)
    05_results/{campaign}/05_dock6/04b_footprint_analysis/  (Phase 2)

Project: molecular_docking
Module: 04b (DOCK6 analysis)
Version: 1.0 (2026-03-22)
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from molecular_docking.m04_dock6_analysis.footprint_rescoring import run_footprint_rescoring
from molecular_docking.m04_dock6_analysis.footprint_analysis import run_footprint_analysis

logger = logging.getLogger(__name__)


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def setup_log_file(log_path, log_level="INFO"):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(str(log_path), mode="w", encoding="utf-8")
    fh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    fh.setFormatter(logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s"))
    logging.getLogger().addHandler(fh)


def main():
    parser = argparse.ArgumentParser(
        description="04b DOCK6 Footprint Analysis — re-score + per-residue energy",
    )
    parser.add_argument("--config", "-c", type=str, help="Module config YAML")
    parser.add_argument("--campaign", type=str, help="Campaign config YAML")

    # Phase control
    parser.add_argument("--rescore-only", action="store_true",
                        help="Only run footprint re-scoring (Phase 1)")
    parser.add_argument("--analyze-only", action="store_true",
                        help="Only run footprint analysis (Phase 2)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Generate re-scoring inputs without executing")

    # Direct mode overrides
    parser.add_argument("--docking-dir", type=str, default=None,
                        help="Path to 01c_dock6_run output")
    parser.add_argument("--reference", type=str, default=None,
                        help="Reference mol2 (e.g., UDX best pose)")
    parser.add_argument("--receptor", type=str, default=None,
                        help="Receptor mol2 (rec_charged.mol2)")
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--name", type=str, default=None,
                        help="Process single molecule (for testing)")
    parser.add_argument("--log-level", type=str, default=None)

    args = parser.parse_args()

    # Defaults
    docking_dir = None
    reference_mol2 = None
    receptor_mol2 = None
    rescore_dir = None
    analysis_dir = None
    campaign_id = "direct"
    pharmacophore_threshold = 0.8
    energy_cutoff = -0.5
    timeout = 120
    log_level = "INFO"

    # --- Campaign config ---
    if args.campaign:
        cc = load_yaml(args.campaign)
        campaign_dir = Path(args.campaign).parent
        campaign_id = cc.get("campaign_id", campaign_dir.name)

        base = Path("05_results") / campaign_id
        docking_dir = str(base / "01c_dock6_run")
        receptor_mol2 = str(base / "00b_receptor_preparation" / "rec_charged.mol2")

        # Reference: best pose of reference ligand from 01d
        ref_config = cc.get("grids", {}).get("binding_site", {})
        ref_name = ref_config.get("reference_name", "UDX")
        # Look for reference in best_poses from 01d
        ref_best = base / "01d_score_collection" / "best_poses"
        if ref_best.exists():
            # Find reference by name pattern
            ref_candidates = list(ref_best.glob(f"*{ref_name}*_scored.mol2"))
            if not ref_candidates:
                ref_candidates = list(ref_best.glob(f"*{ref_name}*.mol2"))
            if ref_candidates:
                reference_mol2 = str(ref_candidates[0])

        rescore_dir = str(base / "05_dock6" / "04b_footprint_rescore")
        analysis_dir = str(base / "05_dock6" / "04b_footprint_analysis")

    # --- Module config ---
    if args.config:
        mc = load_yaml(args.config)
        params = mc.get("parameters", {})
        pharmacophore_threshold = params.get("pharmacophore_threshold", pharmacophore_threshold)
        energy_cutoff = params.get("energy_cutoff", energy_cutoff)
        timeout = params.get("timeout", timeout)
        log_level = params.get("log_level", log_level)
        if params.get("reference_mol2"):
            reference_mol2 = params["reference_mol2"]

    # --- CLI overrides ---
    if args.docking_dir:
        docking_dir = args.docking_dir
    if args.reference:
        reference_mol2 = args.reference
    if args.receptor:
        receptor_mol2 = args.receptor
    if args.output:
        rescore_dir = str(Path(args.output) / "04b_footprint_rescore")
        analysis_dir = str(Path(args.output) / "04b_footprint_analysis")
    if args.log_level:
        log_level = args.log_level

    # --- Validation ---
    if not docking_dir or not Path(docking_dir).exists():
        logger.error(f"Docking dir not found: {docking_dir}")
        logger.error("Run module 01c first, or pass --docking-dir")
        return 1

    if not args.analyze_only:
        if not reference_mol2 or not Path(reference_mol2).exists():
            logger.error(f"Reference mol2 not found: {reference_mol2}")
            logger.error("Pass --reference or set reference in campaign config")
            return 1
        if not receptor_mol2 or not Path(receptor_mol2).exists():
            logger.error(f"Receptor mol2 not found: {receptor_mol2}")
            return 1

    # --- Logging ---
    if log_level:
        logging.getLogger().setLevel(getattr(logging, log_level.upper(), logging.INFO))
    log_base = Path(analysis_dir or rescore_dir or "05_results/05_dock6/04b_footprint_analysis")
    setup_log_file(log_base / "04b_footprint.log", log_level)

    logger.info("=" * 60)
    logger.info("  MOLECULAR_DOCKING — Module 04b: DOCK6 Footprint Analysis")
    logger.info("=" * 60)
    logger.info(f"  Campaign:  {campaign_id}")
    logger.info(f"  Reference: {Path(reference_mol2).name if reference_mol2 else 'N/A'}")

    # =========================================================================
    # PHASE 1: Footprint Re-scoring
    # =========================================================================
    if not args.analyze_only:
        logger.info("")
        logger.info("  Phase 1: Footprint Re-scoring")
        logger.info("  " + "-" * 50)

        rescore_result = run_footprint_rescoring(
            docking_dir=docking_dir,
            output_dir=rescore_dir,
            reference_mol2=reference_mol2,
            receptor_mol2=receptor_mol2,
            timeout=timeout,
            dry_run=args.dry_run,
        )

        if not rescore_result.get("success"):
            logger.error(f"Re-scoring failed: {rescore_result.get('error')}")
            return 1

        if args.rescore_only or args.dry_run:
            logger.info("\nRe-scoring complete. Run with --analyze-only for Phase 2.")
            return 0

    # =========================================================================
    # PHASE 2: Footprint Analysis
    # =========================================================================
    if not args.rescore_only:
        logger.info("")
        logger.info("  Phase 2: Footprint Analysis")
        logger.info("  " + "-" * 50)

        # Use rescore output as input for analysis
        footprint_input = rescore_dir
        if args.analyze_only and not Path(rescore_dir).exists():
            # Try default location
            logger.warning(f"  Rescore dir not found: {rescore_dir}")
            logger.warning("  Run Phase 1 first (--rescore-only)")
            return 1

        analysis_result = run_footprint_analysis(
            footprint_dir=footprint_input,
            output_dir=analysis_dir,
            pharmacophore_threshold=pharmacophore_threshold,
            energy_cutoff=energy_cutoff,
        )

        if not analysis_result.get("success"):
            logger.error(f"Analysis failed: {analysis_result.get('error')}")
            return 1

    logger.info(f"\nNext: python 02_scripts/04c_binding_modes.py "
                f"--config 03_configs/04c_binding_modes.yaml "
                f"--campaign {args.campaign or '<campaign_config.yaml>'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

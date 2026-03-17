#!/usr/bin/env python3
"""
02b GNINA Runner - CLI (Module 02b)
=====================================
Runs GNINA docking (Vina + CNN scoring) for each prepared ligand.

Hardcoded upstream paths:
    - 05_results/{campaign_id}/02a_gnina_preparation/gnina_inputs.json
    - Output: 05_results/{campaign_id}/02b_gnina_run/

Project: molecular_docking
Module: 02b (GNINA engine)
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

from molecular_docking.m02_gnina.gnina_runner import (
    run_gnina_docking, find_gnina, check_gnina_available,
)

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


def main():
    parser = argparse.ArgumentParser(
        description='02b GNINA Runner — Vina + CNN docking (Module 02b)',
    )
    parser.add_argument('--config', '-c', type=str, help='Module YAML config')
    parser.add_argument('--campaign', type=str, help='Campaign config YAML')
    parser.add_argument('--inputs', type=str, default=None, help='GNINA inputs JSON (direct mode)')
    parser.add_argument('--output', '-o', type=str, default=None)
    parser.add_argument('--mode', type=str, choices=['dock', 'rescore'], default=None)
    parser.add_argument('--cnn-scoring', type=str, choices=['none', 'rescore', 'refinement', 'all'], default=None)
    parser.add_argument('--exhaustiveness', type=int, default=None)
    parser.add_argument('--num-modes', type=int, default=None)
    parser.add_argument('--seed', type=int, default=None)
    parser.add_argument('--n-workers', type=int, default=None)
    parser.add_argument('--gnina-path', type=str, default=None)
    parser.add_argument('--log-level', type=str, default=None)
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('--check', action='store_true')

    args = parser.parse_args()

    if args.check:
        gnina = find_gnina(args.gnina_path)
        if gnina:
            ok, msg = check_gnina_available(gnina)
            print(f"✓ {msg}")
            return 0
        print("✗ GNINA not found")
        return 1

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Defaults
    mode = "dock"
    cnn_scoring = "rescore"
    exhaustiveness = 8
    num_modes = 9
    seed = 42
    n_workers = 2
    timeout = 600
    gnina_path = None
    ligand_format = "auto"
    use_autobox = False
    autobox_ligand = None
    autobox_add = 6.0
    log_level = "INFO"

    inputs_json = None
    campaign_id = "direct"

    if args.campaign:
        cc = load_yaml(args.campaign)
        campaign_id = cc.get("campaign_id", Path(args.campaign).parent.name)
        inputs_json = str(Path("05_results") / campaign_id / "02a_gnina_preparation" / "gnina_inputs.json")
        output_dir = str(Path("05_results") / campaign_id / "02b_gnina_run")
    elif args.inputs:
        inputs_json = args.inputs
        output_dir = args.output or "05_results/02b_gnina_run"
    else:
        parser.error("Requires --campaign or --inputs gnina_inputs.json")

    if args.config:
        mc = load_yaml(args.config)
        params = mc.get("parameters", {})
        mode = params.get("mode", mode)
        cnn_scoring = params.get("cnn_scoring", cnn_scoring)
        exhaustiveness = params.get("exhaustiveness", exhaustiveness)
        num_modes = params.get("num_modes", num_modes)
        seed = params.get("seed", seed)
        n_workers = params.get("n_workers", n_workers)
        timeout = params.get("timeout", timeout)
        gnina_path = params.get("gnina_path", gnina_path)
        ligand_format = params.get("ligand_format", ligand_format)
        use_autobox = params.get("use_autobox", use_autobox)
        autobox_ligand = params.get("autobox_ligand", autobox_ligand)
        autobox_add = params.get("autobox_add", autobox_add)
        log_level = params.get("log_level", log_level)
        subdir = mc.get("outputs", {}).get("subdir", "02b_gnina_run")
        if args.campaign and not args.output:
            output_dir = str(Path("05_results") / campaign_id / subdir)

    # CLI overrides
    if args.output: output_dir = args.output
    if args.mode: mode = args.mode
    if args.cnn_scoring: cnn_scoring = args.cnn_scoring
    if args.exhaustiveness: exhaustiveness = args.exhaustiveness
    if args.num_modes: num_modes = args.num_modes
    if args.seed: seed = args.seed
    if args.n_workers: n_workers = args.n_workers
    if args.gnina_path: gnina_path = args.gnina_path
    if args.log_level: log_level = args.log_level

    if not inputs_json or not Path(inputs_json).exists():
        logger.error(f"GNINA inputs not found: {inputs_json}. Run 02a first.")
        return 1

    gnina = find_gnina(gnina_path)
    if not gnina:
        logger.error("GNINA not found! Install or specify --gnina-path")
        return 1

    log_path = Path(output_dir) / '02b_gnina_run.log'
    setup_log_file(log_path, log_level)

    try:
        logger.info("=" * 60)
        logger.info("M02b: GNINA DOCKING (v1.0)")
        logger.info("=" * 60)
        logger.info(f"Campaign:         {campaign_id}")
        logger.info(f"Inputs JSON:      {inputs_json}")
        logger.info(f"Mode:             {mode}")
        logger.info(f"CNN scoring:      {cnn_scoring}")
        logger.info(f"Exhaustiveness:   {exhaustiveness}")
        logger.info(f"Num modes:        {num_modes}")
        logger.info(f"Seed:             {seed}")
        logger.info(f"Workers:          {n_workers}")
        logger.info(f"Timeout:          {timeout}s")
        logger.info(f"Autobox:          {use_autobox}")
        logger.info(f"GNINA:            {gnina}")
        logger.info("-" * 60)

        result = run_gnina_docking(
            inputs_json=inputs_json, output_dir=output_dir,
            mode=mode, cnn_scoring=cnn_scoring,
            exhaustiveness=exhaustiveness, num_modes=num_modes,
            seed=seed, n_workers=n_workers, timeout=timeout,
            gnina_path=gnina, ligand_format=ligand_format,
            use_autobox=use_autobox, autobox_ligand=autobox_ligand,
            autobox_add=autobox_add,
        )

        if not result.get("success"):
            logger.error(f"Failed: {result.get('error')}")
            return 1

        logger.info("")
        logger.info("OUTPUT FILES")
        logger.info(f"  CSV:     {result['results_csv']}")
        logger.info(f"  Poses:   {result['poses_dir']}")
        logger.info("=" * 60)
        logger.info("M02b COMPLETE — Next: 02c_gnina_score_collector.py")
        logger.info("=" * 60)

        if result['n_success'] == 0:
            logger.error("No molecules docked! Check GNINA output.")
            return 2
        return 0

    except Exception as e:
        logger.error(f"Failed: {e}")
        if args.verbose:
            import traceback; traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())

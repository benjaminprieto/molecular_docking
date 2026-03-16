#!/usr/bin/env python3
"""
02b Vina Runner - CLI (Module 02b)
====================================
Runs AutoDock Vina or Vina-GPU docking for each prepared ligand.

Supports dual engine: Vina 1.2.x (CPU) and Vina-GPU (GPU-accelerated),
selectable via YAML config or --engine CLI flag.

Reads campaign_config.yaml to determine:
    - campaign_id       → path construction

Hardcoded upstream paths:
    - 05_results/{campaign_id}/02a_vina_preparation/vina_inputs.json
    - Output: 05_results/{campaign_id}/02b_vina_run/

Preconditions:
    - Module 02a must have been run (vina_inputs.json)
    - Vina or Vina-GPU installed

Usage:
    python 02_scripts/02b_vina_runner.py --config 03_configs/02b_vina_runner.yaml --campaign 04_data/campaigns/example_campaign/campaign_config.yaml

Project: molecular_docking
Module: 02b
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

from molecular_docking.m02_vina.vina_runner import (
    run_vina_docking,
    find_vina,
    check_vina_available,
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
        description='02b Vina Runner — AutoDock Vina / Vina-GPU docking (Module 02b)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Standard campaign mode (Vina CPU)
  python 02_scripts/02b_vina_runner.py --config 03_configs/02b_vina_runner.yaml --campaign 04_data/campaigns/example_campaign/campaign_config.yaml

  # Vina-GPU
  python 02_scripts/02b_vina_runner.py --config 03_configs/02b_vina_runner.yaml --campaign 04_data/campaigns/example_campaign/campaign_config.yaml --engine vina-gpu

  # Direct mode
  python 02_scripts/02b_vina_runner.py --inputs vina_inputs.json --output results/02b_vina_run/

  # Check Vina installation
  python 02_scripts/02b_vina_runner.py --check
  python 02_scripts/02b_vina_runner.py --check --engine vina-gpu
        '''
    )

    parser.add_argument('--config', '-c', type=str, help='Module YAML config')
    parser.add_argument('--campaign', type=str, help='Campaign config YAML')
    parser.add_argument('--inputs', type=str, default=None,
                        help='Vina inputs JSON from 02a (direct mode)')

    # Overrides
    parser.add_argument('--output', '-o', type=str, default=None)
    parser.add_argument('--engine', type=str, choices=['vina', 'vina-gpu'], default=None,
                        help='Docking engine override')
    parser.add_argument('--exhaustiveness', type=int, default=None)
    parser.add_argument('--num-modes', type=int, default=None)
    parser.add_argument('--seed', type=int, default=None)
    parser.add_argument('--n-workers', type=int, default=None)
    parser.add_argument('--vina-path', type=str, default=None)
    parser.add_argument('--log-level', type=str, default=None)
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('--check', action='store_true', help='Check Vina installation')

    args = parser.parse_args()

    # =========================================================================
    # CHECK MODE
    # =========================================================================
    if args.check:
        engine = args.engine or "vina"
        vina = find_vina(args.vina_path, engine)
        if vina:
            ok, msg = check_vina_available(vina, engine)
            print(f"✓ {msg}")
            return 0
        else:
            print(f"✗ {engine} not found")
            print(f"  Tried common locations. Use --vina-path to specify.")
            return 1

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.log_level:
        logging.getLogger().setLevel(getattr(logging, args.log_level.upper(), logging.INFO))

    # =========================================================================
    # DEFAULTS
    # =========================================================================
    engine = "vina"
    exhaustiveness = 8
    num_modes = 9
    energy_range = 3.0
    seed = 42
    n_workers = 4
    cpu_per_worker = 1
    timeout = 600
    vina_path = None
    gpu_batch_size = 0
    use_autobox = False
    autobox_ligand = None
    autobox_add = 6.0
    log_level = "INFO"

    # =========================================================================
    # RESOLVE PARAMETERS
    # =========================================================================
    inputs_json = None
    campaign_id = "direct"

    if args.campaign:
        cc = load_yaml(args.campaign)
        campaign_dir = Path(args.campaign).parent
        campaign_id = cc.get("campaign_id", campaign_dir.name)

        # Hardcoded upstream path
        inputs_json = str(
            Path("05_results") / campaign_id
            / "02a_vina_preparation" / "vina_inputs.json"
        )
        output_dir = str(Path("05_results") / campaign_id / "02b_vina_run")

    elif args.inputs:
        inputs_json = args.inputs
        output_dir = args.output or "05_results/02b_vina_run"

    else:
        parser.error("Requires --campaign or --inputs vina_inputs.json")

    # Module config
    if args.config:
        mc = load_yaml(args.config)
        params = mc.get("parameters", {})
        engine = params.get("engine", engine)
        exhaustiveness = params.get("exhaustiveness", exhaustiveness)
        num_modes = params.get("num_modes", num_modes)
        energy_range = params.get("energy_range", energy_range)
        seed = params.get("seed", seed)
        n_workers = params.get("n_workers", n_workers)
        cpu_per_worker = params.get("cpu_per_worker", cpu_per_worker)
        timeout = params.get("timeout", timeout)
        vina_path = params.get("vina_path", vina_path)
        gpu_batch_size = params.get("gpu_batch_size", gpu_batch_size)
        use_autobox = params.get("use_autobox", use_autobox)
        autobox_ligand = params.get("autobox_ligand", autobox_ligand)
        autobox_add = params.get("autobox_add", autobox_add)
        log_level = params.get("log_level", log_level)

        subdir = mc.get("outputs", {}).get("subdir", "02b_vina_run")
        if args.campaign and not args.output:
            output_dir = str(Path("05_results") / campaign_id / subdir)

    # CLI overrides
    if args.output:
        output_dir = args.output
    if args.engine is not None:
        engine = args.engine
    if args.exhaustiveness is not None:
        exhaustiveness = args.exhaustiveness
    if args.num_modes is not None:
        num_modes = args.num_modes
    if args.seed is not None:
        seed = args.seed
    if args.n_workers is not None:
        n_workers = args.n_workers
    if args.vina_path is not None:
        vina_path = args.vina_path
    if args.log_level:
        log_level = args.log_level

    # Resolve autobox ligand from campaign config
    if use_autobox and not autobox_ligand and args.campaign:
        cc = load_yaml(args.campaign)
        campaign_dir = Path(args.campaign).parent
        bs = cc.get("grids", {}).get("binding_site", {})
        ref = bs.get("reference_mol2")
        if ref:
            ref_path = Path(ref) if Path(ref).is_absolute() else campaign_dir / ref
            if ref_path.exists():
                autobox_ligand = str(ref_path)

    # Validate
    if not inputs_json or not Path(inputs_json).exists():
        logger.error(f"Vina inputs not found: {inputs_json}")
        logger.error("Run module 02a first.")
        return 1

    # Check engine availability
    vina = find_vina(vina_path, engine)
    if not vina:
        logger.error(f"{engine} not found!")
        logger.error("Install Vina or specify path with --vina-path")
        return 1

    # =========================================================================
    # SETUP LOG FILE
    # =========================================================================
    log_path = Path(output_dir) / '02b_vina_run.log'
    setup_log_file(log_path, log_level)

    # =========================================================================
    # RUN
    # =========================================================================
    try:
        logger.info("=" * 60)
        logger.info("M02b: VINA DOCKING (v1.0)")
        logger.info("=" * 60)
        logger.info(f"Timestamp:        {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info(f"Campaign:         {campaign_id}")
        logger.info(f"Inputs JSON:      {inputs_json}")
        logger.info(f"Output dir:       {output_dir}")
        logger.info(f"Engine:           {engine}")
        logger.info(f"Exhaustiveness:   {exhaustiveness}")
        logger.info(f"Num modes:        {num_modes}")
        logger.info(f"Energy range:     {energy_range}")
        logger.info(f"Seed:             {seed}")
        logger.info(f"Workers:          {n_workers}")
        logger.info(f"CPU/worker:       {cpu_per_worker}")
        logger.info(f"Timeout:          {timeout}s")
        logger.info(f"Autobox:          {use_autobox}")
        if use_autobox:
            logger.info(f"Autobox ligand:   {autobox_ligand or 'from JSON'}")
            logger.info(f"Autobox padding:  {autobox_add}Å")
        logger.info(f"Vina path:        {vina}")
        logger.info(f"Log file:         {log_path}")
        logger.info("-" * 60)

        result = run_vina_docking(
            inputs_json=inputs_json,
            output_dir=output_dir,
            exhaustiveness=exhaustiveness,
            num_modes=num_modes,
            energy_range=energy_range,
            seed=seed,
            n_workers=n_workers,
            cpu_per_worker=cpu_per_worker,
            timeout=timeout,
            vina_path=vina,
            engine=engine,
            gpu_batch_size=gpu_batch_size,
            use_autobox=use_autobox,
            autobox_ligand=autobox_ligand,
            autobox_add=autobox_add,
        )

        if not result.get("success"):
            logger.error(f"Docking failed: {result.get('error')}")
            return 1

        # Summary
        logger.info("")
        logger.info("-" * 60)
        logger.info("OUTPUT FILES")
        logger.info("-" * 60)
        logger.info(f"  CSV:     {result['results_csv']}")
        logger.info(f"  JSON:    {result['results_json']}")
        logger.info(f"  Summary: {result['summary_txt']}")
        logger.info(f"  Poses:   {result['poses_dir']}")
        logger.info(f"  Log:     {log_path}")
        logger.info("=" * 60)
        logger.info("M02b COMPLETE")
        logger.info("Next: 02c_vina_score_collector.py")
        logger.info("=" * 60)

        if result['n_success'] == 0:
            logger.error("No molecules were docked! Check Vina output and ligand files.")
            return 2

        return 0

    except Exception as e:
        logger.error(f"Failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
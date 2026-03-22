#!/usr/bin/env python3
"""
05d Binding Site Hotspots - CLI
=================================
Cross-molecule analysis: where in the binding site do fragments
from different molecules converge?

Reads:  05a + 05b + 05c (optional)
Output: 05_results/{campaign_id}/05d_binding_site_hotspots/

Usage:
    python 02_scripts/05d_binding_site_hotspots.py \
        --config 03_configs/05d_binding_site_hotspots.yaml \
        --campaign 04_data/campaigns/UDX_pharmit_pH63/campaign_config.yaml

Project: molecular_docking
Module: 05d
Version: 2.0
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)-8s | %(message)s")

sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from molecular_docking.m05_gnina_analysis.binding_site_hotspots import run_binding_site_hotspots

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
        description="05d Binding Site Hotspots — cross-molecule fragment convergence",
    )
    parser.add_argument("--config", "-c", type=str, required=True)
    parser.add_argument("--campaign", type=str, required=True)
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--receptor", type=str, default=None,
                        help="Receptor PDB for ChimeraX script")
    parser.add_argument("--eps", type=float, default=None,
                        help="Hotspot DBSCAN eps (Angstrom)")
    parser.add_argument("--min-molecules", type=int, default=None,
                        help="Min molecules to form a hotspot")
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    args = parser.parse_args()

    cc = load_yaml(args.campaign)
    mc = load_yaml(args.config)
    params = mc.get("parameters", {})

    campaign_dir = Path(args.campaign).parent
    campaign_id = cc.get("campaign_id", campaign_dir.name)
    results_base = Path("05_results") / campaign_id / "m05_gnina_analysis"

    output_subdir = mc.get("outputs", {}).get("subdir", "05d_binding_site_hotspots")
    output_dir = Path(args.output) if args.output else results_base / output_subdir

    log_level = args.log_level or params.get("log_level", "INFO")
    setup_log_file(output_dir / "05d_binding_site_hotspots.log", log_level)

    hotspot_eps = args.eps or params.get("hotspot_eps", 2.0)
    min_molecules = args.min_molecules or params.get("min_molecules", 3)
    min_dominant_fraction = params.get("min_dominant_fraction", 0.15)
    include_non_dominant = params.get("include_non_dominant", False)

    # Receptor path
    receptor_path = args.receptor
    if not receptor_path:
        rp = params.get("receptor_path")
        if rp:
            receptor_path = rp
        else:
            rec = cc.get("receptor", {})
            rp = rec.get("pdb") or rec.get("path")
            if rp:
                rpath = Path(rp) if Path(rp).is_absolute() else campaign_dir / rp
                if rpath.exists():
                    receptor_path = str(rpath)

    parsed_dir = str(results_base / "05a_parse_and_fragment")
    cluster_dir = str(results_base / "05b_fragment_clustering")
    score_dir = str(results_base / "05c_score_decomposition")

    logger.info(f"Campaign:       {campaign_id}")
    logger.info(f"Hotspot eps:    {hotspot_eps} A")
    logger.info(f"Min molecules:  {min_molecules}")
    logger.info(f"Receptor:       {receptor_path or 'not specified'}")
    logger.info(f"Output:         {output_dir}")

    result = run_binding_site_hotspots(
        parsed_dir=parsed_dir,
        cluster_dir=cluster_dir,
        output_dir=str(output_dir),
        score_dir=score_dir,
        receptor_path=receptor_path,
        hotspot_eps=hotspot_eps,
        min_molecules=min_molecules,
        min_dominant_fraction=min_dominant_fraction,
        include_non_dominant=include_non_dominant,
    )

    if not result.get("success"):
        logger.error(f"Failed: {result.get('error')}")
        sys.exit(1)


if __name__ == "__main__":
    main()

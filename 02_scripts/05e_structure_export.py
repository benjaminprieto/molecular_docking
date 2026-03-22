#!/usr/bin/env python3
"""
05e Structure Export - CLI
============================
Export fragment cluster representatives as mol2 + ChimeraX scripts,
informed by hotspot analysis from 05d.

Reads:  05a + 05b + 05c (optional) + 05d (optional)
Output: 05_results/{campaign_id}/05e_structure_export/

Usage:
    python 02_scripts/05e_structure_export.py \
        --config 03_configs/05e_structure_export.yaml \
        --campaign 04_data/campaigns/UDX_pharmit_pH63/campaign_config.yaml

Project: molecular_docking
Module: 05e
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

from molecular_docking.m05_gnina_analysis.structure_export import run_structure_export

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
        description="05e Structure Export — mol2 representatives + ChimeraX scripts",
    )
    parser.add_argument("--config", "-c", type=str, required=True)
    parser.add_argument("--campaign", type=str, required=True)
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--receptor", type=str, default=None)
    parser.add_argument("--name", type=str, default=None)
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    args = parser.parse_args()

    cc = load_yaml(args.campaign)
    mc = load_yaml(args.config)
    params = mc.get("parameters", {})

    campaign_dir = Path(args.campaign).parent
    campaign_id = cc.get("campaign_id", campaign_dir.name)
    results_base = Path("05_results") / campaign_id / "m05_gnina_analysis"

    output_subdir = mc.get("outputs", {}).get("subdir", "05e_structure_export")
    output_dir = Path(args.output) if args.output else results_base / output_subdir

    log_level = args.log_level or params.get("log_level", "INFO")
    setup_log_file(output_dir / "05e_structure_export.log", log_level)

    max_clusters = params.get("max_clusters_per_fragment", 3)

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
    hotspot_dir = str(results_base / "05d_binding_site_hotspots")
    molecule_names = [args.name] if args.name else None

    logger.info(f"Campaign:  {campaign_id}")
    logger.info(f"Receptor:  {receptor_path or 'not specified'}")
    logger.info(f"Output:    {output_dir}")

    result = run_structure_export(
        parsed_dir=parsed_dir,
        cluster_dir=cluster_dir,
        output_dir=str(output_dir),
        score_dir=score_dir,
        hotspot_dir=hotspot_dir,
        receptor_path=receptor_path,
        max_clusters_per_fragment=max_clusters,
        molecule_names=molecule_names,
    )

    if not result.get("success"):
        logger.error(f"Failed: {result.get('error')}")
        sys.exit(1)


if __name__ == "__main__":
    main()

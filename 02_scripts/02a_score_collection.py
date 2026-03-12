#!/usr/bin/env python3
"""02a Score Collection - CLI"""
import argparse, logging, sys, yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))
from molecular_docking.m02_collection.score_collector import run_score_collection

logger = logging.getLogger(__name__)

def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f: return yaml.safe_load(f)

def main():
    parser = argparse.ArgumentParser(description="Collect scores → dock2profile Excel")
    parser.add_argument("--config", type=str, required=True)
    parser.add_argument("--campaign", type=str, required=True)
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args()

    cc = load_yaml(args.campaign)
    campaign_id = cc.get("campaign_id", Path(args.campaign).parent.name)
    oc = cc.get("output", {})

    docking_dir = str(Path("05_results") / campaign_id / "01b_dock6_run")
    molecules_csv = str(Path("05_results") / campaign_id / "00a_molecule_parser" / "unique_molecules.csv")
    output_dir = args.output or str(Path("05_results") / campaign_id / "02a_score_collection")
    params = load_yaml(args.config).get("parameters", {})

    logger.info("=" * 60)
    logger.info("  MOLECULAR_DOCKING - Module 02a: Score Collection")
    logger.info("=" * 60)

    result = run_score_collection(
        docking_dir=docking_dir, molecules_csv=molecules_csv,
        output_dir=output_dir,
        score_key=params.get("score_key", "Grid_Score"),
        max_molecules=params.get("max_molecules", 500),
        scores_filename=oc.get("scores_filename", "01_top_500_molecules.xlsx"),
        mol2_dirname=oc.get("mol2_dirname", "docked_molecules"),
        source_label=cc.get("metadata", {}).get("source"),
    )
    return 0

if __name__ == "__main__":
    sys.exit(main())

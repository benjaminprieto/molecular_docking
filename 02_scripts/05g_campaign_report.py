#!/usr/bin/env python3
"""
05g Campaign Report - CLI
===========================
Generates consolidated HTML report with figures, tables, and text summary.

Usage:
    python 02_scripts/05g_campaign_report.py -c 03_configs/05g_campaign_report.yaml \
        --campaign 04_data/campaigns/UDX_pharmit_pH63/campaign_config.yaml

Project: molecular_docking
Module: 05g
Version: 1.0
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)-8s | %(message)s')
sys.path.insert(0, str(Path(__file__).parent.parent / '01_src'))

from molecular_docking.m05_gnina_analysis.campaign_report import run_campaign_report

logger = logging.getLogger(__name__)


def load_yaml(path):
    with open(path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser(description='05g Campaign Report')
    parser.add_argument('--config', '-c', type=str, help='Module YAML config')
    parser.add_argument('--campaign', type=str, required=True, help='Campaign config YAML')
    parser.add_argument('--output', '-o', type=str, default=None)
    parser.add_argument('--reference', type=str, default=None,
                        help='Reference mol2 for RMSD validation (e.g. crystallographic ligand)')
    parser.add_argument('--control-name', type=str, default='UDX',
                        help='Name of control molecule to compute RMSD against reference')
    parser.add_argument('--name', type=str, default=None, help='Filter single molecule')
    args = parser.parse_args()

    cc = load_yaml(args.campaign)
    campaign_id = cc.get('campaign_id', Path(args.campaign).parent.name)

    engine = 'gnina'

    # Module config
    params = {}
    if args.config:
        mc = load_yaml(args.config)
        params = mc.get('parameters', {})

    results_base = Path('05_results') / campaign_id / 'm05_gnina_analysis'
    output_dir = args.output or str(results_base / '05g_campaign_report')

    # Reference mol2 for RMSD validation
    reference_mol2 = args.reference
    if not reference_mol2:
        # Auto-detect from campaign
        campaign_dir = Path(args.campaign).parent
        auto_ref = campaign_dir / 'reference' / 'UDX.mol2'
        if auto_ref.exists():
            reference_mol2 = str(auto_ref)

    logger.info("=" * 60)
    logger.info(f"M05g: CAMPAIGN REPORT")
    logger.info("=" * 60)
    logger.info(f"Campaign:  {campaign_id}")
    logger.info(f"Engine:    {engine}")
    logger.info(f"Reference: {reference_mol2 or 'none'}")
    logger.info(f"Output:    {output_dir}")

    result = run_campaign_report(
        parsed_dir=str(results_base / '05a_parse_and_fragment'),
        cluster_dir=str(results_base / '05b_fragment_clustering'),
        score_dir=str(results_base / '05c_score_decomposition'),
        hotspot_dir=str(results_base / '05d_binding_site_hotspots'),
        contact_dir=str(results_base / '05f_contact_mapping'),
        output_dir=output_dir,
        campaign_id=campaign_id,
        engine=engine,
        reference_mol2=reference_mol2,
        control_name=args.control_name,
    )

    if not result.get('success'):
        logger.error(f"Failed: {result.get('error')}")
        return 1

    logger.info(f"\nOpen in browser: {result['report_path']}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
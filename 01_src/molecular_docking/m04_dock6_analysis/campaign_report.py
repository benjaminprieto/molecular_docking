"""
Campaign Report - Core Module (04e)
=======================================
Generates an HTML report for DOCK6 analysis answering 5 key questions:

    1. What molecules bind best?       → Grid_Score ranking (vdW + ES)
    2. Where do they bind?             → binding modes
    3. Why do they bind?               → footprint (which residues, vdW vs ES)
    4. How confident?                  → n_modes, score spread
    5. What residues matter?           → footprint consensus

Input:
    04a: molecule_ranking.csv, score_components.csv
    04b: residue_consensus.csv, pharmacophore_residues.json, molecule_footprint_summary.csv
    04c: binding_modes_summary.csv
    04d: contact_summary.csv, contact_vs_footprint.csv

Output:
    campaign_report.html        — full HTML report
    composite_ranking.csv       — combined ranking

Location: 01_src/molecular_docking/m04_dock6_analysis/campaign_report.py
Project: molecular_docking
Module: 04e (DOCK6 analysis)
Version: 1.0 (2026-03-22)
"""

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# COMPOSITE RANKING
# =============================================================================

def _build_composite_ranking(
        score_csv: str,
        modes_csv: Optional[str] = None,
        fps_summary_csv: Optional[str] = None,
        weights: Optional[Dict[str, float]] = None,
) -> pd.DataFrame:
    """
    Build composite ranking from multiple DOCK6 analysis metrics.

    Weights (default):
        grid_score:  0.40   (binding energy)
        n_modes:     0.15   (binding mode diversity)
        fps_score:   0.30   (footprint similarity to reference)
        score_std:   0.15   (reproducibility across modes)
    """
    if weights is None:
        weights = {
            "grid_score": 0.40,
            "n_modes": 0.15,
            "fps_score": 0.30,
            "score_std": 0.15,
        }

    df = pd.read_csv(score_csv)

    # Merge modes data
    if modes_csv and Path(modes_csv).exists():
        df_modes = pd.read_csv(modes_csv)
        df = df.merge(df_modes[["Name", "n_modes", "mean_pairwise_rmsd"]],
                       on="Name", how="left")
    else:
        df["n_modes"] = 1
        df["mean_pairwise_rmsd"] = 0.0

    # Merge footprint summary
    if fps_summary_csv and Path(fps_summary_csv).exists():
        df_fps = pd.read_csv(fps_summary_csv)
        merge_cols = ["Name"]
        if "fps_score" in df_fps.columns:
            merge_cols.append("fps_score")
        if "n_residues_contributing" in df_fps.columns:
            merge_cols.append("n_residues_contributing")
        df = df.merge(df_fps[merge_cols], on="Name", how="left")

    # Normalize each metric to 0-1 (higher = better)
    def normalize_min_better(series):
        """Lower is better → 1.0 = best (most negative)"""
        valid = series.dropna()
        if len(valid) == 0 or valid.max() == valid.min():
            return pd.Series(0.5, index=series.index)
        return (series.max() - series) / (series.max() - series.min())

    def normalize_max_better(series):
        """Higher is better → 1.0 = best"""
        valid = series.dropna()
        if len(valid) == 0 or valid.max() == valid.min():
            return pd.Series(0.5, index=series.index)
        return (series - series.min()) / (series.max() - series.min())

    # Grid_Score: more negative = better
    df["norm_grid_score"] = normalize_min_better(df["Grid_Score"])

    # n_modes: more = more confident
    df["norm_n_modes"] = normalize_max_better(df["n_modes"])

    # FPS score: lower = more similar to reference = better
    if "fps_score" in df.columns:
        df["norm_fps_score"] = normalize_min_better(df["fps_score"].fillna(df["fps_score"].max()))
    else:
        df["norm_fps_score"] = 0.5

    # Score std: lower = more reproducible = better
    if "std_score" in df.columns:
        df["norm_score_std"] = normalize_min_better(df["std_score"].fillna(df["std_score"].max()))
    else:
        df["norm_score_std"] = 0.5

    # Composite score
    df["composite_score"] = (
        weights["grid_score"] * df["norm_grid_score"] +
        weights["n_modes"] * df["norm_n_modes"] +
        weights["fps_score"] * df["norm_fps_score"] +
        weights["score_std"] * df["norm_score_std"]
    ).round(4)

    df.sort_values("composite_score", ascending=False, inplace=True)
    df.reset_index(drop=True, inplace=True)
    df.insert(0, "Composite_Rank", range(1, len(df) + 1))

    return df


# =============================================================================
# HTML REPORT GENERATION
# =============================================================================

def _generate_html(
        df_ranking: pd.DataFrame,
        df_residues: Optional[pd.DataFrame],
        pharma_data: Optional[Dict],
        df_modes: Optional[pd.DataFrame],
        df_contacts: Optional[pd.DataFrame],
        campaign_id: str,
        n_molecules: int,
) -> str:
    """Generate complete HTML report."""

    # Top 20 molecules table
    top20 = df_ranking.head(20)
    ranking_table = _df_to_html_table(top20, [
        "Composite_Rank", "Name", "Grid_Score", "Grid_vdw_energy",
        "Grid_es_energy", "n_modes", "composite_score",
    ])

    # Residue consensus table
    residue_table = ""
    if df_residues is not None and len(df_residues) > 0:
        top_res = df_residues.head(20)
        residue_table = _df_to_html_table(top_res, [
            "residue_id", "n_molecules_contributing", "frac_contributing",
            "mean_vdw", "mean_es", "mean_total",
        ])

    # Pharmacophore residues
    pharma_section = ""
    if pharma_data and pharma_data.get("residues"):
        pharma_items = ""
        for r in pharma_data["residues"]:
            pharma_items += (f"<li><strong>{r['residue_id']}</strong>: "
                           f"{r['frac_contributing']*100:.0f}% of molecules, "
                           f"mean energy = {r['mean_total']:.2f} kcal/mol</li>\n")
        pharma_section = f"""
        <h3>Pharmacophore Residues (contacted by ≥{pharma_data['threshold']*100:.0f}% of molecules)</h3>
        <ul>{pharma_items}</ul>
        """

    # Binding modes summary
    modes_section = ""
    if df_modes is not None and len(df_modes) > 0:
        modes_stats = (f"<p>Modes per molecule: {df_modes['n_modes'].min()}–"
                      f"{df_modes['n_modes'].max()} "
                      f"(mean {df_modes['n_modes'].mean():.1f})</p>")
        modes_table = _df_to_html_table(df_modes.head(15), [
            "Name", "n_modes", "best_Grid_Score", "score_spread",
            "mean_pairwise_rmsd",
        ])
        modes_section = f"<h3>Binding Mode Diversity</h3>\n{modes_stats}\n{modes_table}"

    # Contact summary
    contacts_section = ""
    if df_contacts is not None and len(df_contacts) > 0:
        top_contacts = df_contacts.head(15)
        contacts_table = _df_to_html_table(top_contacts, [
            "residue_id", "n_molecules", "frac_molecules",
            "mean_min_distance", "mean_n_contacts",
        ])
        contacts_section = f"<h3>Contact Residues</h3>\n{contacts_table}"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>DOCK6 Campaign Report — {campaign_id}</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif;
               max-width: 1200px; margin: 0 auto; padding: 20px; color: #333; background: #fafafa; }}
        h1 {{ color: #1a5276; border-bottom: 3px solid #2980b9; padding-bottom: 10px; }}
        h2 {{ color: #2c3e50; margin-top: 30px; border-bottom: 1px solid #bdc3c7; padding-bottom: 5px; }}
        h3 {{ color: #34495e; }}
        table {{ border-collapse: collapse; width: 100%; margin: 15px 0; font-size: 13px; }}
        th {{ background: #2c3e50; color: white; padding: 8px 12px; text-align: left; }}
        td {{ padding: 6px 12px; border-bottom: 1px solid #ecf0f1; }}
        tr:hover {{ background: #eaf2f8; }}
        tr:nth-child(even) {{ background: #f8f9fa; }}
        .summary-box {{ background: #eaf2f8; border-left: 4px solid #2980b9;
                        padding: 15px; margin: 15px 0; border-radius: 4px; }}
        .question {{ font-weight: bold; color: #2980b9; margin-top: 10px; }}
        .metric {{ display: inline-block; background: #2980b9; color: white;
                   padding: 3px 8px; border-radius: 3px; font-size: 12px; margin: 2px; }}
        ul {{ line-height: 1.8; }}
        .footer {{ margin-top: 40px; padding-top: 20px; border-top: 1px solid #bdc3c7;
                   font-size: 12px; color: #7f8c8d; }}
    </style>
</head>
<body>
    <h1>DOCK6 Campaign Report</h1>
    <div class="summary-box">
        <p><strong>Campaign:</strong> {campaign_id}</p>
        <p><strong>Molecules:</strong> {n_molecules}</p>
        <p><strong>Engine:</strong> DOCK6 (Grid-based scoring + Footprint analysis)</p>
        <p><strong>Date:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
    </div>

    <h2>1. What molecules bind best?</h2>
    <p class="question">→ Composite ranking: Grid_Score (40%) + Footprint similarity (30%)
       + Mode diversity (15%) + Score reproducibility (15%)</p>
    {ranking_table}

    <h2>2. Where do they bind?</h2>
    <p class="question">→ Each DOCK6 pose is a genuine binding mode (pre-clustered at 2.0 Å RMSD)</p>
    {modes_section}

    <h2>3. Why do they bind?</h2>
    <p class="question">→ DOCK6 footprint: per-residue vdW + ES energy decomposition</p>
    {residue_table}
    {pharma_section}

    <h2>4. How confident are we?</h2>
    <p class="question">→ More binding modes + tighter score spread = more confident prediction</p>
    <p>Molecules with multiple modes that all score well are the most reliable hits.
       A single mode with a great score could be a false positive.</p>

    <h2>5. What residues define the binding site?</h2>
    <p class="question">→ Geometric contacts + energetic footprint cross-reference</p>
    {contacts_section}

    <div class="footer">
        <p>Generated by molecular_docking v2.1 — Module 04e (DOCK6 analysis)</p>
        <p>Reference: Balius et al. J Chem Inf Model 2011, 51(8):1942-56</p>
    </div>
</body>
</html>"""

    return html


def _df_to_html_table(df: pd.DataFrame, columns: list) -> str:
    """Convert a DataFrame to an HTML table string."""
    cols = [c for c in columns if c in df.columns]
    if not cols:
        return "<p><em>No data available.</em></p>"

    header = "".join(f"<th>{c}</th>" for c in cols)
    rows = ""
    for _, row in df.iterrows():
        cells = ""
        for c in cols:
            val = row[c]
            if isinstance(val, float):
                cells += f"<td>{val:.3f}</td>"
            else:
                cells += f"<td>{val}</td>"
        rows += f"<tr>{cells}</tr>\n"

    return f"<table>\n<thead><tr>{header}</tr></thead>\n<tbody>\n{rows}</tbody>\n</table>"


# =============================================================================
# MAIN REPORT GENERATOR
# =============================================================================

def run_campaign_report(
        results_base: Union[str, Path],
        output_dir: Union[str, Path],
        campaign_id: str = "",
        composite_weights: Optional[Dict[str, float]] = None,
) -> Dict[str, Any]:
    """
    Generate DOCK6 campaign HTML report.

    Reads outputs from 04a-04d and produces a composite ranking + HTML report.

    Args:
        results_base: Base directory containing 04a-04d outputs
                      (05_results/{campaign}/05_dock6/)
        output_dir:   Output directory for report
        campaign_id:  Campaign identifier for report title
        composite_weights: Weights for composite ranking

    Returns:
        Dict with: success, output paths
    """
    results_base = Path(results_base)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("  04e DOCK6 Campaign Report v1.0")
    logger.info("=" * 60)

    # --- Locate input files ---
    score_csv = results_base / "04a_score_ranking" / "score_components.csv"
    consensus_csv = results_base / "04b_footprint_analysis" / "residue_consensus.csv"
    pharma_json = results_base / "04b_footprint_analysis" / "pharmacophore_residues.json"
    fps_summary = results_base / "04b_footprint_analysis" / "molecule_footprint_summary.csv"
    modes_csv = results_base / "04c_binding_modes" / "binding_modes_summary.csv"
    contacts_csv = results_base / "04d_contact_mapping" / "contact_summary.csv"

    if not score_csv.exists():
        return {"success": False,
                "error": f"04a score_components.csv not found: {score_csv}"}

    # --- Build composite ranking ---
    df_ranking = _build_composite_ranking(
        str(score_csv),
        modes_csv=str(modes_csv) if modes_csv.exists() else None,
        fps_summary_csv=str(fps_summary) if fps_summary.exists() else None,
        weights=composite_weights,
    )

    composite_csv = output_dir / "composite_ranking.csv"
    df_ranking.to_csv(composite_csv, index=False, encoding="utf-8")
    logger.info(f"  Saved: {composite_csv}")

    # --- Load supplementary data ---
    df_residues = pd.read_csv(consensus_csv) if consensus_csv.exists() else None
    pharma_data = None
    if pharma_json.exists():
        with open(pharma_json) as f:
            pharma_data = json.load(f)
    df_modes = pd.read_csv(modes_csv) if modes_csv.exists() else None
    df_contacts = pd.read_csv(contacts_csv) if contacts_csv.exists() else None

    # --- Generate HTML ---
    html = _generate_html(
        df_ranking=df_ranking,
        df_residues=df_residues,
        pharma_data=pharma_data,
        df_modes=df_modes,
        df_contacts=df_contacts,
        campaign_id=campaign_id,
        n_molecules=len(df_ranking),
    )

    report_html = output_dir / "campaign_report.html"
    with open(report_html, "w", encoding="utf-8") as f:
        f.write(html)
    logger.info(f"  Saved: {report_html}")

    # --- Top 5 summary ---
    logger.info("")
    logger.info("  Top 5 Composite Ranking:")
    for _, row in df_ranking.head(5).iterrows():
        logger.info(f"    {int(row['Composite_Rank'])}. {row['Name']}  "
                     f"Grid={row['Grid_Score']:.2f}  "
                     f"composite={row['composite_score']:.3f}")

    logger.info("=" * 60)

    return {
        "success": True,
        "n_molecules": len(df_ranking),
        "composite_ranking_csv": str(composite_csv),
        "campaign_report_html": str(report_html),
        "output_dir": str(output_dir),
    }

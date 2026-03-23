"""
Score Ranking - Core Module (04a)
====================================
Ranks DOCK6 docked molecules by Grid_Score and decomposes into vdW + ES.

Input:
    dock6_all_poses.csv from 01d (1 row per pose per molecule)

Output:
    molecule_ranking.csv   — 1 row per molecule, ranked by best Grid_Score
    score_components.csv   — vdW + ES breakdown, dominance flag, per-pose stats

Location: 01_src/molecular_docking/m04_dock6_analysis/score_ranking.py
Project: molecular_docking
Module: 04a (DOCK6 analysis)
Version: 1.0 (2026-03-22)
"""

import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional, Union

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def run_score_ranking(
        all_poses_csv: Union[str, Path],
        output_dir: Union[str, Path],
        score_key: str = "Grid_Score",
        vdw_key: str = "Grid_vdw_energy",
        es_key: str = "Grid_es_energy",
        top_n: int = 0,
) -> Dict[str, Any]:
    """
    Rank molecules by DOCK6 Grid_Score and decompose into vdW + ES.

    For each molecule:
      - Best pose Grid_Score (most negative)
      - vdW and ES components of best pose
      - vdW/ES dominance classification
      - Per-pose statistics: mean, std, n_poses, score range

    Args:
        all_poses_csv: Path to dock6_all_poses.csv (from 01d)
        output_dir:    Output directory for ranking files
        score_key:     Column to rank by
        vdw_key:       vdW energy column
        es_key:        ES energy column
        top_n:         Limit to top N molecules (0 = all)

    Returns:
        Dict with: success, n_molecules, output paths
    """
    all_poses_csv = Path(all_poses_csv)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not all_poses_csv.exists():
        return {"success": False, "error": f"Not found: {all_poses_csv}"}

    logger.info("=" * 60)
    logger.info("  04a DOCK6 Score Ranking v1.0")
    logger.info("=" * 60)

    # --- Load all poses ---
    df = pd.read_csv(all_poses_csv)
    logger.info(f"  Loaded: {len(df)} poses from {all_poses_csv.name}")

    required = [score_key, "Name"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        return {"success": False,
                "error": f"Missing columns: {missing}. Available: {list(df.columns)}"}

    # --- Per-molecule statistics ---
    stats = []
    for name, group in df.groupby("Name"):
        scores = group[score_key].values
        best_idx = scores.argmin()
        best_row = group.iloc[best_idx]

        best_score = float(best_row[score_key])
        best_vdw = float(best_row[vdw_key]) if vdw_key in group.columns else np.nan
        best_es = float(best_row[es_key]) if es_key in group.columns else np.nan

        # vdW/ES dominance: which component contributes more to binding?
        # More negative = stronger contribution
        if not np.isnan(best_vdw) and not np.isnan(best_es):
            vdw_frac = best_vdw / best_score if best_score != 0 else 0
            es_frac = best_es / best_score if best_score != 0 else 0
            if abs(best_vdw) > 2 * abs(best_es):
                dominance = "vdW_dominant"
            elif abs(best_es) > 2 * abs(best_vdw):
                dominance = "ES_dominant"
            else:
                dominance = "balanced"
        else:
            vdw_frac = np.nan
            es_frac = np.nan
            dominance = "unknown"

        stats.append({
            "Name": name,
            "Grid_Score": best_score,
            "Grid_vdw_energy": best_vdw,
            "Grid_es_energy": best_es,
            "vdW_fraction": round(vdw_frac, 3) if not np.isnan(vdw_frac) else np.nan,
            "ES_fraction": round(es_frac, 3) if not np.isnan(es_frac) else np.nan,
            "dominance": dominance,
            "n_poses": len(group),
            "mean_score": round(float(scores.mean()), 3),
            "std_score": round(float(scores.std()), 3) if len(scores) > 1 else 0.0,
            "worst_score": round(float(scores.max()), 3),
            "score_range": round(float(scores.max() - scores.min()), 3),
        })

    df_rank = pd.DataFrame(stats)
    df_rank.sort_values("Grid_Score", ascending=True, inplace=True)
    df_rank.reset_index(drop=True, inplace=True)
    df_rank.insert(0, "Rank", range(1, len(df_rank) + 1))

    if top_n > 0:
        df_rank = df_rank.head(top_n)

    n_molecules = len(df_rank)

    # --- Save molecule_ranking.csv ---
    ranking_csv = output_dir / "molecule_ranking.csv"
    rank_cols = ["Rank", "Name", "Grid_Score", "Grid_vdw_energy", "Grid_es_energy",
                 "n_poses", "mean_score", "std_score"]
    df_rank[rank_cols].to_csv(ranking_csv, index=False, encoding="utf-8")
    logger.info(f"  Saved: {ranking_csv}")

    # --- Save score_components.csv ---
    components_csv = output_dir / "score_components.csv"
    df_rank.to_csv(components_csv, index=False, encoding="utf-8")
    logger.info(f"  Saved: {components_csv}")

    # --- Summary ---
    n_vdw = sum(1 for s in stats if s["dominance"] == "vdW_dominant")
    n_es = sum(1 for s in stats if s["dominance"] == "ES_dominant")
    n_bal = sum(1 for s in stats if s["dominance"] == "balanced")

    logger.info("")
    logger.info(f"  Molecules ranked: {n_molecules}")
    logger.info(f"  Best: {df_rank.iloc[0]['Name']} ({df_rank.iloc[0]['Grid_Score']:.2f})")
    logger.info(f"  Worst: {df_rank.iloc[-1]['Name']} ({df_rank.iloc[-1]['Grid_Score']:.2f})")
    logger.info(f"  Dominance: {n_vdw} vdW, {n_es} ES, {n_bal} balanced")
    logger.info("=" * 60)

    return {
        "success": True,
        "n_molecules": n_molecules,
        "n_vdw_dominant": n_vdw,
        "n_es_dominant": n_es,
        "n_balanced": n_bal,
        "molecule_ranking_csv": str(ranking_csv),
        "score_components_csv": str(components_csv),
        "output_dir": str(output_dir),
    }

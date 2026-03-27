"""
Campaign Report - Core Module (04e)
=======================================
Generates a publication-quality HTML report for DOCK6 analysis.

Answers 5 key screening questions:
    1. What molecules bind best?       → Grid_Score composite ranking
    2. Where do they bind?             → Binding modes (n_modes, diversity)
    3. Why do they bind?               → Footprint (per-residue vdW + ES)
    4. How confident are we?           → n_modes, score spread, convergence
    5. What residues matter?           → Pharmacophore residues, contacts

Reads: 04a-04d outputs
Saves: campaign_report.html + figures/ directory + composite_ranking.csv

Location: 01_src/molecular_docking/m04_dock6_analysis/campaign_report.py
Project: molecular_docking
Module: 04e (DOCK6 analysis)
Version: 4.0 (2026-03-27) — composite = mean Grid_Score across all poses

Reference: Balius et al. J Chem Inf Model 2011, 51(8):1942-56
"""

import json
import logging
import base64
from io import BytesIO
from pathlib import Path
from typing import Dict, Any, Optional, Union
from datetime import datetime

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    logger.warning("matplotlib not available — no figures in report")


# =============================================================================
# FIGURE HELPERS
# =============================================================================

def _fig_to_base64(fig) -> str:
    """Convert matplotlib figure to base64 PNG for HTML embedding."""
    buf = BytesIO()
    fig.savefig(buf, format='png', dpi=120, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode('utf-8')


# =============================================================================
# FIGURE 1: MOLECULE RANKING (Grid_Score bar chart)
# =============================================================================

def fig_molecule_ranking(df: pd.DataFrame, top_n: int = 30) -> str:
    """Horizontal bar chart of top molecules by Grid_Score."""
    if not HAS_MPL or df.empty:
        return ""

    score_col = "Grid_Score" if "Grid_Score" in df.columns else df.columns[1]
    top = df.nsmallest(top_n, score_col)

    fig, ax = plt.subplots(figsize=(10, max(6, len(top) * 0.3)))
    colors = plt.cm.RdYlGn_r(np.linspace(0.2, 0.8, len(top)))
    bars = ax.barh(range(len(top)), top[score_col].values, color=colors, edgecolor='gray', linewidth=0.5)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([n[:30] for n in top["Name"].values], fontsize=7)
    ax.set_xlabel("Grid_Score (kcal/mol)")
    ax.set_title(f"Top {len(top)} Molecules by Grid_Score")
    ax.invert_yaxis()
    ax.axvline(x=0, color='gray', linestyle='-', alpha=0.3)

    for i, (val, name) in enumerate(zip(top[score_col].values, top["Name"].values)):
        ax.text(val - 0.5, i, f"{val:.1f}", va='center', ha='right', fontsize=6, color='white', fontweight='bold')

    fig.tight_layout()
    return _fig_to_base64(fig)


# =============================================================================
# FIGURE 2: vdW vs ES DECOMPOSITION
# =============================================================================

def fig_vdw_vs_es(df: pd.DataFrame) -> str:
    """Scatter plot: vdW vs ES energy per molecule."""
    if not HAS_MPL or df.empty:
        return ""

    vdw_col = next((c for c in df.columns if 'vdw' in c.lower() and 'energy' in c.lower()), None)
    es_col = next((c for c in df.columns if 'es' in c.lower() and 'energy' in c.lower()), None)

    if not vdw_col or not es_col:
        # Try Grid_vdw_energy, Grid_es_energy
        vdw_col = "Grid_vdw_energy" if "Grid_vdw_energy" in df.columns else None
        es_col = "Grid_es_energy" if "Grid_es_energy" in df.columns else None

    if not vdw_col or not es_col:
        return ""

    fig, ax = plt.subplots(figsize=(8, 6))
    score_col = "Grid_Score" if "Grid_Score" in df.columns else None

    if score_col:
        sc = ax.scatter(df[vdw_col], df[es_col], c=df[score_col], cmap='RdYlGn_r',
                        s=40, alpha=0.7, edgecolors='gray', linewidth=0.5)
        plt.colorbar(sc, ax=ax, label='Grid_Score (kcal/mol)')
    else:
        ax.scatter(df[vdw_col], df[es_col], s=40, alpha=0.7, edgecolors='gray', linewidth=0.5)

    ax.set_xlabel("vdW Energy (kcal/mol)")
    ax.set_ylabel("ES Energy (kcal/mol)")
    ax.set_title("Energy Decomposition: vdW vs Electrostatic")
    ax.axhline(y=0, color='gray', linestyle=':', alpha=0.3)
    ax.axvline(x=0, color='gray', linestyle=':', alpha=0.3)

    # Quadrant labels
    ax.text(0.02, 0.02, 'vdW-dominated\n(hydrophobic)', transform=ax.transAxes,
            fontsize=8, color='#2980b9', alpha=0.5)
    ax.text(0.02, 0.95, 'Both unfavorable', transform=ax.transAxes,
            fontsize=8, color='#e74c3c', alpha=0.5, va='top')
    ax.text(0.95, 0.02, 'ES-dominated\n(polar/charged)', transform=ax.transAxes,
            fontsize=8, color='#27ae60', alpha=0.5, ha='right')

    fig.tight_layout()
    return _fig_to_base64(fig)


# =============================================================================
# FIGURE 3: PHARMACOPHORE RESIDUES
# =============================================================================

def fig_pharmacophore_residues(df_consensus: pd.DataFrame, pharma_data: dict,
                               top_n: int = 25) -> str:
    """Bar chart of top contributing residues with pharmacophore threshold."""
    if not HAS_MPL or df_consensus is None or df_consensus.empty:
        return ""

    top = df_consensus.head(top_n).copy()
    threshold = pharma_data.get("threshold", 0.8) if pharma_data else 0.8

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, max(6, len(top) * 0.3)))

    # Left: mean energy per residue
    colors_e = ['#e74c3c' if v < -1.0 else '#f39c12' if v < -0.5 else '#bdc3c7'
                for v in top["mean_total"].values]
    ax1.barh(range(len(top)), top["mean_total"].values, color=colors_e, edgecolor='gray', linewidth=0.5)
    ax1.set_yticks(range(len(top)))
    ax1.set_yticklabels(top["residue_id"].values, fontsize=7)
    ax1.set_xlabel("Mean Total Energy (kcal/mol)")
    ax1.set_title("Top Residues by Energy")
    ax1.invert_yaxis()
    ax1.axvline(x=-0.5, color='orange', linestyle='--', alpha=0.5, label='cutoff=-0.5')
    ax1.legend(fontsize=7)

    # Right: fraction of molecules contributing
    colors_f = ['#27ae60' if v >= threshold else '#3498db' if v >= 0.5 else '#bdc3c7'
                for v in top["frac_contributing"].values]
    ax2.barh(range(len(top)), top["frac_contributing"].values, color=colors_f, edgecolor='gray', linewidth=0.5)
    ax2.set_yticks(range(len(top)))
    ax2.set_yticklabels(top["residue_id"].values, fontsize=7)
    ax2.set_xlabel("Fraction of Molecules Contributing")
    ax2.set_title("Residue Frequency (Pharmacophore)")
    ax2.invert_yaxis()
    ax2.axvline(x=threshold, color='green', linestyle='--', alpha=0.5, label=f'threshold={threshold}')
    ax2.set_xlim(0, 1.05)
    ax2.legend(fontsize=7)

    fig.tight_layout()
    return _fig_to_base64(fig)


# =============================================================================
# FIGURE 4: BINDING MODES
# =============================================================================

def fig_binding_modes(df_modes: pd.DataFrame) -> str:
    """Histogram of n_modes + scatter of score_spread vs n_modes."""
    if not HAS_MPL or df_modes is None or df_modes.empty:
        return ""

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Left: histogram of n_modes
    n_modes = df_modes["n_modes"].values if "n_modes" in df_modes.columns else []
    if len(n_modes) > 0:
        ax1.hist(n_modes, bins=range(1, max(n_modes) + 2), color='#3498db',
                 edgecolor='white', alpha=0.8, align='left')
        ax1.set_xlabel("Number of Binding Modes")
        ax1.set_ylabel("Number of Molecules")
        ax1.set_title("Binding Mode Distribution")
        ax1.axvline(x=np.mean(n_modes), color='red', linestyle='--', alpha=0.5,
                     label=f'mean={np.mean(n_modes):.1f}')
        ax1.legend(fontsize=8)

    # Right: score spread vs n_modes
    if "score_spread" in df_modes.columns and "n_modes" in df_modes.columns:
        sc = ax2.scatter(df_modes["n_modes"], df_modes["score_spread"],
                         c=df_modes.get("best_Grid_Score", df_modes["n_modes"]),
                         cmap='RdYlGn_r', s=40, alpha=0.7, edgecolors='gray', linewidth=0.5)
        ax2.set_xlabel("Number of Binding Modes")
        ax2.set_ylabel("Score Spread (kcal/mol)")
        ax2.set_title("Diversity vs Score Range")
        if "best_Grid_Score" in df_modes.columns:
            plt.colorbar(sc, ax=ax2, label='Best Grid_Score')

    fig.tight_layout()
    return _fig_to_base64(fig)


# =============================================================================
# FIGURE 5: CONTACT HEATMAP
# =============================================================================

def fig_contact_heatmap(df_contacts: pd.DataFrame, top_residues: int = 20,
                        top_molecules: int = 25) -> str:
    """Heatmap: top residues vs top molecules by contact distance or count."""
    if not HAS_MPL or df_contacts is None or df_contacts.empty:
        return ""

    # Find appropriate columns
    res_col = "residue_id" if "residue_id" in df_contacts.columns else None
    name_col = "Name" if "Name" in df_contacts.columns else None
    dist_col = next((c for c in df_contacts.columns if 'distance' in c.lower()), None)
    count_col = next((c for c in df_contacts.columns if 'n_contacts' in c.lower() or 'count' in c.lower()), None)

    if not res_col or not name_col:
        return ""

    # Top residues by frequency
    res_freq = df_contacts.groupby(res_col)[name_col].nunique().sort_values(ascending=False)
    top_res = res_freq.head(top_residues).index.tolist()

    # Top molecules
    mol_counts = df_contacts[df_contacts[res_col].isin(top_res)][name_col].value_counts()
    top_mol = mol_counts.head(top_molecules).index.tolist()

    val_col = dist_col or count_col
    if not val_col:
        return ""

    filt = df_contacts[df_contacts[res_col].isin(top_res) & df_contacts[name_col].isin(top_mol)]
    pivot = filt.pivot_table(values=val_col, index=name_col, columns=res_col, aggfunc='min')
    pivot = pivot.reindex(columns=top_res)

    if pivot.empty:
        return ""

    fig, ax = plt.subplots(figsize=(max(8, top_residues * 0.5),
                                     max(6, len(pivot) * 0.25)))
    im = ax.imshow(pivot.values, cmap='RdYlGn', aspect='auto', vmin=2.0, vmax=5.0)
    plt.colorbar(im, ax=ax, label=f'{val_col} (Å)')

    ax.set_xticks(range(len(top_res)))
    ax.set_xticklabels(top_res, rotation=45, ha='right', fontsize=7)
    ax.set_yticks(range(len(pivot)))
    ax.set_yticklabels([n[:25] for n in pivot.index], fontsize=6)
    ax.set_title(f'Contact Map: Top {top_residues} Residues × Top {len(pivot)} Molecules')

    fig.tight_layout()
    return _fig_to_base64(fig)


# =============================================================================
# FIGURE 6: FOOTPRINT vs REFERENCE
# =============================================================================

def fig_footprint_vs_reference(df_consensus: pd.DataFrame, top_n: int = 30) -> str:
    """Grouped bar chart: mean pose energy vs reference energy for top residues."""
    if not HAS_MPL or df_consensus is None or df_consensus.empty:
        return ""

    if "ref_vdw" not in df_consensus.columns or "ref_es" not in df_consensus.columns:
        return ""

    top = df_consensus.head(top_n).copy()
    top["ref_total"] = top["ref_vdw"] + top["ref_es"]

    fig, ax = plt.subplots(figsize=(14, 6))
    x = np.arange(len(top))
    w = 0.35

    ax.bar(x - w / 2, top["mean_total"].values, w, label='Campaign Mean', color='#3498db', alpha=0.8)
    ax.bar(x + w / 2, top["ref_total"].values, w, label='Reference (UDX)', color='#e67e22', alpha=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(top["residue_id"].values, rotation=45, ha='right', fontsize=7)
    ax.set_ylabel("Energy (kcal/mol)")
    ax.set_title("Per-Residue Footprint: Campaign vs Reference")
    ax.legend()
    ax.axhline(y=0, color='gray', linestyle='-', alpha=0.3)

    fig.tight_layout()
    return _fig_to_base64(fig)


# =============================================================================
# HTML TABLE HELPER
# =============================================================================

def _html_table(df: pd.DataFrame, columns: list = None, max_rows: int = 50) -> str:
    """Convert DataFrame to styled HTML table."""
    if df is None or df.empty:
        return "<p><em>No data available.</em></p>"

    cols = columns if columns else list(df.columns)
    cols = [c for c in cols if c in df.columns]
    if not cols:
        return "<p><em>No matching columns.</em></p>"

    show = df.head(max_rows)
    header = "".join(f"<th>{c}</th>" for c in cols)
    rows = ""
    for _, row in show.iterrows():
        cells = ""
        for c in cols:
            val = row[c]
            if isinstance(val, float):
                if abs(val) < 0.01 and val != 0:
                    cells += f"<td>{val:.4f}</td>"
                else:
                    cells += f"<td>{val:.2f}</td>"
            else:
                cells += f"<td>{val}</td>"
        rows += f"<tr>{cells}</tr>\n"

    more = f"<p style='font-size:11px;color:#999'>Showing {len(show)}/{len(df)} rows.</p>" if len(df) > max_rows else ""

    return f"""<table class="data-table">
<thead><tr>{header}</tr></thead>
<tbody>{rows}</tbody>
</table>{more}"""


# =============================================================================
# COMPOSITE RANKING
# =============================================================================

def _build_composite_ranking(
        score_csv: str,
        all_poses_csv: Optional[str] = None,
) -> pd.DataFrame:
    """
    Build composite ranking: mean Grid_Score across all binding modes.

    The mean captures whether a molecule binds well in ALL modes (consistent)
    or only in one outlier pose. More negative mean = better.
    """
    df = pd.read_csv(score_csv)

    # Get Grid_Score column (best pose)
    score_col = "Grid_Score" if "Grid_Score" in df.columns else df.columns[1]

    # Compute per-molecule stats from all poses
    if all_poses_csv and Path(all_poses_csv).exists():
        df_poses = pd.read_csv(all_poses_csv)
        pose_score_col = "Grid_Score" if "Grid_Score" in df_poses.columns else df_poses.columns[1]

        pose_stats = df_poses.groupby("Name").agg(
            mean_Grid_Score=(pose_score_col, "mean"),
            n_modes=(pose_score_col, "count"),
            best_Grid_Score=(pose_score_col, "min"),
            worst_Grid_Score=(pose_score_col, "max"),
        ).reset_index()
        pose_stats["spread"] = pose_stats["worst_Grid_Score"] - pose_stats["best_Grid_Score"]
        pose_stats["ratio"] = np.where(
            pose_stats["best_Grid_Score"] != 0,
            pose_stats["mean_Grid_Score"] / pose_stats["best_Grid_Score"],
            0.0,
        )

        df = df.merge(pose_stats, on="Name", how="left")
    else:
        # Fallback: use best Grid_Score as mean (no all_poses data)
        df["mean_Grid_Score"] = df[score_col]
        df["n_modes"] = 1
        df["spread"] = 0.0
        df["ratio"] = 1.0
        logger.warning("  No all_poses_csv — composite = best Grid_Score (no multi-mode analysis)")

    df = df.fillna(0)

    # Composite = mean Grid_Score (simple, interpretable, kcal/mol)
    df["composite"] = df["mean_Grid_Score"]

    # Sort by composite (most negative first = best)
    df = df.sort_values("composite", ascending=True).reset_index(drop=True)
    df["Rank"] = range(1, len(df) + 1)
    cols = ["Rank"] + [c for c in df.columns if c != "Rank"]
    df = df[cols]

    return df


# =============================================================================
# HTML REPORT GENERATOR
# =============================================================================

def _generate_html(
        df_ranking: pd.DataFrame,
        df_residues: Optional[pd.DataFrame],
        pharma_data: Optional[dict],
        df_modes: Optional[pd.DataFrame],
        df_contacts: Optional[pd.DataFrame],
        df_fps_summary: Optional[pd.DataFrame],
        figures: Dict[str, str],
        campaign_id: str,
        n_molecules: int,
        df_drug: Optional[pd.DataFrame] = None,
        df_residue_per_pose: Optional[pd.DataFrame] = None,
        df_all_poses: Optional[pd.DataFrame] = None,
) -> str:
    """Generate the complete HTML report with 3 tables + detail cards."""

    # --- Compute summary stats ---
    score_col = "Grid_Score" if "Grid_Score" in df_ranking.columns else "composite"
    best_score = df_ranking[score_col].min() if score_col in df_ranking.columns else 0
    worst_score = df_ranking[score_col].max() if score_col in df_ranking.columns else 0
    mean_score = df_ranking[score_col].mean() if score_col in df_ranking.columns else 0
    best_name = df_ranking.iloc[0]["Name"] if len(df_ranking) > 0 else "N/A"
    best_composite = df_ranking.iloc[0]["composite"] if "composite" in df_ranking.columns else 0

    n_pharmacophore = pharma_data.get("n_pharmacophore_residues", 0) if pharma_data else 0
    n_residues = len(df_residues) if df_residues is not None else 0

    n_modes_mean = df_modes["n_modes"].mean() if df_modes is not None and "n_modes" in df_modes.columns else 0
    n_modes_range = f"{df_modes['n_modes'].min()}-{df_modes['n_modes'].max()}" if df_modes is not None and "n_modes" in df_modes.columns else "N/A"

    # vdW vs ES dominance
    n_vdw_dom = 0
    n_es_dom = 0
    n_balanced = 0
    if "Grid_vdw_energy" in df_ranking.columns and "Grid_es_energy" in df_ranking.columns:
        for _, row in df_ranking.iterrows():
            vdw = abs(row.get("Grid_vdw_energy", 0))
            es = abs(row.get("Grid_es_energy", 0))
            total = vdw + es
            if total > 0:
                if vdw / total > 0.65:
                    n_vdw_dom += 1
                elif es / total > 0.65:
                    n_es_dom += 1
                else:
                    n_balanced += 1

    # Merge drug-likeness into ranking for table display
    df_display = df_ranking.copy()
    if df_drug is not None and not df_drug.empty and "Name" in df_drug.columns:
        drug_cols = [c for c in ["Name", "MW", "QED", "TPSA", "LogP"] if c in df_drug.columns]
        if len(drug_cols) > 1:
            df_display = df_display.merge(df_drug[drug_cols], on="Name", how="left")

    # --- TABLE 1: Top 20 by Grid_Score (raw) ---
    df_table1 = df_display.nsmallest(20, score_col).copy()
    df_table1 = df_table1.reset_index(drop=True)
    df_table1["GS_Rank"] = range(1, len(df_table1) + 1)
    df_table1 = df_table1[["GS_Rank"] + [c for c in df_table1.columns if c != "GS_Rank"]]
    table1_cols = ["GS_Rank", "Name", score_col]
    if "Grid_vdw_energy" in df_table1.columns:
        table1_cols.append("Grid_vdw_energy")
    if "Grid_es_energy" in df_table1.columns:
        table1_cols.append("Grid_es_energy")
    if "n_poses" in df_table1.columns:
        table1_cols.append("n_poses")
    for dc in ["MW", "QED", "TPSA", "LogP"]:
        if dc in df_table1.columns:
            table1_cols.append(dc)
    table1_html = _html_table(df_table1, table1_cols)

    # --- TABLE 2: Top 20 by Composite (mean Grid_Score) ---
    df_table2 = df_display.head(20).copy()
    table2_cols = ["Rank", "Name"]
    if "best_Grid_Score" in df_table2.columns:
        table2_cols.append("best_Grid_Score")
    table2_cols.append("composite")  # = mean_Grid_Score
    if "n_modes" in df_table2.columns:
        table2_cols.append("n_modes")
    if "spread" in df_table2.columns:
        table2_cols.append("spread")
    if "ratio" in df_table2.columns:
        table2_cols.append("ratio")
    for dc in ["MW", "QED"]:
        if dc in df_table2.columns:
            table2_cols.append(dc)
    table2_html = _html_table(df_table2, table2_cols)

    # --- TABLE 3: Top 20 receptor residues ---
    residue_table_html = ""
    if df_residues is not None and not df_residues.empty:
        df_res = df_residues.copy()
        # Compute sorting metric: n_molecules × |mean_energy|
        if "n_molecules" in df_res.columns and "mean_total" in df_res.columns:
            df_res["importance"] = df_res["n_molecules"] * df_res["mean_total"].abs()
            df_res = df_res.sort_values("importance", ascending=False).reset_index(drop=True)
        elif "frac_contributing" in df_res.columns and "mean_total" in df_res.columns:
            df_res["importance"] = df_res["frac_contributing"] * df_res["mean_total"].abs()
            df_res = df_res.sort_values("importance", ascending=False).reset_index(drop=True)
        df_res_top = df_res.head(20).copy()
        df_res_top["Res_Rank"] = range(1, len(df_res_top) + 1)
        df_res_top = df_res_top[["Res_Rank"] + [c for c in df_res_top.columns if c != "Res_Rank"]]
        res_cols = ["Res_Rank", "residue_id"]
        if "n_molecules" in df_res_top.columns:
            res_cols.append("n_molecules")
        if "mean_total" in df_res_top.columns:
            res_cols.append("mean_total")
        if "std_total" in df_res_top.columns:
            res_cols.append("std_total")
        if "mean_vdw" in df_res_top.columns:
            res_cols.append("mean_vdw")
        if "mean_es" in df_res_top.columns:
            res_cols.append("mean_es")
        if "frac_contributing" in df_res_top.columns:
            res_cols.append("frac_contributing")
        residue_table_html = _html_table(df_res_top, res_cols)

    # --- DETAIL CARDS: Top 20 composite molecules ---
    detail_cards_html = ""
    top20_composite = df_display.head(20)
    for _, mol_row in top20_composite.iterrows():
        mol_name = mol_row["Name"]
        best_gs = mol_row.get("best_Grid_Score", mol_row.get(score_col, 0))
        comp = mol_row.get("composite", 0)
        n_m = mol_row.get("n_modes", 0)
        spread_val = mol_row.get("spread", 0)
        ratio_val = mol_row.get("ratio", 0)
        mw_val = mol_row.get("MW", "N/A")
        qed_val = mol_row.get("QED", "N/A")
        tpsa_val = mol_row.get("TPSA", "N/A")

        # Header
        mw_str = f"{mw_val:.1f}" if isinstance(mw_val, (int, float)) and mw_val != "N/A" else str(mw_val)
        qed_str = f"{qed_val:.3f}" if isinstance(qed_val, (int, float)) and qed_val != "N/A" else str(qed_val)
        tpsa_str = f"{tpsa_val:.1f}" if isinstance(tpsa_val, (int, float)) and tpsa_val != "N/A" else str(tpsa_val)

        card_html = f"""
<div style="border:1px solid #bdc3c7;border-radius:8px;padding:15px;margin:15px 0;background:white;">
    <h4 style="margin:0 0 10px 0;color:#1a5276">{mol_name}</h4>
    <p style="font-size:12px;margin:5px 0">
        <b>Mean Grid_Score (composite):</b> {comp:.2f} &nbsp;|&nbsp;
        <b>Best Grid_Score:</b> {best_gs:.2f} &nbsp;|&nbsp;
        <b>n_modes:</b> {int(n_m)} &nbsp;|&nbsp;
        <b>Spread:</b> {spread_val:.1f} &nbsp;|&nbsp;
        <b>Ratio:</b> {ratio_val:.3f} &nbsp;|&nbsp;
        <b>MW:</b> {mw_str} &nbsp;|&nbsp;
        <b>QED:</b> {qed_str} &nbsp;|&nbsp;
        <b>TPSA:</b> {tpsa_str}
    </p>
"""
        # Per-residue table from residue_per_pose
        if df_residue_per_pose is not None and not df_residue_per_pose.empty:
            mol_residues = df_residue_per_pose[df_residue_per_pose["Name"] == mol_name].copy() if "Name" in df_residue_per_pose.columns else pd.DataFrame()
            if not mol_residues.empty:
                # Sort by |mean_energy| descending
                energy_col = "mean_energy" if "mean_energy" in mol_residues.columns else "mean_total" if "mean_total" in mol_residues.columns else None
                if energy_col:
                    mol_residues["abs_energy"] = mol_residues[energy_col].abs()
                    mol_residues = mol_residues.sort_values("abs_energy", ascending=False).head(10)

                    res_table_rows = ""
                    for _, rrow in mol_residues.iterrows():
                        res_id = rrow.get("residue_id", "?")
                        e_val = rrow.get(energy_col, 0)
                        vdw_val = rrow.get("vdw", rrow.get("mean_vdw", 0))
                        es_val = rrow.get("es", rrow.get("mean_es", 0))
                        cons_val = rrow.get("consistency", rrow.get("frac_consistent", 0))
                        n_fav = rrow.get("n_favorable", "")
                        n_pos = rrow.get("n_poses", n_p)

                        # Assessment
                        if e_val < -5 and cons_val > 0.8:
                            badge = '<span style="color:#f39c12">&#9733; Anchor</span>'
                        elif e_val < -2 and cons_val > 0.7:
                            badge = '<span style="color:#27ae60">&#10003; Reliable</span>'
                        elif e_val < -5 and cons_val < 0.5:
                            badge = '<span style="color:#e67e22">&#9889; Strong but variable</span>'
                        elif cons_val < 0.3:
                            badge = '<span style="color:#e74c3c">&#10007; Unreliable</span>'
                        else:
                            badge = ""

                        n_fav_str = f"{int(n_fav)}/{int(n_pos)}" if n_fav != "" and not pd.isna(n_fav) else ""
                        res_table_rows += (
                            f"<tr><td>{res_id}</td><td>{e_val:.2f}</td>"
                            f"<td>{vdw_val:.2f}</td><td>{es_val:.2f}</td>"
                            f"<td>{cons_val:.2f}</td><td>{n_fav_str}</td>"
                            f"<td>{badge}</td></tr>\n"
                        )

                    card_html += f"""
    <table class="data-table" style="font-size:11px">
    <thead><tr><th>Residue</th><th>Mean Energy</th><th>vdW</th><th>ES</th>
               <th>Consistency</th><th>n_favorable/n_poses</th><th>Assessment</th></tr></thead>
    <tbody>{res_table_rows}</tbody>
    </table>
"""
        # Per-mode Grid_Score table
        if df_all_poses is not None and not df_all_poses.empty:
            mol_poses = df_all_poses[df_all_poses["Name"] == mol_name].copy() if "Name" in df_all_poses.columns else pd.DataFrame()
            if not mol_poses.empty:
                pose_score_col = "Grid_Score" if "Grid_Score" in mol_poses.columns else mol_poses.columns[1]
                mol_poses = mol_poses.sort_values(pose_score_col, ascending=True)
                mode_rows = ""
                for pi, (_, prow) in enumerate(mol_poses.iterrows()):
                    p_gs = prow.get(pose_score_col, 0)
                    p_vdw = prow.get("Grid_vdw_energy", 0)
                    p_es = prow.get("Grid_es_energy", 0)
                    mode_rows += f"<tr><td>{pi+1}</td><td>{p_gs:.2f}</td><td>{p_vdw:.2f}</td><td>{p_es:.2f}</td></tr>\n"
                card_html += f"""
    <p style="font-size:11px;margin:10px 0 3px 0;color:#555"><b>All binding modes:</b></p>
    <table class="data-table" style="font-size:11px">
    <thead><tr><th>Mode</th><th>Grid_Score</th><th>vdW</th><th>ES</th></tr></thead>
    <tbody>{mode_rows}</tbody>
    </table>
"""

        card_html += "</div>\n"
        detail_cards_html += card_html

    # --- Pharmacophore section ---
    pharma_section = ""
    if pharma_data and pharma_data.get("residues"):
        pharma_items = ""
        for r in pharma_data["residues"]:
            pharma_items += (f"<li><b>{r['residue_id']}</b>: "
                             f"{r['frac_contributing'] * 100:.0f}% of molecules, "
                             f"mean energy = {r['mean_total']:.2f} kcal/mol</li>\n")
        pharma_section = f"<ul>{pharma_items}</ul>"

    # Top candidates summary
    top3_lines = []
    for i in range(min(3, len(df_ranking))):
        row = df_ranking.iloc[i]
        gs = row.get(score_col, 0)
        cs = row.get("composite", 0)
        top3_lines.append(f"<b>{row['Name']}</b> (GS={gs:.1f}, composite={cs:.2f})")
    top3_str = ", ".join(top3_lines) if top3_lines else "N/A"

    # --- Build HTML ---
    html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>DOCK6 Campaign Report: {campaign_id}</title>
<style>
    body {{ font-family: 'Segoe UI', Arial, sans-serif; max-width: 1100px; margin: 0 auto;
           padding: 20px; background: #fafafa; color: #333; }}
    h1 {{ color: #1a5276; border-bottom: 3px solid #1a5276; padding-bottom: 10px; }}
    h2 {{ color: #2c3e50; margin-top: 40px; border-bottom: 1px solid #bdc3c7; padding-bottom: 5px; }}
    h3 {{ color: #34495e; }}
    .summary-box {{ background: #eaf2f8; border-left: 4px solid #2980b9; padding: 15px;
                    margin: 15px 0; border-radius: 4px; }}
    .key-finding {{ background: #e8f8f5; border-left: 4px solid #27ae60; padding: 15px;
                    margin: 15px 0; border-radius: 4px; }}
    .warning {{ background: #fdf2e9; border-left: 4px solid #e67e22; padding: 15px;
               margin: 15px 0; border-radius: 4px; }}
    .data-table {{ border-collapse: collapse; width: 100%; margin: 10px 0; font-size: 12px; }}
    .data-table th {{ background: #2c3e50; color: white; padding: 8px 12px; text-align: left; }}
    .data-table td {{ padding: 6px 12px; border-bottom: 1px solid #ecf0f1; }}
    .data-table tr:nth-child(even) {{ background: #f8f9fa; }}
    .data-table tr:hover {{ background: #ebf5fb; }}
    .figure {{ text-align: center; margin: 20px 0; }}
    .figure img {{ max-width: 100%; border: 1px solid #ddd; border-radius: 4px; }}
    .figure-caption {{ font-size: 11px; color: #777; margin-top: 5px; }}
    .stats {{ display: flex; gap: 20px; flex-wrap: wrap; margin: 15px 0; }}
    .stat-card {{ background: white; border: 1px solid #ddd; border-radius: 8px;
                  padding: 15px; min-width: 150px; text-align: center; }}
    .stat-value {{ font-size: 24px; font-weight: bold; color: #2980b9; }}
    .stat-label {{ font-size: 11px; color: #777; }}
    .footer {{ margin-top: 40px; padding-top: 10px; border-top: 1px solid #ddd;
              font-size: 11px; color: #999; }}
</style>
</head>
<body>

<h1>DOCK6 Campaign Report: {campaign_id}</h1>
<p><b>Engine:</b> DOCK6 (Grid-based scoring) &nbsp; | &nbsp;
   <b>Generated:</b> {datetime.now().strftime('%Y-%m-%d %H:%M')} &nbsp; | &nbsp;
   <b>Molecules:</b> {n_molecules}</p>

<div class="stats">
    <div class="stat-card">
        <div class="stat-value">{n_molecules}</div>
        <div class="stat-label">Molecules</div>
    </div>
    <div class="stat-card">
        <div class="stat-value">{best_score:.1f}</div>
        <div class="stat-label">Best Grid_Score</div>
    </div>
    <div class="stat-card">
        <div class="stat-value">{n_modes_mean:.1f}</div>
        <div class="stat-label">Mean Binding Modes</div>
    </div>
    <div class="stat-card">
        <div class="stat-value">{n_pharmacophore}</div>
        <div class="stat-label">Pharmacophore Residues</div>
    </div>
    <div class="stat-card">
        <div class="stat-value">{n_residues}</div>
        <div class="stat-label">Total Residues Tracked</div>
    </div>
</div>

<!-- ============================================================ -->
<h2>1. Top 20 by Grid_Score (raw)</h2>

<div class="summary-box">
    <b>Grid_Score range:</b> {best_score:.2f} to {worst_score:.2f} kcal/mol<br>
    <b>Campaign mean:</b> {mean_score:.2f} kcal/mol<br>
    <b>Dominance:</b> {n_vdw_dom} vdW-dominated, {n_es_dom} ES-dominated, {n_balanced} balanced<br>
    <b>Top 3:</b> {top3_str}
</div>

{table1_html}

{f'<div class="figure"><img src="data:image/png;base64,{figures["ranking"]}"><div class="figure-caption">Fig 1. Molecules ranked by Grid_Score (vdW + ES).</div></div>' if figures.get("ranking") else ''}

{f'<div class="figure"><img src="data:image/png;base64,{figures["vdw_es"]}"><div class="figure-caption">Fig 2. Energy decomposition: vdW vs electrostatic per molecule. Color = Grid_Score.</div></div>' if figures.get("vdw_es") else ''}

<!-- ============================================================ -->
<h2>2. Top 20 by Mean Grid_Score (composite)</h2>

<div class="summary-box">
    <b>Composite:</b> mean Grid_Score across all binding modes<br>
    <b>Best composite:</b> {best_name} ({best_composite:.2f} kcal/mol)<br>
    <b>Interpretation:</b> More negative = better. The mean penalizes molecules
    that only bind well in one outlier pose and rewards consistent binders.
</div>

{table2_html}

<!-- ============================================================ -->
<h2>3. Top 20 Receptor Residues</h2>

<div class="summary-box">
    <b>Ranked by:</b> n_molecules &times; |mean_energy|<br>
    <b>Pharmacophore residues:</b> {n_pharmacophore} residues contacted by &ge;{(pharma_data.get('threshold', 0.8) * 100) if pharma_data else 80:.0f}% of molecules.
</div>

{pharma_section}

{residue_table_html}

{f'<div class="figure"><img src="data:image/png;base64,{figures["pharmacophore"]}"><div class="figure-caption">Fig 3. Left: top residues by mean energy. Right: fraction of molecules contributing per residue.</div></div>' if figures.get("pharmacophore") else ''}

{f'<div class="figure"><img src="data:image/png;base64,{figures["footprint_ref"]}"><div class="figure-caption">Fig 4. Per-residue footprint: campaign mean vs reference (UDX). Differences highlight selectivity drivers.</div></div>' if figures.get("footprint_ref") else ''}

<!-- ============================================================ -->
<h2>4. Molecule Detail Cards (Top 20 Composite)</h2>

<div class="summary-box">
    <b>Per-molecule breakdown:</b> Top 10 residues by |mean_energy| for each of the top 20 molecules.<br>
    <b>Assessment key:</b>
    &#9733; Anchor (E &lt; -5 and cons &gt; 0.8) &nbsp;|&nbsp;
    &#10003; Reliable (E &lt; -2 and cons &gt; 0.7) &nbsp;|&nbsp;
    &#9889; Strong but variable (E &lt; -5 and cons &lt; 0.5) &nbsp;|&nbsp;
    &#10007; Unreliable (cons &lt; 0.3)
</div>

{detail_cards_html}

<!-- ============================================================ -->
<h2>5. Binding Modes &amp; Contacts</h2>

<div class="summary-box">
    <b>Binding modes per molecule:</b> {n_modes_range} (mean={n_modes_mean:.1f})<br>
    <b>Molecules with multiple modes:</b> {(df_modes['n_modes'] > 1).sum() if df_modes is not None and 'n_modes' in df_modes.columns else 'N/A'}/{n_molecules}
</div>

{f'<div class="figure"><img src="data:image/png;base64,{figures["modes"]}"><div class="figure-caption">Fig 5. Left: binding mode distribution. Right: score spread vs number of modes.</div></div>' if figures.get("modes") else ''}

{f'<div class="figure"><img src="data:image/png;base64,{figures["contacts"]}"><div class="figure-caption">Fig 6. Contact heatmap: top residues x top molecules. Color = min distance (Angstrom). Green = close contact.</div></div>' if figures.get("contacts") else ''}

<div class="footer">
    <p>Generated by molecular_docking — Module 04e (DOCK6 analysis) v4.0</p>
    <p>Reference: Balius et al. J Chem Inf Model 2011, 51(8):1942-56</p>
</div>

</body>
</html>"""

    return html


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def run_campaign_report(
        results_base: Union[str, Path],
        output_dir: Union[str, Path],
        campaign_id: str = "",
        gnina_scores_csv: Optional[str] = None,
        molecules_csv: Optional[str] = None,
        all_poses_csv: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate DOCK6 campaign HTML report with figures.

    Reads outputs from 04a-04d and produces composite ranking + HTML report.
    Composite = mean Grid_Score across all binding modes.
    """
    results_base = Path(results_base)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_path = output_dir / "figures"
    fig_path.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("  04e DOCK6 Campaign Report v4.0")
    logger.info("=" * 60)

    # --- Locate input files ---
    score_csv = results_base / "04a_score_ranking" / "score_components.csv"
    ranking_csv = results_base / "04a_score_ranking" / "molecule_ranking.csv"
    consensus_csv = results_base / "04b_footprint_analysis" / "residue_consensus.csv"
    pharma_json = results_base / "04b_footprint_analysis" / "pharmacophore_residues.json"
    fps_summary_csv = results_base / "04b_footprint_analysis" / "molecule_footprint_summary.csv"
    consistency_csv = results_base / "04b_footprint_analysis" / "molecule_consistency_summary.csv"
    residue_per_pose_csv = results_base / "04b_footprint_analysis" / "residue_per_pose.csv"
    modes_csv = results_base / "04c_binding_modes" / "binding_modes_summary.csv"
    contacts_csv = results_base / "04d_contact_mapping" / "contact_summary.csv"

    # Auto-detect all_poses_csv from 01e
    if not all_poses_csv:
        # Look in parent of results_base (04_dock6_analysis) → go to 01e
        campaign_results = results_base.parent
        auto_poses = campaign_results / "01e_score_collection" / "dock6_all_poses.csv"
        if auto_poses.exists():
            all_poses_csv = str(auto_poses)
            logger.info(f"  Auto-detected all_poses_csv: {all_poses_csv}")

    # Use score_components if available, else molecule_ranking
    main_csv = str(score_csv) if score_csv.exists() else str(ranking_csv) if ranking_csv.exists() else None

    if not main_csv:
        return {"success": False,
                "error": f"No score data found in {results_base / '04a_score_ranking'}"}

    # --- Build composite ranking ---
    logger.info("  Building composite ranking (mean Grid_Score across modes)...")
    df_ranking = _build_composite_ranking(
        main_csv,
        all_poses_csv=all_poses_csv,
    )

    composite_out = output_dir / "composite_ranking.csv"
    df_ranking.to_csv(composite_out, index=False, encoding="utf-8")
    logger.info(f"  Saved: {composite_out} ({len(df_ranking)} molecules)")

    # --- Load drug-likeness data ---
    df_drug = None
    if gnina_scores_csv and Path(gnina_scores_csv).exists():
        try:
            df_drug = pd.read_csv(gnina_scores_csv)
            logger.info(f"  Loaded drug-likeness from gnina_scores: {gnina_scores_csv}")
        except Exception as e:
            logger.warning(f"  Could not load gnina_scores_csv: {e}")
    if df_drug is None and molecules_csv and Path(molecules_csv).exists():
        try:
            df_drug = pd.read_csv(molecules_csv)
            logger.info(f"  Loaded drug-likeness from molecules_csv: {molecules_csv}")
        except Exception as e:
            logger.warning(f"  Could not load molecules_csv: {e}")

    # --- Load supplementary data ---
    df_residues = pd.read_csv(consensus_csv) if consensus_csv.exists() else None
    pharma_data = None
    if pharma_json.exists():
        with open(pharma_json) as f:
            pharma_data = json.load(f)
    df_modes = pd.read_csv(modes_csv) if modes_csv.exists() else None
    df_contacts = pd.read_csv(contacts_csv) if contacts_csv.exists() else None
    df_fps_summary = pd.read_csv(fps_summary_csv) if fps_summary_csv.exists() else None
    df_residue_per_pose = pd.read_csv(residue_per_pose_csv) if residue_per_pose_csv.exists() else None
    df_all_poses = pd.read_csv(all_poses_csv) if all_poses_csv and Path(all_poses_csv).exists() else None

    # --- Generate figures ---
    logger.info("  Generating figures...")
    figures = {}

    if HAS_MPL:
        figures["ranking"] = fig_molecule_ranking(df_ranking)
        figures["vdw_es"] = fig_vdw_vs_es(df_ranking)
        figures["pharmacophore"] = fig_pharmacophore_residues(df_residues, pharma_data)
        figures["modes"] = fig_binding_modes(df_modes)
        figures["contacts"] = fig_contact_heatmap(df_contacts)
        figures["footprint_ref"] = fig_footprint_vs_reference(df_residues)

        # Save as PNG files too
        for name, b64 in figures.items():
            if b64:
                img_data = base64.b64decode(b64)
                with open(fig_path / f"{name}.png", "wb") as f:
                    f.write(img_data)
        logger.info(f"  Figures: {sum(1 for v in figures.values() if v)}/{len(figures)} generated")
    else:
        logger.warning("  matplotlib not available — no figures")

    # --- Generate HTML ---
    logger.info("  Generating HTML report...")
    html = _generate_html(
        df_ranking=df_ranking,
        df_residues=df_residues,
        pharma_data=pharma_data,
        df_modes=df_modes,
        df_contacts=df_contacts,
        df_fps_summary=df_fps_summary,
        figures=figures,
        campaign_id=campaign_id,
        n_molecules=len(df_ranking),
        df_drug=df_drug,
        df_residue_per_pose=df_residue_per_pose,
        df_all_poses=df_all_poses,
    )

    report_html = output_dir / "campaign_report.html"
    with open(report_html, "w", encoding="utf-8") as f:
        f.write(html)
    logger.info(f"  Saved: {report_html}")

    # --- Log top 5 ---
    logger.info("")
    logger.info("  Top 5 Composite Ranking:")
    for _, row in df_ranking.head(5).iterrows():
        score_col = "Grid_Score" if "Grid_Score" in df_ranking.columns else "composite"
        logger.info(f"    {int(row['Rank'])}. {row['Name']} "
                     f"(GS={row.get(score_col, 0):.2f}, composite={row['composite']:.2f})")

    logger.info(f"\n{'=' * 60}")
    logger.info("  DOCK6 CAMPAIGN REPORT COMPLETE")
    logger.info(f"{'=' * 60}")
    logger.info(f"  Report: {report_html}")
    logger.info(f"  Figures: {fig_path}/ ({sum(1 for v in figures.values() if v)} images)")

    return {
        "success": True,
        "report_path": str(report_html),
        "composite_csv": str(composite_out),
        "n_figures": sum(1 for v in figures.values() if v),
        "output_dir": str(output_dir),
    }
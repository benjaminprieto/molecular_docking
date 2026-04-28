"""
Cross-Engine Comparison (07a) — DOCK6 vs GNINA consensus analysis
====================================================================
Merges results from both docking engines to identify molecules validated
by independent scoring functions. Molecules that rank well in BOTH engines
are the most confident candidates for purchase/testing.

Rationale:
    DOCK6 (Grid_Score): force-field based, in-vacuo electrostatics, ~88% ES-driven.
        Strength: per-residue energy decomposition (footprint).
        Weakness: overvalues charged species without solvent screening.
    GNINA (Vina): empirical, trained on binding data, balanced terms.
        Strength: fragment-level decomposition via atom_terms.
        Weakness: CNN score unreliable for nucleotide-sugar pockets.

    Molecules in the intersection of both rankings are less likely to be
    artifacts of a single scoring function's biases.

Inputs:
    - gnina_scores.csv (02c): per-molecule Vina + CNN scores + drug-likeness
    - gnina composite_ranking.csv (05g): Vina × weighted_convergence
    - dock6_scores.csv (01e): per-molecule best Grid_Score
    - dock6_all_poses.csv (01e): all poses for mean Grid_Score
    - dock6 composite_ranking.csv (04e): mean Grid_Score ranking

Outputs:
    - cross_engine_report.html: comparison report with tables + figures
    - cross_engine_ranking.csv: master table with all molecules and both engine ranks
    - consensus_hits.csv: molecules appearing in top N of both engines

Location: 01_src/molecular_docking/m07_cross_engine/cross_engine_comparison.py
Project: molecular_docking
Module: 07a (core)
Version: 1.0 (2026-03-27)
"""

import json
import logging
import base64
from io import BytesIO
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from datetime import datetime

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

try:
    from scipy.stats import spearmanr
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


# =============================================================================
# FIGURE HELPERS
# =============================================================================

def _fig_to_base64(fig) -> str:
    buf = BytesIO()
    fig.savefig(buf, format='png', dpi=120, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode('utf-8')


# =============================================================================
# DATA LOADING
# =============================================================================

def _load_gnina_data(
        gnina_scores_csv: str,
        gnina_composite_csv: Optional[str] = None,
) -> pd.DataFrame:
    """
    Load GNINA results into a unified DataFrame.
    
    Returns DataFrame with columns:
        Name, vina_affinity, cnn_score, cnn_affinity,
        gnina_composite, gnina_weighted_conv,
        MW, QED, TPSA, LogP, ...
    """
    df = pd.read_csv(gnina_scores_csv)
    
    # Normalize name column
    if 'name' in df.columns and 'Name' not in df.columns:
        df = df.rename(columns={'name': 'Name'})
    
    # Filter successful only
    if 'success' in df.columns:
        df = df[df['success'] == True].copy()
    
    # Add Vina rank
    df = df.sort_values('vina_affinity', ascending=True).reset_index(drop=True)
    df['gnina_vina_rank'] = range(1, len(df) + 1)
    
    # Merge composite ranking if available
    if gnina_composite_csv and Path(gnina_composite_csv).exists():
        df_comp = pd.read_csv(gnina_composite_csv)
        # Find the composite column (could be named differently)
        comp_col = None
        for c in ['composite', 'composite_score']:
            if c in df_comp.columns:
                comp_col = c
                break
        
        name_col = 'name' if 'name' in df_comp.columns else 'Name'
        if comp_col and name_col in df_comp.columns:
            # Rank by composite
            df_comp = df_comp.sort_values(comp_col, ascending=True).reset_index(drop=True)
            df_comp['gnina_composite_rank'] = range(1, len(df_comp) + 1)
            
            merge_cols = [name_col, comp_col, 'gnina_composite_rank']
            # Also grab weighted_conv if present
            for extra in ['weighted_conv', 'weighted_convergence']:
                if extra in df_comp.columns:
                    merge_cols.append(extra)
            
            df_comp = df_comp[merge_cols].rename(columns={
                name_col: 'Name',
                comp_col: 'gnina_composite',
            })
            df = df.merge(df_comp, on='Name', how='left')
    
    if 'gnina_composite_rank' not in df.columns:
        df['gnina_composite_rank'] = np.nan
    if 'gnina_composite' not in df.columns:
        df['gnina_composite'] = np.nan
    
    return df


def _load_dock6_data(
        dock6_scores_csv: str,
        dock6_all_poses_csv: Optional[str] = None,
        dock6_composite_csv: Optional[str] = None,
) -> pd.DataFrame:
    """
    Load DOCK6 results into a unified DataFrame.
    
    Returns DataFrame with columns:
        Name, best_Grid_Score, mean_Grid_Score, n_modes, spread, ratio,
        dock6_best_rank, dock6_mean_rank
    """
    df = pd.read_csv(dock6_scores_csv)
    
    # Normalize
    if 'name' in df.columns and 'Name' not in df.columns:
        df = df.rename(columns={'name': 'Name'})
    
    score_col = 'Grid_Score' if 'Grid_Score' in df.columns else df.columns[1]
    df = df.rename(columns={score_col: 'best_Grid_Score'})
    
    # Rank by best
    df = df.sort_values('best_Grid_Score', ascending=True).reset_index(drop=True)
    df['dock6_best_rank'] = range(1, len(df) + 1)
    
    # Compute mean from all_poses
    if dock6_all_poses_csv and Path(dock6_all_poses_csv).exists():
        df_poses = pd.read_csv(dock6_all_poses_csv)
        pose_score = 'Grid_Score' if 'Grid_Score' in df_poses.columns else df_poses.columns[1]
        
        pose_stats = df_poses.groupby('Name').agg(
            mean_Grid_Score=(pose_score, 'mean'),
            n_modes=(pose_score, 'count'),
            best_pose=(pose_score, 'min'),
            worst_Grid_Score=(pose_score, 'max'),
        ).reset_index()

        pose_stats['spread'] = pose_stats['worst_Grid_Score'] - pose_stats['best_pose']
        pose_stats['ratio'] = np.where(
            pose_stats['best_pose'] != 0,
            pose_stats['mean_Grid_Score'] / pose_stats['best_pose'],
            0.0,
        )
        pose_stats = pose_stats.drop(columns=['best_pose', 'worst_Grid_Score'])
        
        df = df.merge(pose_stats[['Name', 'mean_Grid_Score', 'n_modes', 'spread', 'ratio']],
                      on='Name', how='left')
        
        # Rank by mean
        df = df.sort_values('mean_Grid_Score', ascending=True).reset_index(drop=True)
        df['dock6_mean_rank'] = range(1, len(df) + 1)
    else:
        df['mean_Grid_Score'] = df['best_Grid_Score']
        df['n_modes'] = 1
        df['spread'] = 0.0
        df['ratio'] = 1.0
        df['dock6_mean_rank'] = df['dock6_best_rank']
    
    # Merge composite ranking if available
    if dock6_composite_csv and Path(dock6_composite_csv).exists():
        df_comp = pd.read_csv(dock6_composite_csv)
        name_col = 'name' if 'name' in df_comp.columns else 'Name'
        comp_col = None
        for c in ['composite', 'composite_score', 'mean_Grid_Score']:
            if c in df_comp.columns:
                comp_col = c
                break
        if comp_col and name_col in df_comp.columns:
            rank_col_name = 'Rank' if 'Rank' in df_comp.columns else 'rank'
            if rank_col_name in df_comp.columns:
                df_comp_merge = df_comp[[name_col, rank_col_name]].rename(columns={
                    name_col: 'Name',
                    rank_col_name: 'dock6_composite_rank',
                })
                df = df.merge(df_comp_merge, on='Name', how='left')
    
    if 'dock6_composite_rank' not in df.columns:
        df['dock6_composite_rank'] = df['dock6_mean_rank']
    
    return df


# =============================================================================
# CROSS-ENGINE MERGE
# =============================================================================

def merge_engines(
        df_gnina: pd.DataFrame,
        df_dock6: pd.DataFrame,
) -> pd.DataFrame:
    """
    Merge GNINA and DOCK6 results into a master DataFrame.
    
    All molecules from both engines, with ranks from each.
    Molecules only in one engine get NaN for the other's ranks.
    """
    # Select GNINA columns
    gnina_cols = ['Name', 'vina_affinity', 'gnina_vina_rank', 'gnina_composite',
                  'gnina_composite_rank']
    # Add drug-likeness if present
    for col in ['MW', 'QED', 'TPSA', 'LogP', 'cnn_score', 'cnn_affinity',
                'weighted_conv', 'weighted_convergence']:
        if col in df_gnina.columns:
            gnina_cols.append(col)
    gnina_cols = [c for c in gnina_cols if c in df_gnina.columns]
    
    # Select DOCK6 columns
    dock6_cols = ['Name', 'best_Grid_Score', 'mean_Grid_Score', 'n_modes',
                  'spread', 'ratio', 'dock6_best_rank', 'dock6_mean_rank']
    dock6_cols = [c for c in dock6_cols if c in df_dock6.columns]
    
    # Outer merge
    df = df_gnina[gnina_cols].merge(df_dock6[dock6_cols], on='Name', how='outer')
    
    return df


def compute_consensus(
        df_merged: pd.DataFrame,
        top_n: int = 20,
) -> pd.DataFrame:
    """
    Compute consensus score: how many top-N lists each molecule appears in.
    
    Lists checked:
        1. GNINA Vina top N (gnina_vina_rank <= top_n)
        2. GNINA Composite top N (gnina_composite_rank <= top_n)
        3. DOCK6 Best Grid_Score top N (dock6_best_rank <= top_n)
        4. DOCK6 Mean Grid_Score top N (dock6_mean_rank <= top_n)
    
    Returns df_merged with additional columns:
        n_lists: 0-4, how many top-N lists the molecule appears in
        in_gnina_vina, in_gnina_composite, in_dock6_best, in_dock6_mean: bool
        consensus_rank: ordered by n_lists desc, then by average rank
    """
    df = df_merged.copy()
    
    df['in_gnina_vina'] = df['gnina_vina_rank'] <= top_n
    df['in_gnina_composite'] = df['gnina_composite_rank'] <= top_n
    df['in_dock6_best'] = df['dock6_best_rank'] <= top_n
    df['in_dock6_mean'] = df['dock6_mean_rank'] <= top_n
    
    df['n_lists'] = (
        df['in_gnina_vina'].astype(int) +
        df['in_gnina_composite'].astype(int) +
        df['in_dock6_best'].astype(int) +
        df['in_dock6_mean'].astype(int)
    )
    
    # Average rank across available engines (lower = better)
    rank_cols = ['gnina_vina_rank', 'gnina_composite_rank',
                 'dock6_best_rank', 'dock6_mean_rank']
    df['avg_rank'] = df[rank_cols].mean(axis=1, skipna=True)
    
    # Best rank in any list
    df['best_rank'] = df[rank_cols].min(axis=1, skipna=True)
    
    # Consensus rank: sort by n_lists (desc), then avg_rank (asc)
    df = df.sort_values(['n_lists', 'avg_rank'], ascending=[False, True]).reset_index(drop=True)
    df['consensus_rank'] = range(1, len(df) + 1)
    
    return df


# =============================================================================
# FIGURES
# =============================================================================

def fig_vina_vs_gridscore(df: pd.DataFrame) -> str:
    """Scatter plot: GNINA Vina affinity vs DOCK6 mean Grid_Score."""
    if not HAS_MPL:
        return ""
    
    # Only molecules with both scores
    both = df.dropna(subset=['vina_affinity', 'mean_Grid_Score'])
    if len(both) < 5:
        return ""
    
    fig, ax = plt.subplots(figsize=(8, 7))
    
    # Color by n_lists
    colors = both['n_lists'].values
    sc = ax.scatter(both['vina_affinity'], both['mean_Grid_Score'],
                    c=colors, cmap='RdYlGn', s=30, alpha=0.6,
                    edgecolors='gray', linewidth=0.3, vmin=0, vmax=4)
    plt.colorbar(sc, ax=ax, label='# Top-20 lists (0-4)', ticks=[0, 1, 2, 3, 4])
    
    # Label consensus hits (n_lists >= 2)
    consensus = both[both['n_lists'] >= 2]
    for _, row in consensus.iterrows():
        name = row['Name']
        short = name[:20] if len(name) > 20 else name
        ax.annotate(short, (row['vina_affinity'], row['mean_Grid_Score']),
                    fontsize=6, alpha=0.8, ha='left',
                    xytext=(5, 5), textcoords='offset points')
    
    # Correlation
    corr = both[['vina_affinity', 'mean_Grid_Score']].corr().iloc[0, 1]
    
    ax.set_xlabel('GNINA Vina Affinity (kcal/mol)')
    ax.set_ylabel('DOCK6 Mean Grid_Score (kcal/mol)')
    ax.set_title(f'Cross-Engine Score Correlation (r={corr:.3f}, n={len(both)})')
    
    # Trend line
    z = np.polyfit(both['vina_affinity'], both['mean_Grid_Score'], 1)
    x_line = np.linspace(both['vina_affinity'].min(), both['vina_affinity'].max(), 50)
    ax.plot(x_line, np.polyval(z, x_line), 'k--', alpha=0.3)
    
    return _fig_to_base64(fig)


def fig_rank_comparison(df: pd.DataFrame) -> str:
    """Scatter: GNINA Vina rank vs DOCK6 mean rank. Consensus = diagonal."""
    if not HAS_MPL:
        return ""
    
    both = df.dropna(subset=['gnina_vina_rank', 'dock6_mean_rank'])
    if len(both) < 5:
        return ""
    
    fig, ax = plt.subplots(figsize=(8, 7))
    
    colors = both['n_lists'].values
    sc = ax.scatter(both['gnina_vina_rank'], both['dock6_mean_rank'],
                    c=colors, cmap='RdYlGn', s=25, alpha=0.5,
                    edgecolors='gray', linewidth=0.3, vmin=0, vmax=4)
    plt.colorbar(sc, ax=ax, label='# Top-20 lists')
    
    # Diagonal = perfect agreement
    max_rank = max(both['gnina_vina_rank'].max(), both['dock6_mean_rank'].max())
    ax.plot([0, max_rank], [0, max_rank], 'k:', alpha=0.3, label='Perfect agreement')
    
    # Top-20 box
    ax.axhline(y=20, color='red', linestyle='--', alpha=0.3)
    ax.axvline(x=20, color='red', linestyle='--', alpha=0.3)
    ax.fill_between([0, 20], [0, 0], [20, 20], alpha=0.05, color='green')
    ax.text(10, 10, 'CONSENSUS\nZONE', ha='center', va='center',
            fontsize=10, color='green', alpha=0.4, weight='bold')
    
    # Label consensus
    consensus = both[both['n_lists'] >= 3]
    for _, row in consensus.iterrows():
        short = row['Name'][:18]
        ax.annotate(short, (row['gnina_vina_rank'], row['dock6_mean_rank']),
                    fontsize=6, alpha=0.8)
    
    # Rank correlation
    if HAS_SCIPY:
        try:
            rho, pval = spearmanr(both['gnina_vina_rank'], both['dock6_mean_rank'])
            ax.set_title(f'Rank Agreement (Spearman ρ={rho:.3f}, p={pval:.2e})')
        except Exception:
            ax.set_title('Rank Agreement')
    else:
        ax.set_title('Rank Agreement')
    
    ax.set_xlabel('GNINA Vina Rank')
    ax.set_ylabel('DOCK6 Mean Grid_Score Rank')
    ax.set_xlim(0, max_rank * 1.05)
    ax.set_ylim(0, max_rank * 1.05)
    
    return _fig_to_base64(fig)


def fig_overlap_venn(df: pd.DataFrame, top_n: int = 20) -> str:
    """Bar chart showing overlap between the 4 top-N lists."""
    if not HAS_MPL:
        return ""
    
    lists = {
        'GNINA Vina': set(df[df['in_gnina_vina']]['Name']),
        'GNINA Composite': set(df[df['in_gnina_composite']]['Name']),
        'DOCK6 Best': set(df[df['in_dock6_best']]['Name']),
        'DOCK6 Mean': set(df[df['in_dock6_mean']]['Name']),
    }
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Left: molecules by n_lists
    counts = df['n_lists'].value_counts().sort_index()
    ax = axes[0]
    bars = ax.bar(counts.index, counts.values, color=['#e74c3c', '#f39c12', '#f1c40f', '#2ecc71', '#27ae60'])
    ax.set_xlabel('Number of Top-20 Lists')
    ax.set_ylabel('Number of Molecules')
    ax.set_title(f'Distribution of Consensus (top {top_n})')
    ax.set_xticks([0, 1, 2, 3, 4])
    for bar, count in zip(bars, counts.values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                str(count), ha='center', fontsize=10)
    
    # Right: pairwise overlaps
    ax = axes[1]
    pairs = [
        ('GNINA Vina', 'DOCK6 Mean'),
        ('GNINA Comp.', 'DOCK6 Best'),
        ('GNINA Vina', 'GNINA Comp.'),
        ('DOCK6 Best', 'DOCK6 Mean'),
        ('GNINA Vina', 'DOCK6 Best'),
        ('GNINA Comp.', 'DOCK6 Mean'),
    ]
    pair_keys = [
        ('GNINA Vina', 'DOCK6 Mean'),
        ('GNINA Composite', 'DOCK6 Best'),
        ('GNINA Vina', 'GNINA Composite'),
        ('DOCK6 Best', 'DOCK6 Mean'),
        ('GNINA Vina', 'DOCK6 Best'),
        ('GNINA Composite', 'DOCK6 Mean'),
    ]
    overlaps = []
    for (n1, n2), (k1, k2) in zip(pairs, pair_keys):
        overlap = len(lists.get(k1, set()) & lists.get(k2, set()))
        overlaps.append((f'{n1}\n∩\n{n2}', overlap))
    
    labels, values = zip(*overlaps)
    colors = ['#3498db' if 'GNINA' in l and 'DOCK6' in l else '#95a5a6' for l in labels]
    bars = ax.bar(range(len(values)), values, color=colors)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel('Overlap (molecules)')
    ax.set_title(f'Pairwise Overlap (top {top_n})')
    ax.set_ylim(0, top_n + 2)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                str(val), ha='center', fontsize=10, weight='bold')
    
    fig.tight_layout()
    return _fig_to_base64(fig)


# =============================================================================
# HTML REPORT
# =============================================================================

def _html_table(df: pd.DataFrame, max_rows: int = 50) -> str:
    if df.empty:
        return "<p><i>No data</i></p>"
    display_df = df.head(max_rows)
    html = display_df.to_html(index=False, classes='data-table', border=0,
                               float_format=lambda x: f'{x:.2f}' if isinstance(x, float) else str(x))
    if len(df) > max_rows:
        html += f"<p><i>Showing {max_rows} of {len(df)} rows</i></p>"
    return html


def _stability_badge(val, label="") -> str:
    if pd.isna(val):
        return '<span style="color:#999">—</span>'
    if isinstance(val, (int, float)):
        return f'{val:.0f}'
    return str(val)


def generate_html_report(
        df_consensus: pd.DataFrame,
        figures: Dict[str, str],
        campaign_id: str,
        top_n: int = 20,
        n_gnina: int = 0,
        n_dock6: int = 0,
) -> str:
    """Generate cross-engine comparison HTML report."""
    
    n_total = len(df_consensus)
    n_both = len(df_consensus.dropna(subset=['vina_affinity', 'best_Grid_Score']))
    n_consensus_2plus = len(df_consensus[df_consensus['n_lists'] >= 2])
    n_consensus_3plus = len(df_consensus[df_consensus['n_lists'] >= 3])
    n_consensus_4 = len(df_consensus[df_consensus['n_lists'] == 4])
    
    # Top consensus molecules
    consensus_hits = df_consensus[df_consensus['n_lists'] >= 2].head(20)
    
    # Build consensus table
    consensus_cols = ['consensus_rank', 'Name', 'n_lists', 'vina_affinity',
                      'gnina_vina_rank', 'gnina_composite_rank',
                      'mean_Grid_Score', 'dock6_mean_rank', 'best_Grid_Score',
                      'dock6_best_rank']
    # Add drug-likeness
    for dc in ['MW', 'QED', 'TPSA']:
        if dc in consensus_hits.columns:
            consensus_cols.append(dc)
    consensus_cols = [c for c in consensus_cols if c in consensus_hits.columns]
    df_consensus_display = consensus_hits[consensus_cols].copy()
    
    # Rename for readability
    rename_map = {
        'consensus_rank': 'Rank',
        'n_lists': 'Lists',
        'vina_affinity': 'Vina',
        'gnina_vina_rank': 'G.Vina Rk',
        'gnina_composite_rank': 'G.Comp Rk',
        'mean_Grid_Score': 'Mean GS',
        'dock6_mean_rank': 'D6.Mean Rk',
        'best_Grid_Score': 'Best GS',
        'dock6_best_rank': 'D6.Best Rk',
    }
    df_consensus_display = df_consensus_display.rename(columns=rename_map)
    
    # Build engine-specific top 20 tables
    # GNINA only (not in any DOCK6 top 20)
    gnina_only = df_consensus[
        (df_consensus['in_gnina_vina'] | df_consensus['in_gnina_composite']) &
        ~df_consensus['in_dock6_best'] & ~df_consensus['in_dock6_mean']
    ].head(20)
    
    gnina_only_cols = ['Name', 'vina_affinity', 'gnina_vina_rank', 'gnina_composite_rank',
                       'mean_Grid_Score', 'dock6_mean_rank']
    gnina_only_cols = [c for c in gnina_only_cols if c in gnina_only.columns]
    
    # DOCK6 only (not in any GNINA top 20)
    dock6_only = df_consensus[
        (df_consensus['in_dock6_best'] | df_consensus['in_dock6_mean']) &
        ~df_consensus['in_gnina_vina'] & ~df_consensus['in_gnina_composite']
    ].head(20)
    
    dock6_only_cols = ['Name', 'best_Grid_Score', 'mean_Grid_Score', 'dock6_mean_rank',
                       'dock6_best_rank', 'vina_affinity', 'gnina_vina_rank']
    dock6_only_cols = [c for c in dock6_only_cols if c in dock6_only.columns]
    
    # List membership detail for top molecules
    def _list_badges(row):
        badges = []
        if row.get('in_gnina_vina', False):
            badges.append('<span style="background:#3498db;color:white;padding:2px 6px;border-radius:3px;font-size:10px">G.Vina</span>')
        if row.get('in_gnina_composite', False):
            badges.append('<span style="background:#2980b9;color:white;padding:2px 6px;border-radius:3px;font-size:10px">G.Comp</span>')
        if row.get('in_dock6_best', False):
            badges.append('<span style="background:#e67e22;color:white;padding:2px 6px;border-radius:3px;font-size:10px">D6.Best</span>')
        if row.get('in_dock6_mean', False):
            badges.append('<span style="background:#d35400;color:white;padding:2px 6px;border-radius:3px;font-size:10px">D6.Mean</span>')
        return ' '.join(badges)
    
    # Build consensus table HTML manually for badges
    consensus_table_html = """
<table class="data-table">
<thead><tr>
    <th>Rank</th><th>Molecule</th><th>Lists</th><th>In</th>
    <th>Vina</th><th>G.Vina Rk</th><th>G.Comp Rk</th>
    <th>Mean GS</th><th>D6.Mean Rk</th>
    <th>MW</th><th>QED</th>
</tr></thead>
<tbody>
"""
    for _, row in consensus_hits.iterrows():
        badges = _list_badges(row)
        vina = f"{row['vina_affinity']:.2f}" if pd.notna(row.get('vina_affinity')) else "—"
        gvr = f"{int(row['gnina_vina_rank'])}" if pd.notna(row.get('gnina_vina_rank')) else "—"
        gcr = f"{int(row['gnina_composite_rank'])}" if pd.notna(row.get('gnina_composite_rank')) else "—"
        mgs = f"{row['mean_Grid_Score']:.1f}" if pd.notna(row.get('mean_Grid_Score')) else "—"
        dmr = f"{int(row['dock6_mean_rank'])}" if pd.notna(row.get('dock6_mean_rank')) else "—"
        mw = f"{row['MW']:.0f}" if pd.notna(row.get('MW')) else "—"
        qed = f"{row['QED']:.3f}" if pd.notna(row.get('QED')) else "—"
        
        consensus_table_html += f"""<tr>
    <td>{int(row['consensus_rank'])}</td>
    <td><b>{row['Name']}</b></td>
    <td style="text-align:center"><b>{int(row['n_lists'])}</b></td>
    <td>{badges}</td>
    <td>{vina}</td><td>{gvr}</td><td>{gcr}</td>
    <td>{mgs}</td><td>{dmr}</td>
    <td>{mw}</td><td>{qed}</td>
</tr>
"""
    consensus_table_html += "</tbody></table>"
    
    # Correlation stats
    both_scores = df_consensus.dropna(subset=['vina_affinity', 'mean_Grid_Score'])
    if len(both_scores) > 5:
        pearson_r = both_scores[['vina_affinity', 'mean_Grid_Score']].corr().iloc[0, 1]
        if HAS_SCIPY:
            try:
                spearman_rho, spearman_p = spearmanr(
                    both_scores['gnina_vina_rank'].dropna(),
                    both_scores['dock6_mean_rank'].dropna()
                )
            except Exception:
                spearman_rho, spearman_p = 0, 1
        else:
            spearman_rho, spearman_p = 0, 1
    else:
        pearson_r = 0
        spearman_rho, spearman_p = 0, 1
    
    html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Cross-Engine Comparison: {campaign_id}</title>
<style>
    body {{ font-family: 'Segoe UI', Arial, sans-serif; max-width: 1200px; margin: 0 auto;
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
    .data-table th {{ background: #2c3e50; color: white; padding: 8px 10px; text-align: left; }}
    .data-table td {{ padding: 6px 10px; border-bottom: 1px solid #ecf0f1; }}
    .data-table tr:nth-child(even) {{ background: #f8f9fa; }}
    .data-table tr:hover {{ background: #ebf5fb; }}
    .figure {{ text-align: center; margin: 20px 0; }}
    .figure img {{ max-width: 100%; border: 1px solid #ddd; border-radius: 4px; }}
    .stats {{ display: flex; gap: 20px; flex-wrap: wrap; margin: 15px 0; }}
    .stat-card {{ background: white; border: 1px solid #ddd; border-radius: 8px;
                  padding: 15px; min-width: 130px; text-align: center; }}
    .stat-value {{ font-size: 24px; font-weight: bold; color: #2980b9; }}
    .stat-label {{ font-size: 11px; color: #777; }}
    .footer {{ margin-top: 40px; padding-top: 10px; border-top: 1px solid #ddd;
              font-size: 11px; color: #999; }}
</style>
</head>
<body>

<h1>Cross-Engine Comparison: {campaign_id}</h1>
<p><b>Generated:</b> {datetime.now().strftime('%Y-%m-%d %H:%M')} &nbsp;|&nbsp;
   <b>GNINA molecules:</b> {n_gnina} &nbsp;|&nbsp;
   <b>DOCK6 molecules:</b> {n_dock6} &nbsp;|&nbsp;
   <b>In both engines:</b> {n_both}</p>

<div class="stats">
    <div class="stat-card">
        <div class="stat-value">{n_both}</div>
        <div class="stat-label">Molecules in Both</div>
    </div>
    <div class="stat-card">
        <div class="stat-value">{n_consensus_2plus}</div>
        <div class="stat-label">In 2+ Top-{top_n} Lists</div>
    </div>
    <div class="stat-card">
        <div class="stat-value">{n_consensus_3plus}</div>
        <div class="stat-label">In 3+ Top-{top_n} Lists</div>
    </div>
    <div class="stat-card">
        <div class="stat-value">{n_consensus_4}</div>
        <div class="stat-label">In All 4 Lists</div>
    </div>
    <div class="stat-card">
        <div class="stat-value">{pearson_r:.3f}</div>
        <div class="stat-label">Score Correlation (r)</div>
    </div>
    <div class="stat-card">
        <div class="stat-value">{spearman_rho:.3f}</div>
        <div class="stat-label">Rank Correlation (ρ)</div>
    </div>
</div>

<!-- ============================================================ -->
<h2>1. Consensus Hits (molecules in 2+ top-{top_n} lists)</h2>

<div class="summary-box">
    <b>Cross-engine validation:</b> {n_consensus_2plus} molecules appear in at least 2 of the 4 top-{top_n} lists
    (GNINA Vina, GNINA Composite, DOCK6 Best, DOCK6 Mean).<br>
    <b>Strongest consensus:</b> {n_consensus_3plus} molecules in 3+ lists, {n_consensus_4} in all 4.<br>
    <b>Interpretation:</b> Molecules in more lists are validated by independent scoring functions
    and are less likely to be artifacts of a single engine's biases.
</div>

{consensus_table_html}

{f'<div class="figure"><img src="data:image/png;base64,{figures["overlap"]}"><div class="figure-caption">Fig 1. Left: distribution of molecules by number of top-{top_n} lists. Right: pairwise overlap between lists. Blue = cross-engine pairs.</div></div>' if figures.get("overlap") else ''}

<!-- ============================================================ -->
<h2>2. Score Correlation</h2>

<div class="summary-box">
    <b>Pearson r (scores):</b> {pearson_r:.3f} — {"strong" if abs(pearson_r) > 0.7 else "moderate" if abs(pearson_r) > 0.4 else "weak"} correlation<br>
    <b>Spearman ρ (ranks):</b> {spearman_rho:.3f} (p={spearman_p:.2e})<br>
    <b>Note:</b> Low correlation is expected — DOCK6 Grid_Score is force-field based (in-vacuo, ES-dominated)
    while GNINA Vina is empirical (balanced terms). They evaluate different molecular properties.
</div>

{f'<div class="figure"><img src="data:image/png;base64,{figures["scatter"]}"><div class="figure-caption">Fig 2. GNINA Vina affinity vs DOCK6 mean Grid_Score. Color = number of top-{top_n} lists. Labeled = consensus hits (2+ lists).</div></div>' if figures.get("scatter") else ''}

{f'<div class="figure"><img src="data:image/png;base64,{figures["ranks"]}"><div class="figure-caption">Fig 3. Rank comparison. Green zone = consensus (both engines rank in top {top_n}). Diagonal = perfect agreement.</div></div>' if figures.get("ranks") else ''}

<!-- ============================================================ -->
<h2>3. Engine-Specific Hits</h2>

<div class="warning">
    <b>Divergence analysis:</b> Molecules in the top {top_n} of one engine but not the other.
    These may represent genuine chemical diversity (different properties recognized by different
    scoring functions) or artifacts of scoring biases.
</div>

<h3>GNINA-only hits (not in any DOCK6 top {top_n})</h3>
<p style="font-size:12px;color:#777">Molecules favored by Vina/CNN but not by force-field Grid_Score.
Likely: balanced interactions, H-bond/hydrophobic dominated.</p>
{_html_table(gnina_only[gnina_only_cols]) if not gnina_only.empty else '<p><i>None</i></p>'}

<h3>DOCK6-only hits (not in any GNINA top {top_n})</h3>
<p style="font-size:12px;color:#777">Molecules favored by Grid_Score but not Vina.
Likely: charged species, strong electrostatics, potentially solvent-sensitivity artifacts.</p>
{_html_table(dock6_only[dock6_only_cols]) if not dock6_only.empty else '<p><i>None</i></p>'}

<!-- ============================================================ -->
<h2>4. Key Findings</h2>

<div class="key-finding">
<h3>Scoring Function Agreement</h3>
<p>Pearson correlation between Vina affinity and mean Grid_Score is <b>{pearson_r:.3f}</b>
({"indicating the two engines largely agree on relative binding strength" if abs(pearson_r) > 0.5 else "indicating significant divergence between engines — expected for force-field vs empirical scoring"}).
Spearman rank correlation is <b>{spearman_rho:.3f}</b>.</p>
</div>

<div class="key-finding">
<h3>Consensus Candidates</h3>
<p><b>{n_consensus_2plus} molecules</b> appear in 2+ top-{top_n} lists. These are the most reliable
candidates for experimental validation — recognized as good binders by independent scoring approaches.
{"The strongest consensus hit(s) appear in " + str(n_consensus_3plus) + "+ lists." if n_consensus_3plus > 0 else "No molecule achieves 3+ list consensus."}</p>
</div>

<div class="key-finding">
<h3>Interpretation Guide</h3>
<p><b>For purchase decisions:</b> Prioritize consensus hits (2+ lists). Between two consensus hits,
prefer the one with better Vina score (more balanced scoring function).<br>
<b>For pharmacophore:</b> Use the top consensus hit as template — its binding pose is
validated by two independent scoring functions.<br>
<b>For reference_docking:</b> Run the top 15-20 consensus hits through GB/SA Hawkins
to resolve remaining divergence with solvated scoring.</p>
</div>

<div class="footer">
    <p>Generated by molecular_docking — Module 07a (Cross-Engine Comparison) v1.0</p>
</div>

</body>
</html>"""
    
    return html


# =============================================================================
# PIPELINE
# =============================================================================

def run_cross_engine_comparison(
        gnina_scores_csv: str,
        dock6_scores_csv: str,
        output_dir: str,
        campaign_id: str = "",
        gnina_composite_csv: Optional[str] = None,
        dock6_all_poses_csv: Optional[str] = None,
        dock6_composite_csv: Optional[str] = None,
        top_n: int = 20,
) -> Dict[str, Any]:
    """
    Run cross-engine comparison: GNINA vs DOCK6.
    
    Reads: 02c gnina_scores, 01e dock6_scores, 05g/04e composite rankings
    Saves: cross_engine_report.html, cross_engine_ranking.csv, consensus_hits.csv
    """
    logger.info("=" * 60)
    logger.info("07a CROSS-ENGINE COMPARISON v1.0")
    logger.info("=" * 60)
    
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    fig_path = out_path / "figures"
    fig_path.mkdir(exist_ok=True)
    
    # Load data
    logger.info("  Loading GNINA data...")
    df_gnina = _load_gnina_data(gnina_scores_csv, gnina_composite_csv)
    logger.info(f"    {len(df_gnina)} molecules")
    
    logger.info("  Loading DOCK6 data...")
    df_dock6 = _load_dock6_data(dock6_scores_csv, dock6_all_poses_csv, dock6_composite_csv)
    logger.info(f"    {len(df_dock6)} molecules")
    
    # Merge
    logger.info("  Merging engines...")
    df_merged = merge_engines(df_gnina, df_dock6)
    logger.info(f"    {len(df_merged)} unique molecules")
    
    n_both = len(df_merged.dropna(subset=['vina_affinity', 'best_Grid_Score']))
    logger.info(f"    {n_both} in both engines")
    
    # Compute consensus
    logger.info(f"  Computing consensus (top {top_n})...")
    df_consensus = compute_consensus(df_merged, top_n=top_n)
    
    n_2plus = len(df_consensus[df_consensus['n_lists'] >= 2])
    n_3plus = len(df_consensus[df_consensus['n_lists'] >= 3])
    n_4 = len(df_consensus[df_consensus['n_lists'] == 4])
    logger.info(f"    Consensus 2+: {n_2plus}, 3+: {n_3plus}, 4: {n_4}")
    
    # Save CSV outputs
    ranking_csv = out_path / "cross_engine_ranking.csv"
    df_consensus.to_csv(ranking_csv, index=False)
    logger.info(f"  Saved: {ranking_csv}")
    
    consensus_csv = out_path / "consensus_hits.csv"
    df_consensus[df_consensus['n_lists'] >= 2].to_csv(consensus_csv, index=False)
    logger.info(f"  Saved: {consensus_csv}")
    
    # Generate figures
    logger.info("  Generating figures...")
    figures = {}
    if HAS_MPL:
        figures["scatter"] = fig_vina_vs_gridscore(df_consensus)
        figures["ranks"] = fig_rank_comparison(df_consensus)
        figures["overlap"] = fig_overlap_venn(df_consensus, top_n)
        
        for name, b64 in figures.items():
            if b64:
                img_data = base64.b64decode(b64)
                with open(fig_path / f"{name}.png", "wb") as f:
                    f.write(img_data)
    
    # Generate HTML
    logger.info("  Generating HTML report...")
    html = generate_html_report(
        df_consensus=df_consensus,
        figures=figures,
        campaign_id=campaign_id,
        top_n=top_n,
        n_gnina=len(df_gnina),
        n_dock6=len(df_dock6),
    )
    
    report_path = out_path / "cross_engine_report.html"
    report_path.write_text(html, encoding='utf-8')
    
    logger.info(f"\n{'=' * 60}")
    logger.info("CROSS-ENGINE COMPARISON COMPLETE")
    logger.info(f"{'=' * 60}")
    logger.info(f"  Report: {report_path}")
    logger.info(f"  Molecules: {len(df_consensus)} total, {n_both} in both")
    logger.info(f"  Consensus: {n_2plus} (2+), {n_3plus} (3+), {n_4} (4)")
    
    # Log top consensus
    top5 = df_consensus[df_consensus['n_lists'] >= 2].head(5)
    for _, row in top5.iterrows():
        name = row['Name']
        n = int(row['n_lists'])
        vina = f"Vina={row['vina_affinity']:.2f}" if pd.notna(row.get('vina_affinity')) else ""
        gs = f"MeanGS={row['mean_Grid_Score']:.1f}" if pd.notna(row.get('mean_Grid_Score')) else ""
        logger.info(f"    #{int(row['consensus_rank'])}: {name} ({n} lists, {vina}, {gs})")
    
    return {
        "success": True,
        "report_path": str(report_path),
        "n_molecules": len(df_consensus),
        "n_consensus_2plus": n_2plus,
        "n_consensus_3plus": n_3plus,
    }

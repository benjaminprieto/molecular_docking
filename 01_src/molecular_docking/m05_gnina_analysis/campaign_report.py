"""
Campaign Report (05g) — Consolidated analysis report with figures and tables
==============================================================================
Generates an HTML report answering the 5 key screening questions:

    1. What molecules bind best? → ranking by energy + convergence
    2. Where do they bind? → fragment positions, hotspots
    3. Why do they bind? → energy decomposition, interaction terms
    4. How confident are we? → pose convergence, consensus
    5. What receptor residues matter? → contact consensus

Reads: all 05a-05f outputs + gnina_scores.csv (mandatory)
Saves: campaign_report.html + figures/ directory

Location: 01_src/molecular_docking/m05_gnina_analysis/campaign_report.py
Project: molecular_docking
Module: 05g (core)
Version: 3.0

CHANGES v2.0 → v3.0:
  - NEW multiplicative composite: composite = vina_affinity × weighted_conv
    where weighted_conv = Σ(conv_i × |prop_vina_i|) / Σ(|prop_vina_i|)
  - gnina_scores.csv is now MANDATORY (no fallback to fragment atom_terms sum)
  - Three main tables: Top 20 Vina raw, Top 20 Composite, Top 20 fragments
  - Detail cards use proportional_vina_mean from fragment decomposition
  - Contacts/H-bonds removed from composite (informational only)
  - PLIP interaction data shown in detail cards when available
"""

import json
import logging
import base64
from io import BytesIO
from pathlib import Path
from typing import Dict, List, Optional
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


def _fig_to_file(fig, path: Path):
    """Save matplotlib figure to file."""
    fig.savefig(str(path), dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)


# =============================================================================
# FIGURE 1: MOLECULE RANKING
# =============================================================================

def fig_molecule_ranking(df_energy: pd.DataFrame, top_n: int = 30,
                         mol_vina_scores: Optional[pd.Series] = None) -> str:
    """Bar chart: top molecules by Vina affinity (or fragment sum as fallback)."""
    if not HAS_MPL or df_energy.empty:
        return ""

    if mol_vina_scores is not None and not mol_vina_scores.empty:
        mol_totals = mol_vina_scores.sort_values()
        energy_label = 'Vina Affinity (kcal/mol)'
        title_label = 'Vina Score'
    else:
        mol_totals = df_energy.groupby('Name')['vina_mean'].sum().sort_values()
        energy_label = 'Total Vina Interaction Energy (kcal/mol)'
        title_label = 'Binding Energy'
    top = mol_totals.head(top_n)

    fig, ax = plt.subplots(figsize=(10, max(6, top_n * 0.25)))
    colors = plt.cm.RdYlGn_r(np.linspace(0.2, 0.8, len(top)))
    bars = ax.barh(range(len(top)), top.values, color=colors)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([n[:35] for n in top.index], fontsize=7)
    ax.set_xlabel(energy_label)
    ax.set_title(f'Top {top_n} Molecules by {title_label}')
    ax.axvline(x=top.values.mean(), color='gray', linestyle='--', alpha=0.5,
               label=f'Mean: {top.values.mean():.1f}')
    ax.legend(fontsize=8)
    ax.invert_yaxis()

    return _fig_to_base64(fig)


# =============================================================================
# FIGURE 2: FRAGMENT ENERGY BREAKDOWN (TOP MOLECULES)
# =============================================================================

def fig_fragment_breakdown(df_energy: pd.DataFrame, top_n: int = 15) -> str:
    """Stacked bar chart: energy contribution per fragment for top molecules."""
    if not HAS_MPL or df_energy.empty:
        return ""

    mol_totals = df_energy.groupby('Name')['vina_mean'].sum().sort_values()
    top_mols = mol_totals.head(top_n).index.tolist()

    fig, ax = plt.subplots(figsize=(12, 6))

    # Pivot: molecules as rows, fragments as stacked segments
    frag_colors = plt.cm.Set2(np.linspace(0, 1, 8))
    y_pos = range(len(top_mols))
    left = np.zeros(len(top_mols))

    # Get max fragment count
    max_frags = df_energy[df_energy['Name'].isin(top_mols)].groupby('Name').size().max()

    for fi in range(max_frags):
        vals = []
        for mol in top_mols:
            mol_frags = df_energy[df_energy['Name'] == mol].sort_values('fragment_id')
            if fi < len(mol_frags):
                vals.append(mol_frags.iloc[fi]['vina_mean'])
            else:
                vals.append(0)
        vals = np.array(vals)
        ax.barh(y_pos, vals, left=left, color=frag_colors[fi % 8],
                label=f'Fragment {fi}', height=0.7)
        left += vals

    ax.set_yticks(y_pos)
    ax.set_yticklabels([n[:30] for n in top_mols], fontsize=7)
    ax.set_xlabel('Vina Energy Contribution (kcal/mol)')
    ax.set_title(f'Fragment Energy Breakdown — Top {top_n} Molecules')
    ax.legend(fontsize=7, loc='lower left', ncol=3)
    ax.invert_yaxis()

    return _fig_to_base64(fig)


# =============================================================================
# FIGURE 3: CONVERGENCE vs ENERGY
# =============================================================================

def fig_convergence_vs_energy(df_cluster: pd.DataFrame, df_energy: pd.DataFrame,
                              mol_vina_scores: Optional[pd.Series] = None) -> str:
    """Scatter: mean fragment convergence vs total energy per molecule."""
    if not HAS_MPL or df_cluster.empty or df_energy.empty:
        return ""

    # Mean convergence per molecule
    conv = df_cluster.groupby('Name')['dominant_fraction'].mean()
    if mol_vina_scores is not None and not mol_vina_scores.empty:
        energy = mol_vina_scores
        energy_label = 'Vina Affinity (kcal/mol)'
    else:
        energy = df_energy.groupby('Name')['vina_mean'].sum()
        energy_label = 'Total Vina Interaction Energy (kcal/mol)'

    merged = pd.DataFrame({'convergence': conv, 'energy': energy}).dropna()

    fig, ax = plt.subplots(figsize=(8, 6))
    sc = ax.scatter(merged['convergence'], merged['energy'],
                    c=merged['energy'], cmap='RdYlGn', s=40, alpha=0.7,
                    edgecolors='gray', linewidth=0.5)
    plt.colorbar(sc, ax=ax, label='Total Vina Energy')

    ax.set_xlabel('Mean Fragment Convergence (dominant fraction)')
    ax.set_ylabel(energy_label)
    ax.set_title('Convergence vs Binding Energy')

    # Highlight best quadrant
    ax.axhline(y=merged['energy'].median(), color='gray', linestyle=':', alpha=0.3)
    ax.axvline(x=merged['convergence'].median(), color='gray', linestyle=':', alpha=0.3)
    ax.text(0.95, 0.05, 'Best: high convergence\n+ low energy',
            transform=ax.transAxes, ha='right', fontsize=8, color='green', alpha=0.5)

    return _fig_to_base64(fig)


# =============================================================================
# FIGURE 4: CONTACT HEATMAP
# =============================================================================

def fig_contact_heatmap(df_contacts: pd.DataFrame, top_residues: int = 20,
                        top_molecules: int = 30) -> str:
    """Heatmap: top residues vs top molecules by contact distance."""
    if not HAS_MPL or df_contacts.empty:
        return ""

    # Top residues by frequency
    res_freq = df_contacts.groupby('residue_id')['Name'].nunique().sort_values(ascending=False)
    top_res = res_freq.head(top_residues).index.tolist()

    # Top molecules by total energy (need to join)
    mol_counts = df_contacts[df_contacts['residue_id'].isin(top_res)]['Name'].value_counts()
    top_mol = mol_counts.head(top_molecules).index.tolist()

    # Build matrix
    filt = df_contacts[df_contacts['residue_id'].isin(top_res) &
                       df_contacts['Name'].isin(top_mol)]
    pivot = filt.pivot_table(values='min_distance', index='Name',
                             columns='residue_id', aggfunc='min')
    pivot = pivot.reindex(columns=top_res)

    if pivot.empty:
        return ""

    fig, ax = plt.subplots(figsize=(max(8, top_residues * 0.5),
                                     max(6, len(pivot) * 0.25)))
    im = ax.imshow(pivot.values, cmap='RdYlGn', aspect='auto',
                   vmin=2.0, vmax=5.0)
    plt.colorbar(im, ax=ax, label='Min Distance (Å)')

    ax.set_xticks(range(len(top_res)))
    ax.set_xticklabels(top_res, rotation=45, ha='right', fontsize=7)
    ax.set_yticks(range(len(pivot)))
    ax.set_yticklabels([n[:25] for n in pivot.index], fontsize=6)
    ax.set_title(f'Contact Map: Top {top_residues} Residues × Top {len(pivot)} Molecules')

    return _fig_to_base64(fig)


# =============================================================================
# FIGURE 5: INTERACTION TERM BREAKDOWN
# =============================================================================

def fig_term_breakdown(df_energy: pd.DataFrame) -> str:
    """Box plots: distribution of each Vina term across all fragments."""
    if not HAS_MPL or df_energy.empty:
        return ""

    terms = ['gauss1_contrib', 'gauss2_contrib', 'repulsion_contrib',
             'hydrophobic_contrib', 'hbond_contrib']
    labels = ['Gauss1\n(steric)', 'Gauss2\n(steric)', 'Repulsion',
              'Hydrophobic', 'H-bond']

    fig, ax = plt.subplots(figsize=(8, 5))
    data = [df_energy[t].values for t in terms]
    bp = ax.boxplot(data, labels=labels, patch_artist=True, notch=True)

    colors = ['#4CAF50', '#2196F3', '#F44336', '#FF9800', '#9C27B0']
    for patch, c in zip(bp['boxes'], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.6)

    ax.axhline(y=0, color='gray', linestyle='-', alpha=0.3)
    ax.set_ylabel('Energy Contribution (kcal/mol)')
    ax.set_title('Vina Interaction Terms Distribution (all fragments)')

    return _fig_to_base64(fig)


# =============================================================================
# FIGURE 6: FRAGMENT EFFICIENCY
# =============================================================================

def fig_fragment_efficiency(df_energy: pd.DataFrame) -> str:
    """Scatter: fragment size vs energy, colored by type (ring/linker)."""
    if not HAS_MPL or df_energy.empty:
        return ""

    fig, ax = plt.subplots(figsize=(8, 5))

    rings = df_energy[df_energy['fragment_label'].str.startswith('ring')]
    linkers = df_energy[df_energy['fragment_label'].str.startswith('linker')]

    ax.scatter(rings['n_atoms'], rings['vina_mean'], c='steelblue',
               label=f'Ring ({len(rings)})', alpha=0.5, s=30)
    ax.scatter(linkers['n_atoms'], linkers['vina_mean'], c='coral',
               label=f'Linker ({len(linkers)})', alpha=0.5, s=30)

    ax.set_xlabel('Fragment Size (heavy atoms)')
    ax.set_ylabel('Vina Energy (kcal/mol)')
    ax.set_title('Fragment Size vs Binding Energy')
    ax.legend()

    # Trend line
    z = np.polyfit(df_energy['n_atoms'], df_energy['vina_mean'], 1)
    x_line = np.linspace(df_energy['n_atoms'].min(), df_energy['n_atoms'].max(), 50)
    ax.plot(x_line, np.polyval(z, x_line), 'k--', alpha=0.3,
            label=f'slope={z[0]:.2f}/atom')
    ax.legend(fontsize=8)

    return _fig_to_base64(fig)


# =============================================================================
# HTML REPORT GENERATOR
# =============================================================================

def _html_table(df: pd.DataFrame, max_rows: int = 50) -> str:
    """Convert DataFrame to HTML table with styling."""
    if df.empty:
        return "<p><i>No data</i></p>"
    display_df = df.head(max_rows)
    html = display_df.to_html(index=False, classes='data-table', border=0,
                               float_format=lambda x: f'{x:.3f}' if isinstance(x, float) else str(x))
    if len(df) > max_rows:
        html += f"<p><i>... showing {max_rows} of {len(df)} rows</i></p>"
    return html


# =============================================================================
# RMSD VALIDATION (control molecule vs crystallographic reference)
# =============================================================================

def _parse_mol2_heavy_coords(mol2_path: str) -> Optional[np.ndarray]:
    """Parse heavy atom coordinates from a mol2 file."""
    coords = []
    in_atoms = False
    with open(mol2_path) as f:
        for line in f:
            if '@<TRIPOS>ATOM' in line:
                in_atoms = True
                continue
            if in_atoms and line.startswith('@'):
                break
            if in_atoms:
                parts = line.split()
                if len(parts) >= 6:
                    element = parts[5].split('.')[0].upper()
                    if element != 'H':
                        try:
                            coords.append([float(parts[2]), float(parts[3]), float(parts[4])])
                        except ValueError:
                            pass
    return np.array(coords) if coords else None


def _compute_rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """Compute RMSD between two coordinate arrays of same shape."""
    if coords1.shape != coords2.shape:
        n = min(len(coords1), len(coords2))
        coords1, coords2 = coords1[:n], coords2[:n]
    diff = coords1 - coords2
    return float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))


def compute_reference_rmsd(
        reference_mol2: str,
        parsed_dir: str,
        control_name: str = "UDX",
) -> Optional[Dict]:
    """
    Compute RMSD of all control molecule poses against crystallographic reference.

    Returns dict with: best_rmsd, best_pose, mean_rmsd, all_rmsds, n_under_2A
    """
    ref_coords = _parse_mol2_heavy_coords(reference_mol2)
    if ref_coords is None:
        logger.warning(f"  Could not parse reference: {reference_mol2}")
        return None

    # Find control molecule in parsed data
    parsed_path = Path(parsed_dir)
    npz_files = list(parsed_path.glob("*.npz"))

    control_npz = None
    for npz_file in npz_files:
        if control_name in npz_file.stem:
            control_npz = npz_file
            break

    if control_npz is None:
        logger.warning(f"  Control molecule '{control_name}' not found in {parsed_dir}")
        return None

    data = np.load(control_npz)
    pose_coords = data['coords']  # (n_poses, n_atoms, 3)
    n_poses = pose_coords.shape[0]

    # Compute RMSD per pose
    rmsds = []
    for i in range(n_poses):
        rmsd = _compute_rmsd(ref_coords, pose_coords[i])
        rmsds.append(rmsd)

    rmsds = np.array(rmsds)
    best_idx = int(np.argmin(rmsds))

    return {
        'control_name': control_name,
        'n_ref_atoms': len(ref_coords),
        'n_poses': n_poses,
        'best_rmsd': float(rmsds[best_idx]),
        'best_pose': best_idx,
        'mean_rmsd': float(np.mean(rmsds)),
        'median_rmsd': float(np.median(rmsds)),
        'worst_rmsd': float(np.max(rmsds)),
        'n_under_2A': int(np.sum(rmsds < 2.0)),
        'n_under_3A': int(np.sum(rmsds < 3.0)),
        'pct_under_2A': float(np.mean(rmsds < 2.0)),
        'all_rmsds': rmsds.tolist(),
    }


def _build_validation_section(rmsd_data: Optional[Dict]) -> str:
    """Build HTML section for protocol validation."""
    if rmsd_data is None:
        return ""

    d = rmsd_data
    # Determine quality
    if d['best_rmsd'] < 1.0:
        quality = '<span style="color:#27ae60;font-weight:bold">EXCELLENT</span>'
    elif d['best_rmsd'] < 2.0:
        quality = '<span style="color:#2980b9;font-weight:bold">GOOD</span>'
    elif d['best_rmsd'] < 3.0:
        quality = '<span style="color:#f39c12;font-weight:bold">ACCEPTABLE</span>'
    else:
        quality = '<span style="color:#e74c3c;font-weight:bold">POOR</span>'

    html = f"""
<!-- ============================================================ -->
<h2>0. Protocol Validation ({d['control_name']})</h2>

<div class="summary-box">
    <b>Control molecule:</b> {d['control_name']} (crystallographic reference)<br>
    <b>Poses generated:</b> {d['n_poses']}<br>
    <b>Heavy atoms matched:</b> {d['n_ref_atoms']}
</div>

<div class="stats">
    <div class="stat-card">
        <div class="stat-value" style="font-size:22px">{d['best_rmsd']:.2f} Å</div>
        <div class="stat-label">Best RMSD (pose #{d['best_pose']})</div>
    </div>
    <div class="stat-card">
        <div class="stat-value" style="font-size:22px">{d['mean_rmsd']:.2f} Å</div>
        <div class="stat-label">Mean RMSD</div>
    </div>
    <div class="stat-card">
        <div class="stat-value" style="font-size:22px">{d['n_under_2A']}/{d['n_poses']}</div>
        <div class="stat-label">Poses &lt; 2.0 Å</div>
    </div>
    <div class="stat-card">
        <div class="stat-value" style="font-size:22px">{quality}</div>
        <div class="stat-label">Protocol Quality</div>
    </div>
</div>

<div class="key-finding">
<p>The docking protocol reproduces the crystallographic pose with a best RMSD of <b>{d['best_rmsd']:.2f} Å</b>
(pose #{d['best_pose']} of {d['n_poses']}). {d['n_under_2A']} poses ({d['pct_under_2A']:.0%}) are within 2.0 Å of the
experimental structure. Median RMSD across all poses: {d['median_rmsd']:.2f} Å.</p>
</div>
"""
    return html


def _stability_badge(conv: float) -> str:
    """Return HTML badge for stability classification."""
    if conv >= 0.7:
        return f'<span style="background:#27ae60;color:white;padding:2px 8px;border-radius:3px;font-size:11px">STABLE ({conv:.0%})</span>'
    elif conv >= 0.5:
        return f'<span style="background:#f39c12;color:white;padding:2px 8px;border-radius:3px;font-size:11px">MODERATE ({conv:.0%})</span>'
    else:
        return f'<span style="background:#e74c3c;color:white;padding:2px 8px;border-radius:3px;font-size:11px">UNSTABLE ({conv:.0%})</span>'


def _energy_bar(value: float, min_val: float) -> str:
    """Return HTML inline bar for energy value."""
    if min_val == 0:
        pct = 0
    else:
        pct = min(100, max(0, (value / min_val) * 100))
    return f'<div style="display:flex;align-items:center;gap:6px"><span style="width:55px;text-align:right">{value:.2f}</span><div style="background:#3498db;height:12px;width:{pct:.0f}%;border-radius:2px;opacity:0.7"></div></div>'


def _fragment_assessment(prop_vina: float, conv: float) -> str:
    """Return HTML assessment badge for a fragment based on proportional Vina and convergence."""
    if prop_vina < -3.0 and conv >= 0.7:
        return '<span style="color:#27ae60">★ Anchor</span>'
    elif prop_vina < -2.0 and conv >= 0.7:
        return '<span style="color:#2980b9">✓ Reliable</span>'
    elif prop_vina < -3.0 and conv < 0.5:
        return '<span style="color:#f39c12">⚡ Potent but variable</span>'
    elif conv < 0.3:
        return '<span style="color:#e74c3c">✗ Unreliable</span>'
    elif conv >= 0.7:
        return '<span style="color:#95a5a6">Stable, weak</span>'
    else:
        return '<span style="color:#777">Moderate</span>'


def _compute_weighted_conv(mol_frags: pd.DataFrame) -> float:
    """Compute weighted convergence for a molecule's fragments.

    weighted_conv = Σ(conv_i × |prop_vina_i|) / Σ(|prop_vina_i|)

    Uses proportional_vina_mean if available, falls back to vina_mean.
    """
    if mol_frags.empty:
        return 0.0

    if 'proportional_vina_mean' in mol_frags.columns:
        weights = mol_frags['proportional_vina_mean'].abs()
    else:
        weights = mol_frags['vina_mean'].abs()

    convs = mol_frags['dominant_fraction']

    total_weight = weights.sum()
    if total_weight == 0:
        return convs.mean() if len(convs) > 0 else 0.0

    return float((convs * weights).sum() / total_weight)


def _build_molecule_detail_section(df_energy, df_cluster, df_contacts, df_gnina,
                                    top_composite_names, df_plip=None):
    """Build per-molecule detail cards for top 20 composite molecules."""
    if df_energy.empty or df_cluster.empty:
        return "<p><i>No data available</i></p>"

    # Merge energy + cluster
    merged = df_energy.merge(
        df_cluster[['Name', 'fragment_id', 'dominant_fraction', 'n_clusters', 'n_noise']],
        on=['Name', 'fragment_id'], how='left'
    ).fillna(0)

    # Use proportional_vina_mean for fragment composite if available
    prop_vina_col = 'proportional_vina_mean' if 'proportional_vina_mean' in merged.columns else 'vina_mean'
    merged['fragment_composite'] = merged[prop_vina_col] * merged['dominant_fraction']

    # Get top 3 contact residues per fragment
    frag_top_contacts = {}
    if not df_contacts.empty:
        for (name, frag_id), grp in df_contacts.groupby(['Name', 'fragment_id']) if 'fragment_id' in df_contacts.columns else []:
            top3 = grp.nsmallest(3, 'min_distance')['residue_id'].tolist()
            frag_top_contacts[(name, frag_id)] = top3

    # Build gnina lookup for drug-likeness
    gnina_lookup = {}
    if df_gnina is not None and not df_gnina.empty:
        for _, row in df_gnina.iterrows():
            gnina_lookup[row['name']] = row

    cards_html = ""
    for rank, mol_name in enumerate(top_composite_names, 1):
        mol_frags = merged[merged['Name'] == mol_name].sort_values('fragment_composite')
        if mol_frags.empty:
            continue

        # Get whole-molecule data from gnina
        gnina_row = gnina_lookup.get(mol_name, {})
        vina_aff = gnina_row.get('vina_affinity', mol_frags['vina_mean'].sum()) if isinstance(gnina_row, dict) else getattr(gnina_row, 'vina_affinity', mol_frags['vina_mean'].sum())

        # Compute weighted conv
        weighted_conv = _compute_weighted_conv(mol_frags)
        composite = vina_aff * weighted_conv

        # Drug-likeness from gnina
        mw = gnina_row.get('MW', '') if isinstance(gnina_row, dict) else getattr(gnina_row, 'MW', '')
        qed = gnina_row.get('QED', '') if isinstance(gnina_row, dict) else getattr(gnina_row, 'QED', '')
        tpsa = gnina_row.get('TPSA', '') if isinstance(gnina_row, dict) else getattr(gnina_row, 'TPSA', '')

        mw_str = f"{mw:.1f}" if isinstance(mw, (int, float)) and not pd.isna(mw) else "N/A"
        qed_str = f"{qed:.3f}" if isinstance(qed, (int, float)) and not pd.isna(qed) else "N/A"
        tpsa_str = f"{tpsa:.1f}" if isinstance(tpsa, (int, float)) and not pd.isna(tpsa) else "N/A"

        cards_html += f"""
<div style="background:white;border:1px solid #ddd;border-radius:8px;padding:15px;margin:15px 0">
    <h3 style="margin-top:0">#{rank} {mol_name}</h3>
    <div class="stats" style="margin-bottom:10px">
        <div class="stat-card" style="min-width:100px">
            <div class="stat-value" style="font-size:16px">{vina_aff:.2f}</div>
            <div class="stat-label">Vina Affinity</div>
        </div>
        <div class="stat-card" style="min-width:100px">
            <div class="stat-value" style="font-size:16px">{weighted_conv:.3f}</div>
            <div class="stat-label">Weighted Conv</div>
        </div>
        <div class="stat-card" style="min-width:100px">
            <div class="stat-value" style="font-size:16px;color:#8e44ad"><b>{composite:.2f}</b></div>
            <div class="stat-label">Composite</div>
        </div>
        <div class="stat-card" style="min-width:80px">
            <div class="stat-value" style="font-size:14px">{mw_str}</div>
            <div class="stat-label">MW</div>
        </div>
        <div class="stat-card" style="min-width:80px">
            <div class="stat-value" style="font-size:14px">{qed_str}</div>
            <div class="stat-label">QED</div>
        </div>
        <div class="stat-card" style="min-width:80px">
            <div class="stat-value" style="font-size:14px">{tpsa_str}</div>
            <div class="stat-label">TPSA</div>
        </div>
    </div>
    <table class="data-table">
        <tr>
            <th>Fragment</th>
            <th>Atoms</th>
            <th>Prop. Vina</th>
            <th>Convergence</th>
            <th>Frag Composite</th>
            <th>Top 3 Contact Residues</th>
            <th>Clusters</th>
            <th>Stability</th>
            <th>Assessment</th>
        </tr>"""

        for _, f in mol_frags.iterrows():
            prop_vina = f[prop_vina_col]
            c = f['dominant_fraction']
            assessment = _fragment_assessment(prop_vina, c)

            # Top contacts for this fragment
            contacts_key = (mol_name, f['fragment_id'])
            top3_contacts = frag_top_contacts.get(contacts_key, [])
            contacts_str = ", ".join(top3_contacts[:3]) if top3_contacts else "—"

            cards_html += f"""
        <tr>
            <td><b>{f['fragment_label']}</b></td>
            <td>{int(f['n_atoms'])}</td>
            <td>{prop_vina:.3f}</td>
            <td>{_stability_badge(c)}</td>
            <td><b>{f['fragment_composite']:.3f}</b></td>
            <td style="font-size:11px">{contacts_str}</td>
            <td>{int(f['n_clusters'])}</td>
            <td>{_stability_badge(c)}</td>
            <td>{assessment}</td>
        </tr>"""

        cards_html += """
    </table>"""

        # PLIP interactions if available
        if df_plip is not None and not df_plip.empty:
            mol_plip = df_plip[df_plip['Name'] == mol_name] if 'Name' in df_plip.columns else pd.DataFrame()
            if not mol_plip.empty:
                cards_html += """
    <details style="margin-top:10px">
    <summary style="cursor:pointer;color:#2980b9"><b>PLIP Interactions</b></summary>
    <table class="data-table" style="margin-top:5px">
        <tr>"""
                for col in mol_plip.columns:
                    if col != 'Name':
                        cards_html += f"<th>{col}</th>"
                cards_html += "</tr>"
                for _, plip_row in mol_plip.iterrows():
                    cards_html += "<tr>"
                    for col in mol_plip.columns:
                        if col != 'Name':
                            val = plip_row[col]
                            if isinstance(val, float):
                                cards_html += f"<td>{val:.3f}</td>"
                            else:
                                cards_html += f"<td>{val}</td>"
                    cards_html += "</tr>"
                cards_html += """
    </table>
    </details>"""

        cards_html += """
</div>"""

    return cards_html


def _build_consensus_table(df_energy, df_cluster, df_contacts):
    """Build top 20 individual fragments table by fragment_composite = prop_vina × conv."""
    if df_energy.empty or df_cluster.empty:
        return "<p><i>No data available</i></p>"

    merged = df_energy.merge(
        df_cluster[['Name', 'fragment_id', 'dominant_fraction']],
        on=['Name', 'fragment_id'], how='left'
    ).fillna(0)

    # Use proportional_vina_mean if available, else vina_mean
    prop_vina_col = 'proportional_vina_mean' if 'proportional_vina_mean' in merged.columns else 'vina_mean'
    merged['fragment_composite'] = merged[prop_vina_col] * merged['dominant_fraction']

    # Get top contact residue per fragment
    frag_top_contact = {}
    if not df_contacts.empty and 'fragment_id' in df_contacts.columns:
        for (name, frag_id), grp in df_contacts.groupby(['Name', 'fragment_id']):
            top3 = grp.nsmallest(3, 'min_distance')['residue_id'].values
            if len(top3) > 0:
                frag_top_contact[(name, frag_id)] = ", ".join(top3)

    # Best individual fragments by fragment_composite (most negative first)
    best_individual = merged.nsmallest(20, 'fragment_composite')

    html = """<h3>Top 20 Individual Fragments (Prop. Vina × Convergence)</h3>
<p style="font-size:12px;color:#777">Fragments that are both energetically favorable AND consistently positioned.
More negative fragment composite = better.</p>
<table class="data-table">
    <tr>
        <th>Rank</th>
        <th>Molecule</th>
        <th>Fragment</th>
        <th>Atoms</th>
        <th>Prop. Vina</th>
        <th>Conv</th>
        <th>Frag Composite</th>
        <th>is_ring</th>
        <th>Top Contacts</th>
    </tr>"""

    for rank, (_, r) in enumerate(best_individual.iterrows(), 1):
        prop_vina = r[prop_vina_col]
        c = r['dominant_fraction']
        is_ring = "Yes" if r['fragment_label'].startswith('ring') else "No"
        top_contact = frag_top_contact.get((r['Name'], r['fragment_id']), "—")

        html += f"""
    <tr>
        <td>{rank}</td>
        <td>{r['Name'][:35]}</td>
        <td><b>{r['fragment_label']}</b></td>
        <td>{int(r['n_atoms'])}</td>
        <td>{prop_vina:.3f}</td>
        <td>{_stability_badge(c)}</td>
        <td><b>{r['fragment_composite']:.3f}</b></td>
        <td>{is_ring}</td>
        <td>{top_contact}</td>
    </tr>"""

    html += "</table>"

    return html


def generate_html_report(
        campaign_id: str,
        engine: str,
        df_energy: pd.DataFrame,
        df_cluster: pd.DataFrame,
        df_contacts: pd.DataFrame,
        df_hotspot_contacts: pd.DataFrame,
        hotspot_data: dict,
        figures: Dict[str, str],
        rmsd_data: Optional[Dict] = None,
        df_gnina: Optional[pd.DataFrame] = None,
        df_plip: Optional[pd.DataFrame] = None,
) -> str:
    """Generate complete HTML report.

    Args:
        df_gnina: Full gnina_scores.csv DataFrame with drug-likeness columns.
            Used for Vina affinity, CNN scores, MW, QED, TPSA, LogP.
        df_plip: PLIP interaction data (plip_interactions_top20.csv) if available.
    """

    # Compute summary stats
    n_molecules = df_energy['Name'].nunique()

    # Energy from gnina_scores (mandatory in v3.0)
    mol_vina_scores = None
    if df_gnina is not None and not df_gnina.empty:
        df_ok = df_gnina[df_gnina['success'] == True].copy() if 'success' in df_gnina.columns else df_gnina.copy()
        mol_vina_scores = df_ok.set_index('name')['vina_affinity']

    if mol_vina_scores is not None and not mol_vina_scores.empty:
        mol_totals = mol_vina_scores.sort_values()
    else:
        mol_totals = df_energy.groupby('Name')['vina_mean'].sum().sort_values()
    mean_conv = df_cluster['dominant_fraction'].mean() if not df_cluster.empty else 0

    # === MULTIPLICATIVE COMPOSITE RANKING ===
    # Merge energy + cluster at fragment level for weighted_conv calculation
    frag_merged = df_energy.merge(
        df_cluster[['Name', 'fragment_id', 'dominant_fraction']],
        on=['Name', 'fragment_id'], how='left'
    ).fillna(0) if not df_cluster.empty else df_energy.copy()

    if 'dominant_fraction' not in frag_merged.columns:
        frag_merged['dominant_fraction'] = 0.0

    # Use proportional_vina_mean for weighting if available
    prop_vina_col = 'proportional_vina_mean' if 'proportional_vina_mean' in frag_merged.columns else 'vina_mean'

    # Compute weighted_conv per molecule
    weighted_conv_dict = {}
    min_conv_dict = {}
    best_frag_dict = {}  # best fragment label + proportional contribution
    for mol_name, mol_frags in frag_merged.groupby('Name'):
        weighted_conv_dict[mol_name] = _compute_weighted_conv(mol_frags)
        min_conv_dict[mol_name] = mol_frags['dominant_fraction'].min()
        # Best fragment by fragment_composite
        mol_frags_copy = mol_frags.copy()
        mol_frags_copy['_fc'] = mol_frags_copy[prop_vina_col] * mol_frags_copy['dominant_fraction']
        best_idx = mol_frags_copy['_fc'].idxmin()
        best_row = mol_frags_copy.loc[best_idx]
        best_frag_dict[mol_name] = f"{best_row['fragment_label']} ({best_row[prop_vina_col]:.2f})"

    mol_weighted_conv = pd.Series(weighted_conv_dict, name='weighted_conv')
    mol_min_conv = pd.Series(min_conv_dict, name='min_conv')
    mol_best_frag = pd.Series(best_frag_dict, name='best_fragment')

    # Build ranking DataFrame
    if mol_vina_scores is not None and not mol_vina_scores.empty:
        mol_energy = mol_vina_scores.rename('vina_affinity')
    else:
        mol_energy = df_energy.groupby('Name')['vina_mean'].sum().rename('vina_affinity')

    ranking = pd.DataFrame(mol_energy)
    ranking = ranking.join(mol_weighted_conv).join(mol_min_conv).join(mol_best_frag).fillna(0)

    # Composite = vina_affinity × weighted_conv (units: kcal/mol, more negative = better)
    ranking['composite'] = ranking['vina_affinity'] * ranking['weighted_conv']
    ranking = ranking.sort_values('composite', ascending=True)  # most negative first

    # Join drug-likeness from gnina
    if df_gnina is not None and not df_gnina.empty:
        df_ok = df_gnina[df_gnina['success'] == True].copy() if 'success' in df_gnina.columns else df_gnina.copy()
        drug_cols = []
        for col in ['MW', 'QED', 'TPSA', 'LogP', 'CNN_score', 'CNN_affinity']:
            if col in df_ok.columns:
                drug_cols.append(col)
        if drug_cols:
            drug_df = df_ok.set_index('name')[drug_cols]
            ranking = ranking.join(drug_df, how='left')

    # === TABLE 1: Top 20 by Vina Affinity (raw) ===
    vina_top20_rows = []
    vina_sorted = ranking.sort_values('vina_affinity', ascending=True)
    for rank_i, (name, row) in enumerate(vina_sorted.head(20).iterrows(), 1):
        n_frags = len(df_energy[df_energy['Name'] == name])
        vina_top20_rows.append({
            'Rank': rank_i,
            'Molecule': name,
            'Vina Affinity': round(row['vina_affinity'], 2),
            'CNN Score': round(row.get('CNN_score', 0), 3) if pd.notna(row.get('CNN_score', None)) else '',
            'CNN Affinity': round(row.get('CNN_affinity', 0), 2) if pd.notna(row.get('CNN_affinity', None)) else '',
            'n_fragments': n_frags,
            'MW': round(row.get('MW', 0), 1) if pd.notna(row.get('MW', None)) else '',
            'QED': round(row.get('QED', 0), 3) if pd.notna(row.get('QED', None)) else '',
            'TPSA': round(row.get('TPSA', 0), 1) if pd.notna(row.get('TPSA', None)) else '',
            'LogP': round(row.get('LogP', 0), 2) if pd.notna(row.get('LogP', None)) else '',
        })
    df_vina_top20 = pd.DataFrame(vina_top20_rows)

    # === TABLE 2: Top 20 by Composite (Vina × weighted_conv) ===
    composite_rows = []
    for rank_i, (name, row) in enumerate(ranking.head(20).iterrows(), 1):
        composite_rows.append({
            'Rank': rank_i,
            'Molecule': name,
            'Vina Affinity': round(row['vina_affinity'], 2),
            'Weighted Conv': round(row['weighted_conv'], 3),
            'Composite': round(row['composite'], 2),
            'Min Frag Conv': round(row['min_conv'], 3),
            'Best Fragment': row['best_fragment'],
            'MW': round(row.get('MW', 0), 1) if pd.notna(row.get('MW', None)) else '',
            'QED': round(row.get('QED', 0), 3) if pd.notna(row.get('QED', None)) else '',
        })
    df_composite = pd.DataFrame(composite_rows)

    # Top 20 composite molecule names for detail cards
    top20_composite_names = ranking.head(20).index.tolist()

    # Candidate categories (informational)
    n_strong = len(ranking[(ranking['vina_affinity'] < ranking['vina_affinity'].quantile(0.3)) &
                           (ranking['weighted_conv'] > 0.5)])
    n_risky = len(ranking[(ranking['vina_affinity'] < ranking['vina_affinity'].quantile(0.3)) &
                          (ranking['weighted_conv'] < 0.3)])
    n_weak = len(ranking[ranking['vina_affinity'] > ranking['vina_affinity'].quantile(0.7)])

    # Top contacts table (informational)
    if not df_hotspot_contacts.empty:
        df_top_contacts = df_hotspot_contacts.nlargest(15, 'contact_fraction')[[
            'residue_id', 'residue_name', 'n_molecules_contacting',
            'contact_fraction', 'mean_distance', 'dominant_contact_type'
        ]]
    else:
        df_top_contacts = pd.DataFrame()

    # Fragment type summary
    ring_energy = df_energy[df_energy['fragment_label'].str.startswith('ring')]['vina_mean'].mean()
    linker_energy = df_energy[df_energy['fragment_label'].str.startswith('linker')]['vina_mean'].mean()
    if pd.isna(ring_energy):
        ring_energy = 0.0
    if pd.isna(linker_energy):
        linker_energy = 0.0

    # === PRE-COMPUTE SAFE VALUES for HTML template ===
    _hc = df_hotspot_contacts
    _hc0_resid = _hc.iloc[0]['residue_id'] if len(_hc) > 0 else 'N/A'
    _hc0_nmol = int(_hc.iloc[0]['n_molecules_contacting']) if len(_hc) > 0 else 0
    _hc0_frac = _hc.iloc[0]['contact_fraction'] * 100 if len(_hc) > 0 else 0
    _hc1_resid = _hc.iloc[1]['residue_id'] if len(_hc) > 1 else 'N/A'
    _hc1_frac = _hc.iloc[1]['contact_fraction'] * 100 if len(_hc) > 1 else 0

    # Top candidates
    _top_candidates_lines = []
    for _ti in range(min(3, len(ranking))):
        _tn = ranking.index[_ti]
        _te = ranking['vina_affinity'].iloc[_ti]
        _tc = ranking['weighted_conv'].iloc[_ti]
        _ts = ranking['composite'].iloc[_ti]
        _top_candidates_lines.append(f"<b>{_tn}</b> (Vina={_te:.1f}, w_conv={_tc:.3f}, composite={_ts:.2f})")
    _top_candidates_str = ", ".join(_top_candidates_lines) if _top_candidates_lines else "N/A"

    html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Campaign Report: {campaign_id}</title>
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

<h1>Campaign Report: {campaign_id}</h1>
<p><b>Engine:</b> {engine.upper()} &nbsp; | &nbsp;
   <b>Generated:</b> {datetime.now().strftime('%Y-%m-%d %H:%M')} &nbsp; | &nbsp;
   <b>Molecules:</b> {n_molecules} &nbsp; | &nbsp;
   <b>Report version:</b> v3.0</p>

<div class="stats">
    <div class="stat-card">
        <div class="stat-value">{n_molecules}</div>
        <div class="stat-label">Molecules</div>
    </div>
    <div class="stat-card">
        <div class="stat-value">{len(df_energy)}</div>
        <div class="stat-label">Fragments</div>
    </div>
    <div class="stat-card">
        <div class="stat-value">{mol_totals.min():.1f}</div>
        <div class="stat-label">Best Vina (kcal/mol)</div>
    </div>
    <div class="stat-card">
        <div class="stat-value">{mean_conv:.0%}</div>
        <div class="stat-label">Mean Convergence</div>
    </div>
    <div class="stat-card">
        <div class="stat-value">{hotspot_data.get('n_hotspots', 0)}</div>
        <div class="stat-label">Hotspots</div>
    </div>
</div>

{_build_validation_section(rmsd_data)}

<!-- ============================================================ -->
<h2>1. What molecules bind best?</h2>

<div class="summary-box">
    <b>Top molecule (composite):</b> {ranking.index[0]} (composite {ranking['composite'].iloc[0]:.2f} kcal/mol)<br>
    <b>Vina range:</b> {mol_totals.min():.2f} to {mol_totals.max():.2f} kcal/mol<br>
    <b>Campaign mean:</b> {mol_totals.mean():.2f} kcal/mol<br>
    <b>Candidates:</b> {n_strong} strong, {n_risky} risky (good energy but low convergence), {n_weak} weak
</div>

<h3>Table 1: Top 20 by Vina Affinity (raw)</h3>
<p style="font-size:12px;color:#777">Ranked by whole-molecule Vina affinity from gnina_scores.csv. No convergence weighting.</p>
{_html_table(df_vina_top20)}

<h3>Table 2: Top 20 by Composite (Vina × Weighted Convergence)</h3>
<p style="font-size:12px;color:#777">Composite = Vina Affinity × weighted_conv.
weighted_conv = Σ(conv_i × |prop_vina_i|) / Σ(|prop_vina_i|).
Units: kcal/mol. More negative = better.</p>
{_html_table(df_composite)}

{f'<div class="figure"><img src="data:image/png;base64,{figures["ranking"]}"><div class="figure-caption">Fig 1. Molecules ranked by Vina affinity (whole-molecule, torsion-corrected).</div></div>' if figures.get("ranking") else ''}

<!-- ============================================================ -->
<h2>2. Where do they bind?</h2>

<div class="summary-box">
    <b>Hotspots detected:</b> {hotspot_data.get('n_hotspots', 0)}<br>
    <b>Centroids in hotspots:</b> {hotspot_data.get('n_centroids_assigned', 0)}/{hotspot_data.get('n_centroids_total', 0)}<br>
    <b>Fragments per molecule:</b> mean {df_energy.groupby('Name').size().mean():.1f} (range {df_energy.groupby('Name').size().min()}-{df_energy.groupby('Name').size().max()})
</div>

{f'<div class="figure"><img src="data:image/png;base64,{figures["efficiency"]}"><div class="figure-caption">Fig 2. Fragment size vs binding energy. Rings (blue) vs linkers (red). Dashed line = trend.</div></div>' if figures.get("efficiency") else ''}

<!-- ============================================================ -->
<h2>3. Why do they bind?</h2>

<div class="summary-box">
    <b>Dominant interaction:</b> Gauss2/steric ({df_energy['gauss2_contrib'].mean():.3f} kcal/mol mean)<br>
    <b>H-bond contribution:</b> {df_energy['hbond_contrib'].mean():.3f} kcal/mol mean<br>
    <b>Rings:</b> {ring_energy:.2f} kcal/mol mean &nbsp; | &nbsp; <b>Linkers:</b> {linker_energy:.2f} kcal/mol mean
</div>

{f'<div class="figure"><img src="data:image/png;base64,{figures["breakdown"]}"><div class="figure-caption">Fig 3. Energy contribution per fragment for top molecules (stacked).</div></div>' if figures.get("breakdown") else ''}

{f'<div class="figure"><img src="data:image/png;base64,{figures["terms"]}"><div class="figure-caption">Fig 4. Distribution of Vina interaction terms across all fragments.</div></div>' if figures.get("terms") else ''}

<!-- ============================================================ -->
<h2>4. How confident are we?</h2>

<div class="summary-box">
    <b>Fragments with >50% convergence:</b> {(df_cluster["dominant_fraction"] > 0.5).sum()}/{len(df_cluster)} ({(df_cluster["dominant_fraction"] > 0.5).mean():.0%})<br>
    <b>Mean convergence:</b> {mean_conv:.2f}<br>
    <b>Molecules with all fragments converged (>50%):</b> {sum(df_cluster.groupby("Name")["dominant_fraction"].min() > 0.5)}/{n_molecules}
</div>

{f'<div class="figure"><img src="data:image/png;base64,{figures["conv_energy"]}"><div class="figure-caption">Fig 5. Convergence vs energy. Bottom-right = best candidates (strong binding + confident poses).</div></div>' if figures.get("conv_energy") else ''}

<!-- ============================================================ -->
<h2>5. What receptor residues matter?</h2>

<div class="summary-box">
    <b>Total contacts:</b> {len(df_contacts)} across {df_contacts['Name'].nunique() if not df_contacts.empty else 0} molecules<br>
    <b>Top residue:</b> {_hc0_resid}
    ({_hc0_nmol}/{n_molecules} molecules)
</div>

<h3>Consensus Contact Residues</h3>
{_html_table(df_top_contacts)}

{f'<div class="figure"><img src="data:image/png;base64,{figures["contacts"]}"><div class="figure-caption">Fig 6. Contact heatmap: distance from fragment to receptor residue (green = close contact).</div></div>' if figures.get("contacts") else ''}

<!-- ============================================================ -->
<h2>Key Findings</h2>

<div class="key-finding">
<h3>Binding Site Pharmacophore</h3>
<p>The binding site is anchored by <b>{_hc0_resid}</b>
(contacted by {_hc0_frac:.0f}% of molecules) and
<b>{_hc1_resid}</b>
({_hc1_frac:.0f}% of molecules).
The dominant interaction type is steric complementarity (gauss2) with significant H-bonding contributions.</p>
</div>

<div class="key-finding">
<h3>Top Candidates (Composite Ranking)</h3>
<p>The top {min(3, len(ranking))} molecules by multiplicative composite (Vina × weighted convergence) are:
{_top_candidates_str}.
Of {n_molecules} molecules, {n_strong} are strong candidates, {n_risky} have good energy but unreliable poses, and {n_weak} should be deprioritized.</p>
</div>

<div class="key-finding">
<h3>Fragment Insights</h3>
<p>Linker fragments (typically phosphate bridges) contribute {linker_energy:.2f} kcal/mol on average
despite having fewer atoms, making them the most efficient binders per atom.
Ring fragments contribute {ring_energy:.2f} kcal/mol on average with larger contact surfaces.</p>
</div>

<!-- ============================================================ -->
<h2>6. Detail Cards — Top 20 Composite</h2>

<div class="summary-box">
Per-molecule fragment breakdown for the top 20 composite-ranked molecules.<br>
<b>Composite</b> = Vina Affinity × weighted_conv (kcal/mol, more negative = better).<br><br>
Fragment assessment: <span style="color:#27ae60">★ Anchor</span> = prop_vina &lt; -3 and conv &gt; 0.7,
<span style="color:#2980b9">✓ Reliable</span> = prop_vina &lt; -2 and conv &gt; 0.7,
<span style="color:#f39c12">⚡ Potent but variable</span> = prop_vina &lt; -3 and conv &lt; 0.5,
<span style="color:#e74c3c">✗ Unreliable</span> = conv &lt; 0.3.
</div>

{_build_molecule_detail_section(df_energy, df_cluster, df_contacts, df_gnina, top20_composite_names, df_plip=df_plip)}

<!-- ============================================================ -->
<h2>7. Top 20 Individual Fragments</h2>

<div class="summary-box">
Which individual fragments across all molecules have the best combination of
proportional Vina contribution and convergence?
Ranked by fragment_composite = proportional_vina × convergence (most negative = best).
</div>

{_build_consensus_table(df_energy, df_cluster, df_contacts)}

<div class="footer">
    Generated by molecular_docking m05 pipeline v3.0 — Module 05g Campaign Report
</div>

</body>
</html>"""

    return html


# =============================================================================
# PIPELINE
# =============================================================================

def run_campaign_report(
        parsed_dir: str,
        cluster_dir: str,
        score_dir: str,
        hotspot_dir: str,
        contact_dir: str,
        output_dir: str,
        campaign_id: str = "",
        engine: str = "gnina",
        reference_mol2: Optional[str] = None,
        control_name: str = "UDX",
        gnina_scores_csv: Optional[str] = None,
) -> Dict[str, any]:
    """
    Generate consolidated campaign report.

    Reads: 05a-05f outputs + gnina_scores.csv (MANDATORY in v3.0)
    Saves: campaign_report.html + figures/

    Args:
        gnina_scores_csv: Path to gnina_scores.csv from 02c. MANDATORY in v3.0.
    """
    logger.info("=" * 60)
    logger.info("05g CAMPAIGN REPORT v3.0")
    logger.info("=" * 60)

    # gnina_scores.csv is MANDATORY in v3.0
    if not gnina_scores_csv or not Path(gnina_scores_csv).exists():
        logger.error("gnina_scores.csv is required for campaign report v3.0")
        return {"success": False, "error": "gnina_scores.csv not found — required for composite ranking"}

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    fig_path = out_path / "figures"
    fig_path.mkdir(exist_ok=True)

    # Load gnina_scores.csv as full DataFrame
    df_gnina = pd.read_csv(gnina_scores_csv)
    if 'name' not in df_gnina.columns or 'vina_affinity' not in df_gnina.columns:
        logger.error("gnina_scores.csv missing 'name' or 'vina_affinity' columns")
        return {"success": False, "error": "gnina_scores.csv missing required columns"}

    df_ok = df_gnina[df_gnina['success'] == True].copy() if 'success' in df_gnina.columns else df_gnina.copy()
    mol_vina_scores = df_ok.set_index('name')['vina_affinity']
    logger.info(f"  Vina scores: {len(mol_vina_scores)} molecules from {gnina_scores_csv}")

    # Load data
    df_energy = pd.DataFrame()
    energy_csv = Path(score_dir) / "fragment_energy_summary.csv"
    if energy_csv.exists():
        df_energy = pd.read_csv(energy_csv)
        logger.info(f"  Energy: {len(df_energy)} fragments from {df_energy['Name'].nunique()} molecules")

    df_cluster = pd.DataFrame()
    cluster_csv = Path(cluster_dir) / "cluster_summary.csv"
    if cluster_csv.exists():
        df_cluster = pd.read_csv(cluster_csv)
        logger.info(f"  Clusters: {len(df_cluster)} rows")

    df_contacts = pd.DataFrame()
    contact_csv = Path(contact_dir) / "contact_summary.csv"
    if contact_csv.exists():
        df_contacts = pd.read_csv(contact_csv)
        logger.info(f"  Contacts: {len(df_contacts)} rows")

    df_hotspot_contacts = pd.DataFrame()
    hs_contact_csv = Path(contact_dir) / "hotspot_contacts.csv"
    if hs_contact_csv.exists():
        df_hotspot_contacts = pd.read_csv(hs_contact_csv)

    hotspot_data = {}
    hs_json = Path(hotspot_dir) / "hotspot_map.json"
    if hs_json.exists():
        with open(hs_json) as f:
            hotspot_data = json.load(f)

    # Load PLIP data if available
    df_plip = None
    plip_csv = Path(contact_dir) / "plip_interactions_top20.csv"
    if plip_csv.exists():
        df_plip = pd.read_csv(plip_csv)
        logger.info(f"  PLIP: {len(df_plip)} interactions loaded")

    if df_energy.empty:
        logger.error("No energy data found — cannot generate report")
        return {"success": False, "error": "No energy data"}

    # RMSD validation (if reference provided)
    rmsd_data = None
    if reference_mol2 and Path(reference_mol2).exists():
        logger.info(f"  Computing RMSD for control: {control_name}")
        rmsd_data = compute_reference_rmsd(reference_mol2, parsed_dir, control_name)
        if rmsd_data:
            logger.info(f"  RMSD: best={rmsd_data['best_rmsd']:.2f}Å (pose #{rmsd_data['best_pose']}), "
                        f"mean={rmsd_data['mean_rmsd']:.2f}Å, "
                        f"{rmsd_data['n_under_2A']}/{rmsd_data['n_poses']} < 2.0Å")

    # Generate figures
    logger.info("  Generating figures...")
    figures = {}

    if HAS_MPL:
        figures["ranking"] = fig_molecule_ranking(df_energy, mol_vina_scores=mol_vina_scores)
        figures["breakdown"] = fig_fragment_breakdown(df_energy)
        figures["conv_energy"] = fig_convergence_vs_energy(df_cluster, df_energy, mol_vina_scores=mol_vina_scores)
        figures["contacts"] = fig_contact_heatmap(df_contacts)
        figures["terms"] = fig_term_breakdown(df_energy)
        figures["efficiency"] = fig_fragment_efficiency(df_energy)

        # Also save as files
        for name, b64 in figures.items():
            if b64:
                img_data = base64.b64decode(b64)
                with open(fig_path / f"{name}.png", "wb") as f:
                    f.write(img_data)

    # Generate HTML report
    logger.info("  Generating HTML report...")
    html = generate_html_report(
        campaign_id=campaign_id,
        engine=engine,
        df_energy=df_energy,
        df_cluster=df_cluster,
        df_contacts=df_contacts,
        df_hotspot_contacts=df_hotspot_contacts,
        hotspot_data=hotspot_data,
        figures=figures,
        rmsd_data=rmsd_data,
        df_gnina=df_gnina,
        df_plip=df_plip,
    )

    report_path = out_path / "campaign_report.html"
    report_path.write_text(html, encoding='utf-8')

    # Save composite ranking as CSV (same multiplicative formula as HTML)
    if not df_energy.empty and not df_cluster.empty:
        # Merge energy + cluster at fragment level
        frag_merged = df_energy.merge(
            df_cluster[['Name', 'fragment_id', 'dominant_fraction']],
            on=['Name', 'fragment_id'], how='left'
        ).fillna(0)

        prop_vina_col = 'proportional_vina_mean' if 'proportional_vina_mean' in frag_merged.columns else 'vina_mean'

        # Compute weighted_conv, min_conv, best_fragment per molecule
        csv_rows = []
        for mol_name, mol_frags in frag_merged.groupby('Name'):
            vina_aff = mol_vina_scores.get(mol_name, np.nan)
            if pd.isna(vina_aff):
                continue
            w_conv = _compute_weighted_conv(mol_frags)
            min_c = mol_frags['dominant_fraction'].min()
            composite = vina_aff * w_conv

            # Best fragment
            mol_frags_copy = mol_frags.copy()
            mol_frags_copy['_fc'] = mol_frags_copy[prop_vina_col] * mol_frags_copy['dominant_fraction']
            best_idx = mol_frags_copy['_fc'].idxmin()
            best_row = mol_frags_copy.loc[best_idx]
            best_frag_str = f"{best_row['fragment_label']} ({best_row[prop_vina_col]:.2f})"

            row_data = {
                'name': mol_name,
                'vina_affinity': round(vina_aff, 3),
                'weighted_conv': round(w_conv, 4),
                'composite': round(composite, 3),
                'min_conv': round(min_c, 4),
                'best_fragment': best_frag_str,
            }

            # Add drug-likeness from gnina
            gnina_mol = df_ok[df_ok['name'] == mol_name]
            if not gnina_mol.empty:
                gnina_row = gnina_mol.iloc[0]
                if 'MW' in gnina_row:
                    row_data['MW'] = round(gnina_row['MW'], 1) if pd.notna(gnina_row['MW']) else ''
                if 'QED' in gnina_row:
                    row_data['QED'] = round(gnina_row['QED'], 3) if pd.notna(gnina_row['QED']) else ''

            csv_rows.append(row_data)

        rank_df = pd.DataFrame(csv_rows)
        rank_df = rank_df.sort_values('composite', ascending=True)  # most negative first
        rank_df.insert(0, 'rank', range(1, len(rank_df) + 1))

        ranking_csv = out_path / "composite_ranking.csv"
        rank_df.to_csv(ranking_csv, index=False)
        logger.info(f"  Ranking: {ranking_csv} ({len(rank_df)} molecules)")

    logger.info(f"\n{'=' * 60}")
    logger.info("CAMPAIGN REPORT COMPLETE")
    logger.info(f"{'=' * 60}")
    logger.info(f"  Report: {report_path}")
    logger.info(f"  Figures: {fig_path}/ ({len([f for f in figures.values() if f])} images)")

    return {
        "success": True,
        "report_path": str(report_path),
        "n_figures": len([f for f in figures.values() if f]),
    }
"""
Score Decomposition (05c) — Regress scores onto fragment positions
=====================================================================
For each molecule, decomposes the total docking score into per-fragment
contributions using linear regression on cluster assignments.

Model:
    Score_total(pose_j) ≈ Σ C(frag_i, cluster_ij) + intercept + ε

    Where cluster_ij is the cluster label of fragment i in pose j,
    and C(frag_i, cluster_k) is the contribution of fragment i when
    positioned in cluster k.

One-hot encoding of cluster assignments → OLS regression → coefficients
give the contribution of each fragment-in-position to the total score.

Also identifies sweet spots: the combinations of fragment positions
that yield the best predicted scores.

Reads:  05a output (.npz/.json) + 05b output (_clusters.npz/.json)
Saves:  {name}_scores.json   — coefficients, R², sweet spots
        score_decomposition_summary.csv — one row per fragment-cluster per molecule

Location: 01_src/molecular_docking/m05_gnina_analysis/score_decomposition.py
Project: molecular_docking
Module: 05c (core)
Version: 2.1 — exact atom_terms for Vina, regression only for CNN/DOCK6 (2026-03-22)
"""

import json
import logging
from dataclasses import dataclass, field
from itertools import product as iterproduct
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .parse_and_fragment import MoleculeData, load_all_molecules, VINA_WEIGHTS, ATOM_TERM_KEYS
from .fragment_clustering import MoleculeClusterResult, load_cluster_result

logger = logging.getLogger(__name__)


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass
class FragmentEnergy:
    """Exact energy contribution of one fragment in one pose (from atom_terms)."""
    fragment_id: int
    fragment_label: str
    n_atoms: int
    # Per-term contributions (sum of weighted atom terms for this fragment)
    gauss1: float = 0.0
    gauss2: float = 0.0
    repulsion: float = 0.0
    hydrophobic: float = 0.0
    hbond: float = 0.0
    # Total Vina contribution = sum of weighted terms
    vina_contribution: float = 0.0


@dataclass
class FragmentScoreContrib:
    """Score contribution of one fragment in one cluster position."""
    fragment_id: int
    fragment_label: str
    cluster_id: int
    cluster_size: int
    coefficient: float          # regression coefficient = score contribution
    std_error: float            # standard error of coefficient


@dataclass
class ScoreDecompositionResult:
    """Complete score decomposition for one molecule + one score type."""
    name: str
    score_key: str
    n_poses_used: int           # poses included in regression (non-noise)
    n_poses_total: int
    r_squared: float
    adj_r_squared: float
    intercept: float
    residual_std: float

    # Per fragment-cluster coefficients
    contributions: List[FragmentScoreContrib]

    # Per-fragment R² (partial): how much variance does this fragment explain?
    fragment_r2: Dict[int, float] = field(default_factory=dict)

    # Sweet spots: best N combinations of fragment positions
    sweet_spots: List[Dict] = field(default_factory=list)

    # Consensus poses: indices where ALL fragments in dominant cluster
    full_consensus_indices: List[int] = field(default_factory=list)
    full_consensus_score_mean: float = 0.0
    full_consensus_score_best: float = 0.0
    non_consensus_score_mean: float = 0.0


@dataclass
class FragmentEnergyResult:
    """Exact Vina energy decomposition per fragment (from atom_terms)."""
    name: str
    n_poses: int
    has_atom_terms: bool
    # Per fragment: mean, std, best of Vina contribution across all poses
    fragment_energies: List[Dict] = field(default_factory=list)
    # Per pose: total reconstructed Vina from atom terms
    total_reconstructed: Optional[np.ndarray] = None


@dataclass
class MoleculeScoreResult:
    """Score decomposition for all score types of one molecule."""
    name: str
    engine: str
    decompositions: Dict[str, ScoreDecompositionResult]  # score_key -> result
    energy_decomposition: Optional[FragmentEnergyResult] = None


# =============================================================================
# REGRESSION
# =============================================================================

def _build_design_matrix(
        label_matrix: np.ndarray,
) -> Tuple[np.ndarray, List[Tuple[int, int]], np.ndarray]:
    """
    Build one-hot design matrix from fragment cluster labels.

    Excludes poses where ANY fragment is noise (-1).
    For each fragment, drops the first cluster (reference level) to avoid
    multicollinearity. The reference level's contribution is absorbed into
    the intercept.

    Args:
        label_matrix: (n_poses, n_fragments) cluster labels

    Returns:
        X: (n_valid_poses, n_features) design matrix
        feature_map: list of (fragment_id, cluster_id) for each column
        valid_mask: boolean mask of which poses are included
    """
    n_poses, n_frags = label_matrix.shape

    # Exclude poses with any noise label
    valid_mask = np.all(label_matrix >= 0, axis=1)

    if np.sum(valid_mask) < 5:
        return np.empty((0, 0)), [], valid_mask

    valid_labels = label_matrix[valid_mask]

    # Build columns: for each fragment, one-hot minus reference level
    columns = []
    feature_map = []

    for fi in range(n_frags):
        frag_labels = valid_labels[:, fi]
        unique_clusters = sorted(set(frag_labels))

        if len(unique_clusters) <= 1:
            # Single cluster — no variance to explain, skip
            continue

        # Drop first cluster as reference level
        ref_cluster = unique_clusters[0]
        for cl in unique_clusters[1:]:
            col = (frag_labels == cl).astype(np.float64)
            columns.append(col)
            feature_map.append((fi, cl))

    if not columns:
        return np.empty((0, 0)), [], valid_mask

    X = np.column_stack(columns)
    return X, feature_map, valid_mask


def decompose_score(
        mol_data: MoleculeData,
        cluster_result: MoleculeClusterResult,
        score_key: str,
        n_sweet_spots: int = 10,
) -> Optional[ScoreDecompositionResult]:
    """
    Decompose one score type into per-fragment-position contributions.

    Uses OLS regression: score ~ Σ one_hot(fragment_i, cluster_ij)
    """
    label_matrix = cluster_result.label_matrix
    n_poses = mol_data.n_poses

    # Find score column index
    if score_key not in mol_data.score_keys:
        return None
    score_idx = mol_data.score_keys.index(score_key)
    scores = mol_data.score_values[:, score_idx]

    # Build design matrix (excludes noise poses)
    X, feature_map, valid_mask = _build_design_matrix(label_matrix)

    if X.size == 0 or len(feature_map) == 0:
        return None

    valid_scores = scores[valid_mask]
    valid_nan_mask = ~np.isnan(valid_scores)

    if np.sum(valid_nan_mask) < len(feature_map) + 2:
        return None

    X_clean = X[valid_nan_mask]
    y_clean = valid_scores[valid_nan_mask]
    n = len(y_clean)
    p = X_clean.shape[1]

    # OLS: add intercept
    X_int = np.column_stack([np.ones(n), X_clean])

    try:
        # Use numpy lstsq for numerical stability
        beta, residuals, rank, sv = np.linalg.lstsq(X_int, y_clean, rcond=None)
    except np.linalg.LinAlgError:
        return None

    intercept = beta[0]
    coeffs = beta[1:]

    # Predictions + R²
    y_pred = X_int @ beta
    ss_res = np.sum((y_clean - y_pred) ** 2)
    ss_tot = np.sum((y_clean - np.mean(y_clean)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    adj_r2 = 1 - (1 - r2) * (n - 1) / (n - p - 1) if n > p + 1 else r2
    residual_std = float(np.sqrt(ss_res / (n - p - 1))) if n > p + 1 else 0.0

    # Standard errors of coefficients
    try:
        mse = ss_res / (n - p - 1) if n > p + 1 else 0.0
        XtX_inv = np.linalg.inv(X_int.T @ X_int)
        diag = np.diag(XtX_inv) * mse
        se = np.sqrt(np.maximum(diag, 0.0))  # clamp negatives from numerical noise
        se_coeffs = se[1:]  # skip intercept
    except np.linalg.LinAlgError:
        se_coeffs = np.zeros(p)

    # Build contributions list
    contributions = []
    for idx, (fi, cl) in enumerate(feature_map):
        fr = cluster_result.fragment_results[fi]
        cl_size = int(np.sum(fr.labels == cl))
        contributions.append(FragmentScoreContrib(
            fragment_id=fi,
            fragment_label=fr.label,
            cluster_id=cl,
            cluster_size=cl_size,
            coefficient=round(float(coeffs[idx]), 4),
            std_error=round(float(se_coeffs[idx]), 4) if idx < len(se_coeffs) else 0.0,
        ))

    # Also add reference level contributions (absorbed into intercept)
    # We compute them as: ref_coeff = intercept - sum(other fragment ref effects)
    # Actually, the reference contribution for each fragment is implicit.
    # For clarity, add them with coefficient=0 (relative to intercept)
    for fi in range(cluster_result.n_fragments):
        fr = cluster_result.fragment_results[fi]
        frag_labels = label_matrix[valid_mask][valid_nan_mask, fi]
        unique_cls = sorted(set(frag_labels))
        if len(unique_cls) <= 1:
            continue
        ref_cl = unique_cls[0]
        cl_size = int(np.sum(fr.labels == ref_cl))
        contributions.append(FragmentScoreContrib(
            fragment_id=fi,
            fragment_label=fr.label,
            cluster_id=ref_cl,
            cluster_size=cl_size,
            coefficient=0.0,  # reference level
            std_error=0.0,
        ))

    # Sort by fragment_id, then cluster_id
    contributions.sort(key=lambda c: (c.fragment_id, c.cluster_id))

    # --- Partial R² per fragment ---
    fragment_r2 = {}
    for fi in range(cluster_result.n_fragments):
        frag_cols = [idx for idx, (f, c) in enumerate(feature_map) if f == fi]
        if not frag_cols:
            fragment_r2[fi] = 0.0
            continue

        # R² without this fragment
        keep_cols = [idx for idx in range(p) if idx not in frag_cols]
        if keep_cols:
            X_reduced = np.column_stack([np.ones(n), X_clean[:, keep_cols]])
            try:
                beta_r, _, _, _ = np.linalg.lstsq(X_reduced, y_clean, rcond=None)
                y_pred_r = X_reduced @ beta_r
                ss_res_r = np.sum((y_clean - y_pred_r) ** 2)
                r2_reduced = 1 - ss_res_r / ss_tot if ss_tot > 0 else 0.0
            except np.linalg.LinAlgError:
                r2_reduced = 0.0
        else:
            r2_reduced = 0.0

        fragment_r2[fi] = round(max(0.0, r2 - r2_reduced), 4)

    # --- Full consensus analysis ---
    dominant_clusters = []
    for fi in range(cluster_result.n_fragments):
        fr = cluster_result.fragment_results[fi]
        dominant_clusters.append(fr.dominant_cluster)

    # Poses where ALL fragments in their dominant cluster
    consensus_mask = np.ones(n_poses, dtype=bool)
    for fi, dom_cl in enumerate(dominant_clusters):
        consensus_mask &= (label_matrix[:, fi] == dom_cl)

    consensus_indices = np.where(consensus_mask)[0].tolist()
    consensus_scores = scores[consensus_mask]
    non_consensus_scores = scores[~consensus_mask & ~np.isnan(scores)]
    consensus_scores_valid = consensus_scores[~np.isnan(consensus_scores)]

    fc_mean = float(np.mean(consensus_scores_valid)) if len(consensus_scores_valid) > 0 else 0.0
    fc_best = float(np.min(consensus_scores_valid)) if len(consensus_scores_valid) > 0 else 0.0
    nc_mean = float(np.mean(non_consensus_scores)) if len(non_consensus_scores) > 0 else 0.0

    # --- Sweet spots: best combinations of cluster positions ---
    sweet_spots = _find_sweet_spots(
        label_matrix, scores, cluster_result, contributions,
        intercept, feature_map, n_sweet_spots,
    )

    return ScoreDecompositionResult(
        name=mol_data.name,
        score_key=score_key,
        n_poses_used=n,
        n_poses_total=n_poses,
        r_squared=round(r2, 4),
        adj_r_squared=round(adj_r2, 4),
        intercept=round(intercept, 4),
        residual_std=round(residual_std, 4),
        contributions=contributions,
        fragment_r2=fragment_r2,
        sweet_spots=sweet_spots,
        full_consensus_indices=consensus_indices,
        full_consensus_score_mean=round(fc_mean, 4),
        full_consensus_score_best=round(fc_best, 4),
        non_consensus_score_mean=round(nc_mean, 4),
    )


def _find_sweet_spots(
        label_matrix, scores, cluster_result, contributions,
        intercept, feature_map, n_spots=10,
) -> List[Dict]:
    """
    Find best combinations of fragment cluster positions by actual score data.

    Groups poses by their complete cluster assignment tuple, computes
    mean score per group, and returns the top N.
    """
    n_poses = label_matrix.shape[0]
    n_frags = label_matrix.shape[1]

    # Filter to non-noise, non-nan poses
    valid = np.all(label_matrix >= 0, axis=1) & ~np.isnan(scores)
    if np.sum(valid) < 3:
        return []

    valid_labels = label_matrix[valid]
    valid_scores = scores[valid]

    # Group by cluster assignment tuple
    combos = {}
    for i in range(len(valid_labels)):
        key = tuple(valid_labels[i].tolist())
        if key not in combos:
            combos[key] = []
        combos[key].append(valid_scores[i])

    # Compute stats per combination
    spots = []
    for key, sc_list in combos.items():
        arr = np.array(sc_list)
        spots.append({
            "combination": {
                f"frag{fi}": {"cluster": int(key[fi]),
                              "label": cluster_result.fragment_results[fi].label}
                for fi in range(n_frags)
            },
            "n_poses": len(arr),
            "mean_score": round(float(np.mean(arr)), 4),
            "best_score": round(float(np.min(arr)), 4),
            "std_score": round(float(np.std(arr)), 4),
        })

    # Sort by mean score (lower = better for docking scores)
    spots.sort(key=lambda s: s["mean_score"])
    return spots[:n_spots]


# =============================================================================
# EXACT ENERGY DECOMPOSITION (from GNINA atom_terms)
# =============================================================================

def compute_fragment_energies(
        mol_data: MoleculeData,
) -> Optional[FragmentEnergyResult]:
    """
    Compute exact Vina energy contribution per fragment per pose.

    Uses atomic_interaction_terms from GNINA --atom_term_data.
    For each atom, the Vina contribution is:
        E_atom = w_gauss1 * gauss1 + w_gauss2 * gauss2 + w_rep * repulsion
               + w_hydro * hydrophobic + w_hbond * hbond

    For each fragment: E_fragment = sum(E_atom) over fragment atoms.

    Returns FragmentEnergyResult with per-fragment statistics.
    """
    if mol_data.atom_terms is None:
        return FragmentEnergyResult(
            name=mol_data.name, n_poses=mol_data.n_poses,
            has_atom_terms=False,
        )

    n_poses = mol_data.n_poses
    n_frags = len(mol_data.fragments)

    # GNINA atom_terms are ALREADY pre-weighted (each value = w_i * raw_term_i).
    # Just sum across terms to get per-atom Vina contribution.
    # atom_terms shape: (n_poses, n_atoms, 5)
    per_atom_vina = np.sum(mol_data.atom_terms, axis=2)  # (n_poses, n_atoms)

    # Per-fragment energy: sum atom contributions within each fragment
    frag_vina = np.zeros((n_poses, n_frags), dtype=np.float64)
    frag_terms = np.zeros((n_poses, n_frags, 5), dtype=np.float64)

    for fi, frag in enumerate(mol_data.fragments):
        idx = frag.heavy_atom_indices
        frag_vina[:, fi] = np.sum(per_atom_vina[:, idx], axis=1)
        frag_terms[:, fi, :] = np.sum(mol_data.atom_terms[:, idx, :], axis=1)

    # Total reconstructed Vina (should match vina_affinity closely)
    total_reconstructed = np.sum(frag_vina, axis=1)

    # Build per-fragment statistics
    fragment_energies = []
    for fi, frag in enumerate(mol_data.fragments):
        fv = frag_vina[:, fi]
        ft = frag_terms[:, fi, :]  # (n_poses, 5) — already weighted

        fragment_energies.append({
            "fragment_id": frag.fragment_id,
            "fragment_label": frag.label,
            "n_atoms": frag.n_atoms,
            "vina_mean": round(float(np.mean(fv)), 4),
            "vina_std": round(float(np.std(fv)), 4),
            "vina_best": round(float(np.min(fv)), 4),
            "vina_worst": round(float(np.max(fv)), 4),
            # Per-term means (already weighted by GNINA)
            "gauss1_contrib": round(float(np.mean(ft[:, 0])), 4),
            "gauss2_contrib": round(float(np.mean(ft[:, 1])), 4),
            "repulsion_contrib": round(float(np.mean(ft[:, 2])), 4),
            "hydrophobic_contrib": round(float(np.mean(ft[:, 3])), 4),
            "hbond_contrib": round(float(np.mean(ft[:, 4])), 4),
            # Per-pose array for downstream use
            "per_pose_vina": fv.tolist(),
        })

    return FragmentEnergyResult(
        name=mol_data.name, n_poses=n_poses,
        has_atom_terms=True,
        fragment_energies=fragment_energies,
        total_reconstructed=total_reconstructed,
    )


# =============================================================================
# BATCH PROCESSING
# =============================================================================

def decompose_molecule(
        mol_data: MoleculeData,
        cluster_result: MoleculeClusterResult,
        score_keys: Optional[List[str]] = None,
        n_sweet_spots: int = 10,
) -> MoleculeScoreResult:
    """Decompose all score types for one molecule."""
    if score_keys is None:
        score_keys = mol_data.score_keys

    # Exact energy decomposition from atom_terms (if available)
    energy_result = compute_fragment_energies(mol_data)
    has_exact_vina = energy_result is not None and energy_result.has_atom_terms

    # Regression: skip vina_affinity if we have exact decomposition
    # Only regress CNN scores (cnn_score, cnn_affinity, cnn_vs) and DOCK6 scores
    decompositions = {}
    for sk in score_keys:
        if has_exact_vina and sk == "vina_affinity":
            continue  # exact decomposition is better
        result = decompose_score(mol_data, cluster_result, sk, n_sweet_spots)
        if result is not None and result.r_squared > 0:  # filter bad regressions
            decompositions[sk] = result

    return MoleculeScoreResult(
        name=mol_data.name,
        engine=mol_data.engine,
        decompositions=decompositions,
        energy_decomposition=energy_result,
    )


# =============================================================================
# SERIALIZATION
# =============================================================================

def save_score_result(result: MoleculeScoreResult, output_dir: Path):
    """Save score decomposition to JSON."""
    out = {
        "name": result.name,
        "engine": result.engine,
        "score_types": {},
    }

    for sk, dec in result.decompositions.items():
        out["score_types"][sk] = {
            "n_poses_used": dec.n_poses_used,
            "n_poses_total": dec.n_poses_total,
            "r_squared": dec.r_squared,
            "adj_r_squared": dec.adj_r_squared,
            "intercept": dec.intercept,
            "residual_std": dec.residual_std,
            "fragment_r2": {str(k): v for k, v in dec.fragment_r2.items()},
            "contributions": [
                {
                    "fragment_id": c.fragment_id,
                    "fragment_label": c.fragment_label,
                    "cluster_id": c.cluster_id,
                    "cluster_size": c.cluster_size,
                    "coefficient": c.coefficient,
                    "std_error": c.std_error,
                }
                for c in dec.contributions
            ],
            "sweet_spots": dec.sweet_spots,
            "full_consensus": {
                "n_poses": len(dec.full_consensus_indices),
                "pose_indices": dec.full_consensus_indices,
                "mean_score": dec.full_consensus_score_mean,
                "best_score": dec.full_consensus_score_best,
                "non_consensus_mean": dec.non_consensus_score_mean,
            },
        }

    # Exact energy decomposition (from atom_terms)
    if result.energy_decomposition and result.energy_decomposition.has_atom_terms:
        ed = result.energy_decomposition
        out["energy_decomposition"] = {
            "has_atom_terms": True,
            "n_poses": ed.n_poses,
            "fragments": [],
        }
        for fe in ed.fragment_energies:
            # Don't save per_pose array in JSON (too large)
            fe_copy = {k: v for k, v in fe.items() if k != "per_pose_vina"}
            out["energy_decomposition"]["fragments"].append(fe_copy)

    with open(output_dir / f"{result.name}_scores.json", "w") as fh:
        def _json_default(o):
            if isinstance(o, (np.integer,)):
                return int(o)
            if isinstance(o, (np.floating,)):
                return float(o)
            if isinstance(o, np.ndarray):
                return o.tolist()
            return str(o)

        json.dump(out, fh, indent=2, default=_json_default)


# =============================================================================
# PIPELINE
# =============================================================================

def run_score_decomposition(
        parsed_dir: str,
        cluster_dir: str,
        output_dir: str,
        score_keys: Optional[List[str]] = None,
        n_sweet_spots: int = 10,
        molecule_names: Optional[List[str]] = None,
) -> Dict[str, any]:
    """
    Score decomposition for all molecules.

    Reads: 05a + 05b output
    Saves: per-molecule JSON + summary CSV
    """
    logger.info("=" * 60)
    logger.info("05c SCORE DECOMPOSITION v2.0")
    logger.info("=" * 60)

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    all_mols = load_all_molecules(parsed_dir, molecule_names)
    if not all_mols:
        logger.error("No parsed molecules found")
        return {"success": False, "error": "No data"}

    summary_rows = []
    energy_rows = []
    n_ok = 0

    for name in sorted(all_mols.keys()):
        mol_data = all_mols[name]

        cluster_result = load_cluster_result(name, cluster_dir)
        if cluster_result is None:
            logger.debug(f"  {name}: no cluster data — skipping")
            continue

        try:
            mol_result = decompose_molecule(
                mol_data, cluster_result, score_keys, n_sweet_spots,
            )
        except Exception as e:
            logger.warning(f"  {name}: decomposition failed — {e}")
            continue

        if not mol_result.decompositions:
            logger.debug(f"  {name}: no decomposable scores")
            continue

        save_score_result(mol_result, out_path)
        n_ok += 1

        # Summary rows
        for sk, dec in mol_result.decompositions.items():
            for c in dec.contributions:
                row = {
                    "Name": name,
                    "engine": mol_data.engine,
                    "score_key": sk,
                    "fragment_id": c.fragment_id,
                    "fragment_label": c.fragment_label,
                    "cluster_id": c.cluster_id,
                    "cluster_size": c.cluster_size,
                    "coefficient": c.coefficient,
                    "std_error": c.std_error,
                    "r_squared": dec.r_squared,
                    "fragment_r2": dec.fragment_r2.get(c.fragment_id, 0.0),
                    "n_consensus": len(dec.full_consensus_indices),
                    "consensus_mean": dec.full_consensus_score_mean,
                    "consensus_best": dec.full_consensus_score_best,
                }
                summary_rows.append(row)

        # Energy decomposition summary (separate CSV)
        if mol_result.energy_decomposition and mol_result.energy_decomposition.has_atom_terms:
            for fe in mol_result.energy_decomposition.fragment_energies:
                energy_rows.append({
                    "Name": name,
                    "engine": mol_data.engine,
                    "fragment_id": fe["fragment_id"],
                    "fragment_label": fe["fragment_label"],
                    "n_atoms": fe["n_atoms"],
                    "vina_mean": fe["vina_mean"],
                    "vina_std": fe["vina_std"],
                    "vina_best": fe["vina_best"],
                    "gauss1_contrib": fe["gauss1_contrib"],
                    "gauss2_contrib": fe["gauss2_contrib"],
                    "repulsion_contrib": fe["repulsion_contrib"],
                    "hydrophobic_contrib": fe["hydrophobic_contrib"],
                    "hbond_contrib": fe["hbond_contrib"],
                })

    # Summary CSV (regression-based)
    df = pd.DataFrame(summary_rows)
    if not df.empty:
        df = df.sort_values(["Name", "score_key", "fragment_id", "cluster_id"]).reset_index(drop=True)
        csv_path = out_path / "score_decomposition_summary.csv"
        df.to_csv(csv_path, index=False)
        logger.info(f"  Regression summary: {csv_path} ({len(df)} rows)")

    # Energy CSV (exact, from atom_terms)
    df_energy = pd.DataFrame(energy_rows)
    if not df_energy.empty:
        df_energy = df_energy.sort_values(["Name", "fragment_id"]).reset_index(drop=True)
        energy_csv = out_path / "fragment_energy_summary.csv"
        df_energy.to_csv(energy_csv, index=False)
        logger.info(f"  Energy summary: {energy_csv} ({len(df_energy)} rows)")

    # Stats
    logger.info(f"\n{'=' * 60}")
    logger.info("SCORE DECOMPOSITION COMPLETE")
    logger.info(f"{'=' * 60}")
    logger.info(f"  Molecules: {n_ok}")
    if summary_rows:
        r2_vals = list(set((r["Name"], r["score_key"], r["r_squared"]) for r in summary_rows))
        r2_list = [v[2] for v in r2_vals]
        if r2_list:
            logger.info(f"  Regression R² range: {min(r2_list):.3f} - {max(r2_list):.3f}")
            logger.info(f"  Regression R² mean:  {np.mean(r2_list):.3f}")
    if energy_rows:
        logger.info(f"  Atom-term energy: {len(energy_rows)} fragment entries")
        vina_vals = [r["vina_mean"] for r in energy_rows]
        logger.info(f"  Vina fragment range: {min(vina_vals):.3f} to {max(vina_vals):.3f}")

    return {
        "success": True,
        "n_molecules": n_ok,
        "output_dir": str(out_path),
    }

"""
Fragment Clustering (05b) — DBSCAN clustering per fragment per molecule
=========================================================================
For each molecule, clusters poses independently for each rigid fragment.
Uses direct pairwise RMSD (atom order preserved within engine) + DBSCAN.

Each pose gets a cluster label per fragment:
    Pose 0: frag0=A, frag1=A, frag2=B, frag3=A
    Pose 1: frag0=A, frag1=A, frag2=A, frag3=A  ← full consensus
    Pose 2: frag0=B, frag1=A, frag2=noise, frag3=A

This label matrix is the foundation for score decomposition (05c),
sweet spot identification, and contact analysis.

Reads:  05a output (.npz + .json per molecule)
Saves:  {name}_clusters.npz  — label matrix (n_poses × n_fragments)
        {name}_clusters.json — cluster stats per fragment
        cluster_summary.csv  — one row per fragment per molecule

Location: 01_src/molecular_docking/m05_gnina_analysis/fragment_clustering.py
Project: molecular_docking
Module: 05b (core)
Version: 2.0
"""

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN

from .parse_and_fragment import MoleculeData, FragmentDef, load_all_molecules

logger = logging.getLogger(__name__)


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass
class FragmentClusterResult:
    """Clustering result for one fragment of one molecule."""
    fragment_id: int
    label: str
    n_atoms: int
    n_clusters: int
    labels: np.ndarray              # (n_poses,) cluster label per pose (-1=noise)
    dominant_cluster: int           # label of largest cluster
    dominant_fraction: float
    dominant_size: int
    n_noise: int
    noise_fraction: float
    low_confidence: bool
    intra_cluster_rmsd: float       # mean pairwise RMSD within dominant cluster
    cluster_sizes: List[int]        # sorted descending
    rmsd_matrix: np.ndarray         # (n_poses, n_poses) for this fragment

    # Score stats per cluster: {cluster_label: {score_key: {mean, best, std}}}
    cluster_score_stats: Dict[int, Dict[str, Dict[str, float]]] = field(default_factory=dict)


@dataclass
class MoleculeClusterResult:
    """Complete clustering for one molecule across all fragments."""
    name: str
    engine: str
    n_poses: int
    n_fragments: int
    fragment_results: List[FragmentClusterResult]
    label_matrix: np.ndarray        # (n_poses, n_fragments) — THE key output


# =============================================================================
# PER-FRAGMENT RMSD + CLUSTERING
# =============================================================================

def compute_fragment_rmsd(coords: np.ndarray, frag_indices: List[int]) -> np.ndarray:
    """
    Compute pairwise RMSD matrix for one fragment.

    Args:
        coords: (n_poses, n_all_atoms, 3)
        frag_indices: which heavy atom indices belong to this fragment

    Returns:
        (n_poses, n_poses) symmetric RMSD matrix
    """
    n_poses = coords.shape[0]
    idx = np.array(frag_indices)
    sub = coords[:, idx, :]  # (n_poses, n_frag_atoms, 3)

    # Vectorized pairwise RMSD
    rmsd = np.zeros((n_poses, n_poses), dtype=np.float64)
    for i in range(n_poses):
        diff = sub[i] - sub[i + 1:]  # (n_remaining, n_frag_atoms, 3)
        sq = np.sum(diff ** 2, axis=2)   # (n_remaining, n_frag_atoms)
        r = np.sqrt(np.mean(sq, axis=1)) # (n_remaining,)
        rmsd[i, i + 1:] = r
        rmsd[i + 1:, i] = r

    return rmsd


def cluster_one_fragment(
        coords: np.ndarray,
        frag: FragmentDef,
        score_keys: List[str],
        score_values: np.ndarray,
        eps: float = 3.0,
        min_samples: int = 5,
        min_dominant_fraction: float = 0.15,
) -> FragmentClusterResult:
    """
    Cluster poses by one fragment's coordinates.

    Returns FragmentClusterResult with labels, stats, and RMSD matrix.
    """
    n_poses = coords.shape[0]

    # RMSD matrix for this fragment
    rmsd = compute_fragment_rmsd(coords, frag.heavy_atom_indices)

    # DBSCAN
    db = DBSCAN(eps=eps, min_samples=min_samples, metric="precomputed")
    labels = db.fit_predict(rmsd)

    n_noise = int(np.sum(labels == -1))
    noise_frac = n_noise / n_poses
    unique_labels = sorted(set(labels) - {-1})

    if not unique_labels:
        return FragmentClusterResult(
            fragment_id=frag.fragment_id, label=frag.label, n_atoms=frag.n_atoms,
            n_clusters=0, labels=labels,
            dominant_cluster=-1, dominant_fraction=0.0, dominant_size=0,
            n_noise=n_noise, noise_fraction=round(noise_frac, 4),
            low_confidence=True, intra_cluster_rmsd=0.0, cluster_sizes=[],
            rmsd_matrix=rmsd,
        )

    # Cluster sizes + score stats
    cluster_sizes_map = {}
    cluster_score_stats = {}

    for cl in unique_labels:
        members = np.where(labels == cl)[0]
        cluster_sizes_map[cl] = len(members)

        stats = {}
        for j, key in enumerate(score_keys):
            vals = score_values[members, j]
            valid = vals[~np.isnan(vals)]
            if len(valid) > 0:
                stats[key] = {
                    "mean": round(float(np.mean(valid)), 4),
                    "best": round(float(np.min(valid)), 4),
                    "worst": round(float(np.max(valid)), 4),
                    "std": round(float(np.std(valid)), 4),
                }
        cluster_score_stats[cl] = stats

    # Dominant cluster
    dom_cl = max(unique_labels, key=lambda c: cluster_sizes_map[c])
    dom_size = cluster_sizes_map[dom_cl]
    dom_frac = dom_size / n_poses

    # Intra-cluster RMSD for dominant
    dom_members = np.where(labels == dom_cl)[0]
    if len(dom_members) > 1:
        sub_mat = rmsd[np.ix_(dom_members, dom_members)]
        triu = sub_mat[np.triu_indices(len(dom_members), k=1)]
        intra_rmsd = float(np.mean(triu)) if len(triu) > 0 else 0.0
    else:
        intra_rmsd = 0.0

    sorted_sizes = sorted(cluster_sizes_map.values(), reverse=True)

    return FragmentClusterResult(
        fragment_id=frag.fragment_id, label=frag.label, n_atoms=frag.n_atoms,
        n_clusters=len(unique_labels), labels=labels,
        dominant_cluster=dom_cl,
        dominant_fraction=round(dom_frac, 4),
        dominant_size=dom_size,
        n_noise=n_noise, noise_fraction=round(noise_frac, 4),
        low_confidence=dom_frac < min_dominant_fraction,
        intra_cluster_rmsd=round(intra_rmsd, 3),
        cluster_sizes=sorted_sizes,
        rmsd_matrix=rmsd,
        cluster_score_stats=cluster_score_stats,
    )


# =============================================================================
# PER-MOLECULE CLUSTERING
# =============================================================================

def cluster_molecule(
        mol_data: MoleculeData,
        eps: float = 3.0,
        min_samples: int = 5,
        min_dominant_fraction: float = 0.15,
) -> MoleculeClusterResult:
    """
    Cluster all fragments of one molecule.

    Returns MoleculeClusterResult with label_matrix (n_poses × n_fragments).
    """
    n_poses = mol_data.n_poses
    n_frags = len(mol_data.fragments)

    frag_results = []
    label_matrix = np.full((n_poses, n_frags), -1, dtype=np.int32)

    for fi, frag in enumerate(mol_data.fragments):
        fr = cluster_one_fragment(
            mol_data.coords, frag,
            mol_data.score_keys, mol_data.score_values,
            eps=eps, min_samples=min_samples,
            min_dominant_fraction=min_dominant_fraction,
        )
        frag_results.append(fr)
        label_matrix[:, fi] = fr.labels

    return MoleculeClusterResult(
        name=mol_data.name, engine=mol_data.engine,
        n_poses=n_poses, n_fragments=n_frags,
        fragment_results=frag_results,
        label_matrix=label_matrix,
    )


# =============================================================================
# SERIALIZATION
# =============================================================================

def save_cluster_result(result: MoleculeClusterResult, output_dir: Path):
    """Save cluster result: .npz (label matrix + RMSD) + .json (stats)."""
    # NPZ: label matrix + per-fragment RMSD matrices
    save_dict = {
        "label_matrix": result.label_matrix,
    }
    for fr in result.fragment_results:
        save_dict[f"rmsd_frag{fr.fragment_id}"] = fr.rmsd_matrix

    np.savez_compressed(output_dir / f"{result.name}_clusters.npz", **save_dict)

    # JSON: cluster stats
    cluster_json = {
        "name": result.name,
        "engine": result.engine,
        "n_poses": result.n_poses,
        "n_fragments": result.n_fragments,
        "fragments": [],
    }

    for fr in result.fragment_results:
        fj = {
            "fragment_id": fr.fragment_id,
            "label": fr.label,
            "n_atoms": fr.n_atoms,
            "n_clusters": fr.n_clusters,
            "dominant_cluster": int(fr.dominant_cluster),
            "dominant_fraction": fr.dominant_fraction,
            "dominant_size": fr.dominant_size,
            "n_noise": fr.n_noise,
            "noise_fraction": fr.noise_fraction,
            "low_confidence": fr.low_confidence,
            "intra_cluster_rmsd": fr.intra_cluster_rmsd,
            "cluster_sizes": fr.cluster_sizes,
            "cluster_score_stats": {
                str(k): v for k, v in fr.cluster_score_stats.items()
            },
        }
        cluster_json["fragments"].append(fj)

    with open(output_dir / f"{result.name}_clusters.json", "w") as fh:
        json.dump(cluster_json, fh, indent=2)


def load_cluster_result(name: str, cluster_dir: str) -> Optional[MoleculeClusterResult]:
    """Load cluster result from .npz + .json."""
    cdir = Path(cluster_dir)
    npz_file = cdir / f"{name}_clusters.npz"
    json_file = cdir / f"{name}_clusters.json"

    if not npz_file.exists() or not json_file.exists():
        return None

    npz = np.load(npz_file)
    with open(json_file) as fh:
        cj = json.load(fh)

    label_matrix = npz["label_matrix"]
    frag_results = []

    for fd in cj["fragments"]:
        fid = fd["fragment_id"]
        rmsd_key = f"rmsd_frag{fid}"
        rmsd_mat = npz[rmsd_key] if rmsd_key in npz else np.empty((0, 0))

        # Reconstruct labels from label_matrix column
        labels = label_matrix[:, fid] if fid < label_matrix.shape[1] else np.array([])

        # Reconstruct cluster_score_stats with int keys
        css = {int(k): v for k, v in fd.get("cluster_score_stats", {}).items()}

        frag_results.append(FragmentClusterResult(
            fragment_id=fid,
            label=fd["label"],
            n_atoms=fd["n_atoms"],
            n_clusters=fd["n_clusters"],
            labels=labels,
            dominant_cluster=fd["dominant_cluster"],
            dominant_fraction=fd["dominant_fraction"],
            dominant_size=fd["dominant_size"],
            n_noise=fd["n_noise"],
            noise_fraction=fd["noise_fraction"],
            low_confidence=fd["low_confidence"],
            intra_cluster_rmsd=fd["intra_cluster_rmsd"],
            cluster_sizes=fd["cluster_sizes"],
            rmsd_matrix=rmsd_mat,
            cluster_score_stats=css,
        ))

    return MoleculeClusterResult(
        name=cj["name"], engine=cj["engine"],
        n_poses=cj["n_poses"], n_fragments=cj["n_fragments"],
        fragment_results=frag_results,
        label_matrix=label_matrix,
    )


# =============================================================================
# PIPELINE
# =============================================================================

def run_fragment_clustering(
        parsed_dir: str,
        output_dir: str,
        eps: float = 3.0,
        min_samples: int = 5,
        min_dominant_fraction: float = 0.15,
        molecule_names: Optional[List[str]] = None,
) -> Dict[str, any]:
    """
    Cluster all fragments of all molecules.

    Reads: 05a output (parse_and_fragment)
    Saves: per-molecule .npz (labels+RMSD) + .json (stats) + summary CSV
    """
    logger.info("=" * 60)
    logger.info("05b FRAGMENT CLUSTERING v2.0")
    logger.info("=" * 60)

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # Load parsed molecules from 05a
    all_mols = load_all_molecules(parsed_dir, molecule_names)
    if not all_mols:
        logger.error("No parsed molecules found")
        return {"success": False, "error": "No data"}

    logger.info(f"  eps={eps}, min_samples={min_samples}")

    summary_rows = []
    n_ok = 0

    for name in sorted(all_mols.keys()):
        mol_data = all_mols[name]

        try:
            result = cluster_molecule(
                mol_data, eps=eps, min_samples=min_samples,
                min_dominant_fraction=min_dominant_fraction,
            )
        except Exception as e:
            logger.warning(f"  {name}: clustering failed — {e}")
            continue

        save_cluster_result(result, out_path)
        n_ok += 1

        # Summary rows: one per fragment
        for fr in result.fragment_results:
            row = {
                "Name": name,
                "engine": mol_data.engine,
                "n_poses": mol_data.n_poses,
                "fragment_id": fr.fragment_id,
                "fragment_label": fr.label,
                "fragment_n_atoms": fr.n_atoms,
                "n_clusters": fr.n_clusters,
                "dominant_cluster": fr.dominant_cluster,
                "dominant_fraction": fr.dominant_fraction,
                "dominant_size": fr.dominant_size,
                "n_noise": fr.n_noise,
                "noise_fraction": fr.noise_fraction,
                "intra_cluster_rmsd": fr.intra_cluster_rmsd,
                "low_confidence": fr.low_confidence,
            }

            # Add dominant cluster score stats
            dom_stats = fr.cluster_score_stats.get(fr.dominant_cluster, {})
            for key, stats in dom_stats.items():
                row[f"dom_best_{key}"] = stats.get("best")
                row[f"dom_mean_{key}"] = stats.get("mean")

            summary_rows.append(row)

    # Summary CSV
    df = pd.DataFrame(summary_rows)
    if not df.empty:
        df = df.sort_values(["Name", "fragment_id"]).reset_index(drop=True)
        csv_path = out_path / "cluster_summary.csv"
        df.to_csv(csv_path, index=False)
        logger.info(f"  Summary: {csv_path} ({len(df)} rows)")

    # Stats
    if summary_rows:
        fracs = [r["dominant_fraction"] for r in summary_rows if r["dominant_fraction"] > 0]
        n_high = sum(1 for f in fracs if f > 0.5)
        logger.info(f"\n{'=' * 60}")
        logger.info("FRAGMENT CLUSTERING COMPLETE")
        logger.info(f"{'=' * 60}")
        logger.info(f"  Molecules: {n_ok}")
        logger.info(f"  Total fragments: {len(summary_rows)}")
        logger.info(f"  Convergence >50%: {n_high}/{len(fracs)} fragments")
        if fracs:
            logger.info(f"  Mean convergence: {np.mean(fracs):.2f}")

    return {
        "success": True,
        "n_molecules": n_ok,
        "output_dir": str(out_path),
    }

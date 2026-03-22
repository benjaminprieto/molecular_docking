"""
Binding Site Hotspots (05d) — Cross-molecule fragment position analysis
=========================================================================
After per-molecule clustering (05b) and score decomposition (05c), this
module asks: do fragments from DIFFERENT molecules converge to the same
regions of the binding site?

Algorithm:
    1. Collect centroids of dominant clusters across all molecules
       (centroid = mean coordinate of fragment atoms in the medoid pose)
    2. Cluster these centroids using DBSCAN → "hotspots"
       Each hotspot = a region where multiple molecules place fragments
    3. For each hotspot: which molecules, which fragment types, what scores
    4. Score-weighted hotspot ranking: hotspots where high-scoring fragments
       cluster are pharmacophoric anchors

This is a pharmacophore-style analysis derived purely from docking data:
if 45 different molecules independently place an aromatic ring in the
same pocket region, that's a pharmacophoric hotspot.

Reads:  05a (.npz/.json) + 05b (_clusters.npz/.json) + 05c (_scores.json)
Saves:  hotspot_map.json, hotspot_summary.csv, hotspot_overlay.cxc

Location: 01_src/molecular_docking/m05_gnina_analysis/binding_site_hotspots.py
Project: molecular_docking
Module: 05d (core)
Version: 2.0
"""

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN

from .parse_and_fragment import MoleculeData, FragmentDef, load_all_molecules
from .fragment_clustering import MoleculeClusterResult, load_cluster_result

logger = logging.getLogger(__name__)


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass
class FragmentCentroid:
    """Centroid of one fragment's dominant cluster for one molecule."""
    molecule_name: str
    fragment_id: int
    fragment_label: str
    fragment_n_atoms: int
    is_ring: bool
    cluster_id: int
    cluster_size: int
    dominant_fraction: float
    centroid: np.ndarray            # (3,) xyz
    medoid_pose_index: int
    # Score info (from 05c if available)
    score_contribution: float = 0.0   # regression coefficient
    best_score: float = 0.0
    mean_score: float = 0.0


@dataclass
class Hotspot:
    """A region where fragments from multiple molecules converge."""
    hotspot_id: int
    center: np.ndarray              # (3,) mean of member centroids
    radius: float                   # max distance from center to any member
    n_molecules: int                # unique molecules contributing
    n_fragments: int                # total fragment centroids in hotspot
    members: List[FragmentCentroid]
    # Aggregated stats
    molecules: List[str]
    ring_fraction: float            # fraction of members that are ring fragments
    mean_dominant_fraction: float   # mean convergence of contributing clusters
    mean_score_contribution: float  # mean score contribution of members
    best_score_contribution: float  # best (most negative) score contribution


# =============================================================================
# CENTROID EXTRACTION
# =============================================================================

def _find_medoid_index(rmsd_matrix: np.ndarray, member_indices: np.ndarray) -> int:
    """Find medoid: pose with lowest mean RMSD to cluster members."""
    if len(member_indices) <= 1:
        return int(member_indices[0])
    sub = rmsd_matrix[np.ix_(member_indices, member_indices)]
    return int(member_indices[int(np.argmin(np.mean(sub, axis=1)))])


def extract_centroids(
        mol_data: MoleculeData,
        cluster_result: MoleculeClusterResult,
        score_data: Optional[Dict] = None,
        include_non_dominant: bool = False,
        min_dominant_fraction: float = 0.15,
) -> List[FragmentCentroid]:
    """
    Extract fragment centroids from one molecule's clustering results.

    By default, only extracts the dominant cluster per fragment.
    """
    centroids = []

    for fr in cluster_result.fragment_results:
        frag = mol_data.fragments[fr.fragment_id]

        if fr.n_clusters == 0:
            continue

        # Which clusters to process
        if include_non_dominant:
            clusters_to_process = sorted(set(fr.labels) - {-1})
        else:
            clusters_to_process = [fr.dominant_cluster]

        for cl in clusters_to_process:
            members = np.where(fr.labels == cl)[0]
            if len(members) == 0:
                continue

            cl_size = len(members)
            cl_frac = cl_size / mol_data.n_poses

            if cl_frac < min_dominant_fraction and cl != fr.dominant_cluster:
                continue

            # Find medoid
            medoid_idx = _find_medoid_index(fr.rmsd_matrix, members)

            # Compute centroid of fragment atoms in medoid pose
            frag_coords = mol_data.coords[medoid_idx, frag.heavy_atom_indices, :]
            centroid = np.mean(frag_coords, axis=0)

            # Score contribution from 05c (if available)
            score_contrib = 0.0
            best_sc = 0.0
            mean_sc = 0.0
            if score_data:
                for sk, sd in score_data.get("score_types", {}).items():
                    for c in sd.get("contributions", []):
                        if c["fragment_id"] == fr.fragment_id and c["cluster_id"] == cl:
                            score_contrib = c["coefficient"]
                    # Cluster score stats from 05b
                    cs = fr.cluster_score_stats.get(cl, {})
                    for sk2, stats in cs.items():
                        if "vina" in sk2 or "Grid" in sk2:
                            best_sc = stats.get("best", 0.0)
                            mean_sc = stats.get("mean", 0.0)
                            break
                    break  # first score type only

            centroids.append(FragmentCentroid(
                molecule_name=mol_data.name,
                fragment_id=fr.fragment_id,
                fragment_label=fr.label,
                fragment_n_atoms=frag.n_atoms,
                is_ring=frag.is_ring,
                cluster_id=cl,
                cluster_size=cl_size,
                dominant_fraction=round(cl_frac, 4),
                centroid=centroid,
                medoid_pose_index=medoid_idx,
                score_contribution=score_contrib,
                best_score=best_sc,
                mean_score=mean_sc,
            ))

    return centroids


# =============================================================================
# HOTSPOT DETECTION
# =============================================================================

def detect_hotspots(
        all_centroids: List[FragmentCentroid],
        eps: float = 2.0,
        min_molecules: int = 3,
) -> List[Hotspot]:
    """
    Cluster fragment centroids across molecules to find binding site hotspots.

    Args:
        all_centroids: centroids from all molecules
        eps: DBSCAN eps in Angstrom (spatial proximity for same hotspot)
        min_molecules: minimum unique molecules to form a hotspot

    Returns:
        List of Hotspot objects, sorted by n_molecules descending
    """
    if len(all_centroids) < min_molecules:
        return []

    # Build coordinate matrix
    coords = np.array([c.centroid for c in all_centroids])

    # Euclidean distance DBSCAN
    db = DBSCAN(eps=eps, min_samples=min_molecules, metric="euclidean")
    labels = db.fit_predict(coords)

    unique_labels = sorted(set(labels) - {-1})
    if not unique_labels:
        return []

    hotspots = []
    for hi, hl in enumerate(unique_labels):
        member_mask = labels == hl
        member_indices = np.where(member_mask)[0]
        members = [all_centroids[i] for i in member_indices]

        # Center and radius
        member_coords = coords[member_indices]
        center = np.mean(member_coords, axis=0)
        distances = np.linalg.norm(member_coords - center, axis=1)
        radius = float(np.max(distances))

        # Unique molecules
        mol_names = sorted(set(m.molecule_name for m in members))

        # Stats
        ring_frac = sum(1 for m in members if m.is_ring) / len(members)
        mean_dom_frac = float(np.mean([m.dominant_fraction for m in members]))

        score_contribs = [m.score_contribution for m in members if m.score_contribution != 0]
        mean_sc = float(np.mean(score_contribs)) if score_contribs else 0.0
        best_sc = float(np.min(score_contribs)) if score_contribs else 0.0

        hotspots.append(Hotspot(
            hotspot_id=hi,
            center=center,
            radius=round(radius, 2),
            n_molecules=len(mol_names),
            n_fragments=len(members),
            members=members,
            molecules=mol_names,
            ring_fraction=round(ring_frac, 3),
            mean_dominant_fraction=round(mean_dom_frac, 3),
            mean_score_contribution=round(mean_sc, 4),
            best_score_contribution=round(best_sc, 4),
        ))

    hotspots.sort(key=lambda h: h.n_molecules, reverse=True)
    for i, h in enumerate(hotspots):
        h.hotspot_id = i

    return hotspots


# =============================================================================
# CHIMERAX SCRIPT
# =============================================================================

HOTSPOT_COLORS = [
    "cornflower blue", "orange red", "forest green", "gold",
    "dark violet", "deep pink", "dark cyan", "sienna",
    "lime green", "salmon", "steel blue", "olive drab",
]


def generate_hotspot_cxc(
        hotspots: List[Hotspot],
        output_path: Path,
        receptor_path: Optional[str] = None,
        all_centroids: Optional[List[FragmentCentroid]] = None,
):
    """Generate ChimeraX script to visualize hotspots as spheres."""
    lines = []
    lines.append("# Binding Site Hotspot Visualization")
    lines.append("# Auto-generated by 05d binding_site_hotspots")
    lines.append("")

    # Receptor
    if receptor_path and Path(receptor_path).exists():
        lines.append(f"open {receptor_path}")
        lines.append("color #1 light gray")
        lines.append("surface #1")
        lines.append("transparency #1 70")
        lines.append("")

    # Create hotspot markers as BILD objects
    bild_path = output_path.parent / "hotspots.bild"
    bild_lines = []

    for h in hotspots:
        color = HOTSPOT_COLORS[h.hotspot_id % len(HOTSPOT_COLORS)]
        cx, cy, cz = h.center

        # Hotspot center sphere
        # BILD format: .color R G B then .sphere x y z radius
        bild_lines.append(f".comment Hotspot {h.hotspot_id}: {h.n_molecules} molecules, radius {h.radius:.1f} A")

        # Use transparency for hotspot region
        r_sphere = max(h.radius, 1.5)  # minimum visible radius
        bild_lines.append(f".transparency 0.6")

        # Color by hotspot index
        from_hex = {
            "cornflower blue": (0.39, 0.58, 0.93),
            "orange red": (1.0, 0.27, 0.0),
            "forest green": (0.13, 0.55, 0.13),
            "gold": (1.0, 0.84, 0.0),
            "dark violet": (0.58, 0.0, 0.83),
            "deep pink": (1.0, 0.08, 0.58),
            "dark cyan": (0.0, 0.55, 0.55),
            "sienna": (0.63, 0.32, 0.18),
            "lime green": (0.2, 0.8, 0.2),
            "salmon": (0.98, 0.5, 0.45),
            "steel blue": (0.27, 0.51, 0.71),
            "olive drab": (0.42, 0.56, 0.14),
        }
        rgb = from_hex.get(color, (0.5, 0.5, 0.5))
        bild_lines.append(f".color {rgb[0]:.2f} {rgb[1]:.2f} {rgb[2]:.2f}")
        bild_lines.append(f".sphere {cx:.3f} {cy:.3f} {cz:.3f} {r_sphere:.2f}")

        # Individual centroid dots (smaller, solid)
        bild_lines.append(f".transparency 0.0")
        for m in h.members:
            mx, my, mz = m.centroid
            bild_lines.append(f".sphere {mx:.3f} {my:.3f} {mz:.3f} 0.5")

        bild_lines.append("")

    # Noise centroids (gray, small)
    if all_centroids:
        hotspot_members = set()
        for h in hotspots:
            for m in h.members:
                hotspot_members.add((m.molecule_name, m.fragment_id, m.cluster_id))

        bild_lines.append(f".comment Unassigned centroids")
        bild_lines.append(f".color 0.60 0.60 0.60")
        bild_lines.append(f".transparency 0.3")
        for c in all_centroids:
            key = (c.molecule_name, c.fragment_id, c.cluster_id)
            if key not in hotspot_members:
                bild_lines.append(f".sphere {c.centroid[0]:.3f} {c.centroid[1]:.3f} {c.centroid[2]:.3f} 0.3")

    bild_path.write_text("\n".join(bild_lines))

    lines.append(f"open {bild_path}")
    lines.append("")

    # Labels for hotspots
    for h in hotspots:
        cx, cy, cz = h.center
        label_text = f"H{h.hotspot_id}({h.n_molecules}mol)"
        lines.append(f"# Hotspot {h.hotspot_id}: {h.n_molecules} molecules, "
                     f"ring_frac={h.ring_fraction:.0%}, radius={h.radius:.1f}A")

    lines.append("")
    lines.append("view")

    output_path.write_text("\n".join(lines))


# =============================================================================
# SERIALIZATION
# =============================================================================

def save_hotspot_results(
        hotspots: List[Hotspot],
        all_centroids: List[FragmentCentroid],
        output_dir: Path,
):
    """Save hotspot analysis: JSON + CSV."""

    # JSON: full details
    hotspot_json = {
        "n_hotspots": len(hotspots),
        "n_centroids_total": len(all_centroids),
        "n_centroids_assigned": sum(h.n_fragments for h in hotspots),
        "hotspots": [],
    }

    for h in hotspots:
        hj = {
            "hotspot_id": h.hotspot_id,
            "center": h.center.tolist(),
            "radius": h.radius,
            "n_molecules": h.n_molecules,
            "n_fragments": h.n_fragments,
            "molecules": h.molecules,
            "ring_fraction": h.ring_fraction,
            "mean_dominant_fraction": h.mean_dominant_fraction,
            "mean_score_contribution": h.mean_score_contribution,
            "best_score_contribution": h.best_score_contribution,
            "members": [
                {
                    "molecule": m.molecule_name,
                    "fragment_id": m.fragment_id,
                    "fragment_label": m.fragment_label,
                    "is_ring": m.is_ring,
                    "cluster_id": m.cluster_id,
                    "dominant_fraction": m.dominant_fraction,
                    "score_contribution": m.score_contribution,
                    "centroid": m.centroid.tolist(),
                }
                for m in h.members
            ],
        }
        hotspot_json["hotspots"].append(hj)

    with open(output_dir / "hotspot_map.json", "w") as fh:
        json.dump(hotspot_json, fh, indent=2)

    # CSV summary: one row per hotspot-member
    rows = []
    for h in hotspots:
        for m in h.members:
            rows.append({
                "hotspot_id": h.hotspot_id,
                "hotspot_n_molecules": h.n_molecules,
                "hotspot_center_x": round(h.center[0], 2),
                "hotspot_center_y": round(h.center[1], 2),
                "hotspot_center_z": round(h.center[2], 2),
                "hotspot_radius": h.radius,
                "molecule": m.molecule_name,
                "fragment_id": m.fragment_id,
                "fragment_label": m.fragment_label,
                "fragment_n_atoms": m.fragment_n_atoms,
                "is_ring": m.is_ring,
                "cluster_id": m.cluster_id,
                "dominant_fraction": m.dominant_fraction,
                "score_contribution": m.score_contribution,
                "best_score": m.best_score,
                "mean_score": m.mean_score,
                "centroid_x": round(m.centroid[0], 2),
                "centroid_y": round(m.centroid[1], 2),
                "centroid_z": round(m.centroid[2], 2),
            })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["hotspot_id", "molecule"]).reset_index(drop=True)
    df.to_csv(output_dir / "hotspot_summary.csv", index=False)

    # Also save a per-hotspot summary
    hs_rows = []
    for h in hotspots:
        frag_types = {}
        for m in h.members:
            ft = "ring" if m.is_ring else "linker"
            frag_types[ft] = frag_types.get(ft, 0) + 1

        hs_rows.append({
            "hotspot_id": h.hotspot_id,
            "center_x": round(h.center[0], 2),
            "center_y": round(h.center[1], 2),
            "center_z": round(h.center[2], 2),
            "radius": h.radius,
            "n_molecules": h.n_molecules,
            "n_fragments": h.n_fragments,
            "n_ring": frag_types.get("ring", 0),
            "n_linker": frag_types.get("linker", 0),
            "ring_fraction": h.ring_fraction,
            "mean_convergence": h.mean_dominant_fraction,
            "mean_score_contrib": h.mean_score_contribution,
            "best_score_contrib": h.best_score_contribution,
        })

    pd.DataFrame(hs_rows).to_csv(output_dir / "hotspot_overview.csv", index=False)


# =============================================================================
# PIPELINE
# =============================================================================

def run_binding_site_hotspots(
        parsed_dir: str,
        cluster_dir: str,
        output_dir: str,
        score_dir: Optional[str] = None,
        receptor_path: Optional[str] = None,
        hotspot_eps: float = 2.0,
        min_molecules: int = 3,
        min_dominant_fraction: float = 0.15,
        include_non_dominant: bool = False,
        molecule_names: Optional[List[str]] = None,
) -> Dict[str, any]:
    """
    Cross-molecule hotspot analysis.

    Reads: 05a + 05b + 05c (optional)
    Saves: hotspot_map.json, hotspot_summary.csv, hotspot_overlay.cxc
    """
    logger.info("=" * 60)
    logger.info("05d BINDING SITE HOTSPOTS v2.0")
    logger.info("=" * 60)

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # Load molecules
    all_mols = load_all_molecules(parsed_dir, molecule_names)
    if not all_mols:
        logger.error("No parsed molecules found")
        return {"success": False, "error": "No data"}

    # Collect centroids across all molecules
    all_centroids = []
    n_mols_with_centroids = 0

    for name in sorted(all_mols.keys()):
        mol_data = all_mols[name]
        cluster_result = load_cluster_result(name, cluster_dir)
        if cluster_result is None:
            continue

        # Load score data if available
        score_data = None
        if score_dir:
            sp = Path(score_dir) / f"{name}_scores.json"
            if sp.exists():
                try:
                    with open(sp) as fh:
                        score_data = json.load(fh)
                except Exception:
                    pass

        centroids = extract_centroids(
            mol_data, cluster_result, score_data,
            include_non_dominant=include_non_dominant,
            min_dominant_fraction=min_dominant_fraction,
        )

        if centroids:
            all_centroids.extend(centroids)
            n_mols_with_centroids += 1

    logger.info(f"  Collected {len(all_centroids)} centroids from {n_mols_with_centroids} molecules")

    if len(all_centroids) < min_molecules:
        logger.error("Not enough centroids for hotspot detection")
        return {"success": False, "error": "Insufficient data"}

    # Detect hotspots
    hotspots = detect_hotspots(all_centroids, eps=hotspot_eps, min_molecules=min_molecules)
    logger.info(f"  Detected {len(hotspots)} hotspots (eps={hotspot_eps}A, min_mol={min_molecules})")

    # Save results
    save_hotspot_results(hotspots, all_centroids, out_path)

    # ChimeraX script
    cxc_path = out_path / "hotspot_overlay.cxc"
    generate_hotspot_cxc(hotspots, cxc_path, receptor_path, all_centroids)

    # Log summary
    n_assigned = sum(h.n_fragments for h in hotspots)
    logger.info(f"\n{'=' * 60}")
    logger.info("BINDING SITE HOTSPOTS COMPLETE")
    logger.info(f"{'=' * 60}")
    logger.info(f"  Molecules: {n_mols_with_centroids}")
    logger.info(f"  Centroids: {len(all_centroids)} total, {n_assigned} in hotspots")
    logger.info(f"  Hotspots: {len(hotspots)}")

    for h in hotspots[:5]:
        ring_pct = f"{h.ring_fraction:.0%} ring"
        logger.info(f"    H{h.hotspot_id}: {h.n_molecules} molecules, "
                     f"r={h.radius:.1f}A, {ring_pct}, "
                     f"score_contrib={h.mean_score_contribution:.2f}")

    return {
        "success": True,
        "n_hotspots": len(hotspots),
        "n_centroids": len(all_centroids),
        "n_molecules": n_mols_with_centroids,
        "output_dir": str(out_path),
    }

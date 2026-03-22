"""
Structure Export (05e) — Export 3D structures for ChimeraX visualization
=========================================================================
Exports mol2 files of cluster representatives + ChimeraX scripts that
are informed by hotspot analysis (05d). Each exported structure knows
which hotspot it belongs to, and the ChimeraX script colors accordingly.

Per molecule:
  - Fragment representatives: medoid pose for each cluster, fragment atoms only
  - Full-molecule representatives: consensus poses, sweet spot poses
  - Per-hotspot overlay: all molecules contributing to a hotspot

Global:
  - Master ChimeraX script: receptor + all hotspots + top molecules

The medoid is the real pose closest to the cluster centroid (not averaged).

Reads:  05a + 05b + 05c (optional) + 05d (hotspots)
Saves:  mol2 files + .cxc scripts

Location: 01_src/molecular_docking/m05_gnina_analysis/structure_export.py
Project: molecular_docking
Module: 05e (core)
Version: 2.0
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from .parse_and_fragment import MoleculeData, load_all_molecules
from .fragment_clustering import MoleculeClusterResult, load_cluster_result

logger = logging.getLogger(__name__)


# =============================================================================
# FRAGMENT COLORS (consistent with 05d hotspot colors)
# =============================================================================

HOTSPOT_COLORS = [
    "cornflower blue", "orange red", "forest green", "gold",
    "dark violet", "deep pink", "dark cyan", "sienna",
    "lime green", "salmon", "steel blue", "olive drab",
]


# =============================================================================
# MOL2 WRITER
# =============================================================================

def _write_mol2(
        path: Path,
        mol_name: str,
        coords: np.ndarray,          # (n_all_atoms, 3) for ONE pose
        elements: List[str],
        atom_indices: List[int],      # which atoms to include
        label: str = "",
):
    """
    Write mol2 with subset of atoms from one pose.
    Minimal TRIPOS format compatible with ChimeraX/PyMOL.
    """
    n_atoms = len(atom_indices)
    name = f"{mol_name}_{label}" if label else mol_name

    lines = [
        "@<TRIPOS>MOLECULE",
        name,
        f" {n_atoms} 0 0 0 0",
        "SMALL",
        "NO_CHARGES",
        "",
        "@<TRIPOS>ATOM",
    ]

    for local_idx, global_idx in enumerate(atom_indices):
        elem = elements[global_idx] if global_idx < len(elements) else "X"
        x, y, z = coords[global_idx]
        sybyl = f"{elem}.3" if elem in ("C", "N", "O", "S", "P") else elem
        atom_name = f"{elem}{local_idx + 1}"
        lines.append(
            f"  {local_idx + 1:>5} {atom_name:<6} "
            f"{x:>10.4f} {y:>10.4f} {z:>10.4f} "
            f"{sybyl:<6} 1 LIG 0.0000"
        )

    lines.append("")
    path.write_text("\n".join(lines))


def _write_full_mol2(path, mol_name, coords, elements, label=""):
    """Write mol2 with ALL heavy atoms from one pose."""
    _write_mol2(path, mol_name, coords, elements,
                list(range(len(elements))), label)


# =============================================================================
# MEDOID SELECTION
# =============================================================================

def _find_medoid(rmsd_matrix: np.ndarray, members: np.ndarray) -> int:
    """Pose with lowest mean RMSD to cluster members."""
    if len(members) <= 1:
        return int(members[0])
    sub = rmsd_matrix[np.ix_(members, members)]
    return int(members[int(np.argmin(np.mean(sub, axis=1)))])


# =============================================================================
# PER-MOLECULE EXPORT
# =============================================================================

def export_molecule(
        mol_data: MoleculeData,
        cluster_result: MoleculeClusterResult,
        mol_dir: Path,
        hotspot_map: Optional[Dict] = None,
        score_data: Optional[Dict] = None,
        max_clusters_per_fragment: int = 3,
) -> Dict[str, str]:
    """
    Export 3D structures for one molecule.

    Returns {label: filepath} of all exported files.
    """
    mol_dir.mkdir(parents=True, exist_ok=True)
    exported = {}
    elements = mol_data.heavy_atom_elements

    # Build hotspot lookup: (molecule, fragment_id, cluster_id) -> hotspot_id
    hs_lookup = {}
    if hotspot_map:
        for hs in hotspot_map.get("hotspots", []):
            for m in hs.get("members", []):
                key = (m["molecule"], m["fragment_id"], m["cluster_id"])
                hs_lookup[key] = hs["hotspot_id"]

    # --- Fragment cluster representatives ---
    for fr in cluster_result.fragment_results:
        frag = mol_data.fragments[fr.fragment_id]
        frag_indices = frag.heavy_atom_indices

        unique_cls = sorted(set(fr.labels) - {-1},
                           key=lambda c: -int(np.sum(fr.labels == c)))

        for rank, cl in enumerate(unique_cls[:max_clusters_per_fragment]):
            members = np.where(fr.labels == cl)[0]
            if len(members) == 0:
                continue

            medoid_idx = _find_medoid(fr.rmsd_matrix, members)
            medoid_coords = mol_data.coords[medoid_idx]

            # Hotspot tag
            hs_id = hs_lookup.get((mol_data.name, fr.fragment_id, cl))
            hs_tag = f"_H{hs_id}" if hs_id is not None else ""

            # Fragment-only mol2
            label = f"frag{fr.fragment_id}_cl{cl}{hs_tag}"
            fname = f"{mol_data.name}_{label}.mol2"
            _write_mol2(mol_dir / fname, mol_data.name, medoid_coords,
                       elements, frag_indices, label)
            exported[label] = str(mol_dir / fname)

            # Full molecule at this medoid
            full_label = f"frag{fr.fragment_id}_cl{cl}{hs_tag}_full"
            full_fname = f"{mol_data.name}_{full_label}.mol2"
            _write_full_mol2(mol_dir / full_fname, mol_data.name,
                           medoid_coords, elements, full_label)
            exported[full_label] = str(mol_dir / full_fname)

    # --- Full consensus pose ---
    dom_cls = [fr.dominant_cluster for fr in cluster_result.fragment_results]
    consensus_mask = np.ones(mol_data.n_poses, dtype=bool)
    for fi, dc in enumerate(dom_cls):
        consensus_mask &= (cluster_result.label_matrix[:, fi] == dc)

    consensus_idx = np.where(consensus_mask)[0]
    if len(consensus_idx) > 0:
        if len(consensus_idx) == 1:
            best_idx = consensus_idx[0]
        else:
            fr0 = cluster_result.fragment_results[0]
            sub = fr0.rmsd_matrix[np.ix_(consensus_idx, consensus_idx)]
            best_idx = consensus_idx[int(np.argmin(np.mean(sub, axis=1)))]

        fpath = mol_dir / f"{mol_data.name}_full_consensus.mol2"
        _write_full_mol2(fpath, mol_data.name, mol_data.coords[best_idx],
                        elements, "full_consensus")
        exported["full_consensus"] = str(fpath)

    # --- Sweet spot poses (from 05c) ---
    if score_data:
        for sk, sd in score_data.get("score_types", {}).items():
            for si, spot in enumerate(sd.get("sweet_spots", [])[:3]):
                combo = spot.get("combination", {})
                match_mask = np.ones(mol_data.n_poses, dtype=bool)
                for fkey, fval in combo.items():
                    fi = int(fkey.replace("frag", ""))
                    cl = fval["cluster"]
                    if fi < cluster_result.label_matrix.shape[1]:
                        match_mask &= (cluster_result.label_matrix[:, fi] == cl)

                match_idx = np.where(match_mask)[0]
                if len(match_idx) > 0:
                    sk_idx = mol_data.score_keys.index(sk) if sk in mol_data.score_keys else 0
                    match_scores = mol_data.score_values[match_idx, sk_idx]
                    valid = ~np.isnan(match_scores)
                    if np.any(valid):
                        best_local = int(np.argmin(match_scores[valid]))
                        best_pose = match_idx[np.where(valid)[0][best_local]]

                        label = f"sweet{si + 1}_{sk}"
                        fpath = mol_dir / f"{mol_data.name}_{label}.mol2"
                        _write_full_mol2(fpath, mol_data.name,
                                       mol_data.coords[best_pose], elements, label)
                        exported[label] = str(fpath)
            break  # first score type only

    return exported


# =============================================================================
# CHIMERAX SCRIPT — PER MOLECULE
# =============================================================================

def _generate_mol_cxc(
        mol_name: str,
        mol_dir: Path,
        exported: Dict[str, str],
        receptor_path: Optional[str],
        hotspot_map: Optional[Dict],
) -> str:
    """ChimeraX script for one molecule."""
    lines = [
        f"# ChimeraX: {mol_name}",
        f"# Auto-generated by 05e structure_export",
        "",
    ]

    if receptor_path and Path(receptor_path).exists():
        lines.append(f"open {receptor_path}")
        lines.append("color #1 light gray")
        lines.append("surface #1")
        lines.append("transparency #1 70")
        lines.append("")

    model = 2
    for key in sorted(exported.keys()):
        if key.endswith("_full") or "sweet" in key or "consensus" in key:
            continue  # load fragments first
        fp = exported[key]
        lines.append(f"open {fp}")

        # Color by hotspot if tagged
        if "_H" in key:
            try:
                hs_id = int(key.split("_H")[1].split("_")[0])
                color = HOTSPOT_COLORS[hs_id % len(HOTSPOT_COLORS)]
            except (ValueError, IndexError):
                color = "gray"
        else:
            try:
                fi = int(key.split("frag")[1].split("_")[0])
                color = HOTSPOT_COLORS[fi % len(HOTSPOT_COLORS)]
            except (ValueError, IndexError):
                color = "gray"

        lines.append(f"color #{model} {color}")
        lines.append(f"style #{model} ball")
        model += 1

    # Load consensus/sweet spots
    for key in sorted(exported.keys()):
        if "consensus" in key or "sweet" in key:
            lines.append(f"open {exported[key]}")
            lines.append(f"color #{model} tan")
            lines.append(f"style #{model} stick")
            model += 1

    lines.extend(["", "view", ""])
    return "\n".join(lines)


# =============================================================================
# CHIMERAX SCRIPT — GLOBAL (HOTSPOT-CENTRIC)
# =============================================================================

def _generate_global_cxc(
        output_dir: Path,
        all_exported: Dict[str, Dict[str, str]],
        receptor_path: Optional[str],
        hotspot_map: Optional[Dict],
        bild_path: Optional[Path],
) -> str:
    """Master ChimeraX script: receptor + hotspot spheres + top representatives."""
    lines = [
        "# Master visualization: all hotspots + representatives",
        "# Auto-generated by 05e structure_export",
        "",
    ]

    if receptor_path and Path(receptor_path).exists():
        lines.append(f"open {receptor_path}")
        lines.append("color #1 light gray")
        lines.append("surface #1")
        lines.append("transparency #1 70")
        lines.append("")

    # Load hotspot BILD from 05d if it exists
    if bild_path and bild_path.exists():
        lines.append(f"open {bild_path}")
        lines.append("")

    # Load consensus poses for top molecules (by hotspot membership)
    model = 3  # 1=receptor, 2=bild
    loaded = 0
    max_models = 30  # prevent overload

    if hotspot_map:
        # Collect molecules per hotspot, pick top by score
        for hs in hotspot_map.get("hotspots", [])[:5]:
            hs_id = hs["hotspot_id"]
            color = HOTSPOT_COLORS[hs_id % len(HOTSPOT_COLORS)]
            lines.append(f"# --- Hotspot {hs_id}: {hs['n_molecules']} molecules ---")

            for member in hs.get("members", [])[:5]:
                mol_name = member["molecule"]
                mol_exported = all_exported.get(mol_name, {})

                # Load consensus if available, otherwise first fragment
                consensus = mol_exported.get("full_consensus")
                if consensus and loaded < max_models:
                    lines.append(f"open {consensus}")
                    lines.append(f"color #{model} {color}")
                    lines.append(f"style #{model} stick")
                    model += 1
                    loaded += 1

            lines.append("")
    else:
        # No hotspot data — load all consensus poses
        for mol_name, mol_exported in sorted(all_exported.items()):
            consensus = mol_exported.get("full_consensus")
            if consensus and loaded < max_models:
                lines.append(f"open {consensus}")
                lines.append(f"style #{model} stick")
                model += 1
                loaded += 1

    lines.extend(["", "view", ""])
    return "\n".join(lines)


# =============================================================================
# PIPELINE
# =============================================================================

def run_structure_export(
        parsed_dir: str,
        cluster_dir: str,
        output_dir: str,
        score_dir: Optional[str] = None,
        hotspot_dir: Optional[str] = None,
        receptor_path: Optional[str] = None,
        max_clusters_per_fragment: int = 3,
        molecule_names: Optional[List[str]] = None,
) -> Dict[str, any]:
    """
    Export structures for all molecules, informed by hotspot analysis.

    Reads: 05a + 05b + 05c (optional) + 05d (optional)
    Saves: per-molecule mol2 + cxc, global cxc
    """
    logger.info("=" * 60)
    logger.info("05e STRUCTURE EXPORT v2.0")
    logger.info("=" * 60)

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # Load hotspot map from 05d (optional)
    hotspot_map = None
    bild_path = None
    if hotspot_dir:
        hm_path = Path(hotspot_dir) / "hotspot_map.json"
        if hm_path.exists():
            with open(hm_path) as fh:
                hotspot_map = json.load(fh)
            logger.info(f"  Loaded {len(hotspot_map.get('hotspots', []))} hotspots from 05d")
        bp = Path(hotspot_dir) / "hotspots.bild"
        if bp.exists():
            bild_path = bp

    all_mols = load_all_molecules(parsed_dir, molecule_names)
    if not all_mols:
        logger.error("No parsed molecules found")
        return {"success": False, "error": "No data"}

    all_exported = {}
    n_ok = 0
    total_files = 0

    for name in sorted(all_mols.keys()):
        mol_data = all_mols[name]
        cluster_result = load_cluster_result(name, cluster_dir)
        if cluster_result is None:
            continue

        # Load score data (optional)
        score_data = None
        if score_dir:
            sp = Path(score_dir) / f"{name}_scores.json"
            if sp.exists():
                try:
                    with open(sp) as fh:
                        score_data = json.load(fh)
                except Exception:
                    pass

        try:
            mol_dir = out_path / name
            exported = export_molecule(
                mol_data, cluster_result, mol_dir,
                hotspot_map=hotspot_map,
                score_data=score_data,
                max_clusters_per_fragment=max_clusters_per_fragment,
            )
        except Exception as e:
            logger.warning(f"  {name}: export failed — {e}")
            continue

        # Per-molecule ChimeraX script
        cxc = _generate_mol_cxc(name, mol_dir, exported, receptor_path, hotspot_map)
        cxc_path = mol_dir / f"{name}_visualize.cxc"
        cxc_path.write_text(cxc)
        exported["cxc_script"] = str(cxc_path)

        all_exported[name] = exported
        n_ok += 1
        total_files += len(exported)

    # Global ChimeraX script
    global_cxc = _generate_global_cxc(
        out_path, all_exported, receptor_path, hotspot_map, bild_path
    )
    global_cxc_path = out_path / "master_visualization.cxc"
    global_cxc_path.write_text(global_cxc)

    logger.info(f"\n{'=' * 60}")
    logger.info("STRUCTURE EXPORT COMPLETE")
    logger.info(f"{'=' * 60}")
    logger.info(f"  Molecules: {n_ok}")
    logger.info(f"  Total files: {total_files}")
    logger.info(f"  Master script: {global_cxc_path}")

    return {
        "success": True,
        "n_molecules": n_ok,
        "total_files": total_files,
        "output_dir": str(out_path),
    }

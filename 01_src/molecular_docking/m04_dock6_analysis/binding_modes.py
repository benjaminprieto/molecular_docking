"""
Binding Modes - Core Module (04c)
====================================
Characterize DOCK6 binding modes per molecule.

Each DOCK6 pose IS a genuine binding mode (pre-clustered at 2.0Å RMSD).
This module describes each mode by its interacting residues, vdW vs ES
balance, and computes RMSD between modes.

Input:
    01c_dock6_run/{name}/{name}_scored.mol2  — scored poses (2-14 per molecule)
    04b footprint data (optional): footprint_per_molecule.csv

Output:
    binding_modes_summary.csv   — modes per molecule (n_modes, best, spread)
    binding_modes_detail.csv    — one row per mode per molecule
    {name}_modes.json           — per-molecule mode detail (in per_molecule/)

Location: 01_src/molecular_docking/m04_dock6_analysis/binding_modes.py
Project: molecular_docking
Module: 04c (DOCK6 analysis)
Version: 1.0 (2026-03-22)
"""

import json
import logging
import re
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# MOL2 COORDINATE PARSER
# =============================================================================

def _parse_mol2_poses(mol2_path: str) -> List[Dict[str, Any]]:
    """
    Parse a DOCK6 scored mol2 with multiple poses.

    Returns list of dicts with:
        name, pose_id, coords (Nx3 array), header_fields
    """
    with open(mol2_path, "r") as f:
        content = f.read()

    blocks = content.split("@<TRIPOS>MOLECULE")
    if len(blocks) < 2:
        return []

    poses = []

    for block_idx, block in enumerate(blocks[1:], 0):
        full_block = "@<TRIPOS>MOLECULE" + block
        lines = block.strip().split("\n")
        mol_name = lines[0].strip() if lines else f"pose_{block_idx}"

        # Parse header fields
        header = {}
        for line in lines:
            m = re.match(r"##########\s+(\S+)\s*:\s*(.*)", line)
            if m:
                key, val = m.group(1), m.group(2).strip()
                try:
                    header[key] = float(val)
                except ValueError:
                    header[key] = val

        # Parse atom coordinates from ATOM section
        coords = []
        in_atoms = False
        for line in full_block.split("\n"):
            if line.startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                continue
            if line.startswith("@<TRIPOS>") and in_atoms:
                break
            if in_atoms and line.strip():
                parts = line.split()
                if len(parts) >= 5:
                    try:
                        x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                        coords.append([x, y, z])
                    except ValueError:
                        continue

        poses.append({
            "name": mol_name,
            "pose_id": block_idx,
            "coords": np.array(coords) if coords else np.empty((0, 3)),
            "header_fields": header,
            "grid_score": header.get("Grid_Score"),
            "grid_vdw": header.get("Grid_vdw_energy"),
            "grid_es": header.get("Grid_es_energy"),
        })

    return poses


def _compute_rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """Compute RMSD between two coordinate arrays of same shape."""
    if coords1.shape != coords2.shape or len(coords1) == 0:
        return float("inf")
    diff = coords1 - coords2
    return float(np.sqrt(np.mean(np.sum(diff**2, axis=1))))


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def run_binding_modes(
        docking_dir: Union[str, Path],
        output_dir: Union[str, Path],
        footprint_csv: Optional[Union[str, Path]] = None,
) -> Dict[str, Any]:
    """
    Characterize binding modes for each molecule.

    Args:
        docking_dir: Path to 01c_dock6_run output
        output_dir:  Output directory
        footprint_csv: Optional path to footprint_per_molecule.csv (from 04b)

    Returns:
        Dict with: success, n_molecules, output paths
    """
    docking_dir = Path(docking_dir)
    output_dir = Path(output_dir)
    per_mol_dir = output_dir / "per_molecule"
    output_dir.mkdir(parents=True, exist_ok=True)
    per_mol_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("  04c DOCK6 Binding Modes v1.0")
    logger.info("=" * 60)

    # Load footprint data if available
    df_fp = None
    if footprint_csv and Path(footprint_csv).exists():
        df_fp = pd.read_csv(footprint_csv)
        logger.info(f"  Footprint data: {len(df_fp)} rows")

    # --- Find scored mol2 files ---
    mol_dirs = sorted([
        d for d in docking_dir.iterdir()
        if d.is_dir() and not d.name.startswith(".")
    ])

    summary_rows = []
    detail_rows = []
    n_parsed = 0

    for d in mol_dirs:
        name = d.name
        scored_mol2 = d / f"{name}_scored.mol2"
        if not scored_mol2.exists():
            continue

        poses = _parse_mol2_poses(str(scored_mol2))
        if not poses:
            continue

        n_poses = len(poses)
        n_parsed += 1

        # --- Compute pairwise RMSD ---
        rmsd_matrix = np.zeros((n_poses, n_poses))
        for i in range(n_poses):
            for j in range(i + 1, n_poses):
                rmsd = _compute_rmsd(poses[i]["coords"], poses[j]["coords"])
                rmsd_matrix[i, j] = rmsd
                rmsd_matrix[j, i] = rmsd

        # --- Characterize each mode ---
        scores = [p["grid_score"] for p in poses if p["grid_score"] is not None]
        best_idx = int(np.argmin(scores)) if scores else 0

        mode_details = []
        for pidx, pose in enumerate(poses):
            gs = pose.get("grid_score")
            gvdw = pose.get("grid_vdw")
            ges = pose.get("grid_es")

            # vdW/ES balance for this mode
            if gvdw is not None and ges is not None and gs is not None and gs != 0:
                vdw_frac = gvdw / gs
                es_frac = ges / gs
            else:
                vdw_frac = None
                es_frac = None

            # Average RMSD to other modes
            if n_poses > 1:
                avg_rmsd = float(np.mean([rmsd_matrix[pidx, j]
                                          for j in range(n_poses) if j != pidx]))
                min_rmsd = float(np.min([rmsd_matrix[pidx, j]
                                         for j in range(n_poses) if j != pidx]))
            else:
                avg_rmsd = 0.0
                min_rmsd = 0.0

            mode_info = {
                "Name": name,
                "pose_id": pidx,
                "is_best": pidx == best_idx,
                "Grid_Score": gs,
                "Grid_vdw_energy": gvdw,
                "Grid_es_energy": ges,
                "vdw_fraction": round(vdw_frac, 3) if vdw_frac is not None else None,
                "es_fraction": round(es_frac, 3) if es_frac is not None else None,
                "avg_rmsd_to_others": round(avg_rmsd, 2),
                "min_rmsd_to_others": round(min_rmsd, 2),
                "n_atoms": len(pose["coords"]),
            }
            mode_details.append(mode_info)
            detail_rows.append(mode_info)

        # --- Summary per molecule ---
        scores_arr = np.array([s for s in scores if s is not None])
        summary_rows.append({
            "Name": name,
            "n_modes": n_poses,
            "best_Grid_Score": float(scores_arr.min()) if len(scores_arr) > 0 else None,
            "worst_Grid_Score": float(scores_arr.max()) if len(scores_arr) > 0 else None,
            "mean_Grid_Score": round(float(scores_arr.mean()), 3) if len(scores_arr) > 0 else None,
            "std_Grid_Score": round(float(scores_arr.std()), 3) if len(scores_arr) > 1 else 0.0,
            "score_spread": round(float(scores_arr.max() - scores_arr.min()), 3) if len(scores_arr) > 1 else 0.0,
            "mean_pairwise_rmsd": round(float(rmsd_matrix[np.triu_indices(n_poses, k=1)].mean()), 2)
                                  if n_poses > 1 else 0.0,
            "max_pairwise_rmsd": round(float(rmsd_matrix.max()), 2) if n_poses > 1 else 0.0,
        })

        # --- Save per-molecule JSON ---
        mol_json = per_mol_dir / f"{name}_modes.json"
        with open(mol_json, "w") as f:
            # Serialize RMSD matrix as list of lists
            json.dump({
                "name": name,
                "n_modes": n_poses,
                "modes": mode_details,
                "rmsd_matrix": rmsd_matrix.round(2).tolist(),
            }, f, indent=2, default=str)

    logger.info(f"  Parsed: {n_parsed} molecules")

    if not summary_rows:
        return {"success": False, "error": "No scored mol2 files found"}

    # --- Save CSVs ---
    df_summary = pd.DataFrame(summary_rows)
    df_summary.sort_values("best_Grid_Score", ascending=True, inplace=True)
    df_summary.reset_index(drop=True, inplace=True)

    summary_csv = output_dir / "binding_modes_summary.csv"
    df_summary.to_csv(summary_csv, index=False, encoding="utf-8")
    logger.info(f"  Saved: {summary_csv}")

    df_detail = pd.DataFrame(detail_rows)
    detail_csv = output_dir / "binding_modes_detail.csv"
    df_detail.to_csv(detail_csv, index=False, encoding="utf-8")
    logger.info(f"  Saved: {detail_csv}")

    # --- Log summary ---
    modes = df_summary["n_modes"]
    logger.info(f"  Modes per molecule: {modes.min()}-{modes.max()} "
                f"(mean {modes.mean():.1f})")
    logger.info(f"  Total modes: {modes.sum()}")
    logger.info("=" * 60)

    return {
        "success": True,
        "n_molecules": n_parsed,
        "total_modes": int(modes.sum()),
        "binding_modes_summary_csv": str(summary_csv),
        "binding_modes_detail_csv": str(detail_csv),
        "per_molecule_dir": str(per_mol_dir),
        "output_dir": str(output_dir),
    }

"""
Contact Mapping - Core Module (04d)
======================================
Distance-based contacts cross-referenced with footprint energies.

Key insight: a residue can be geometrically close (contact) but contribute
nothing energetically, or far but contribute via long-range electrostatics.
This module captures both perspectives.

Input:
    01c_dock6_run/{name}/{name}_scored.mol2  — pose coordinates
    00b receptor: rec_charged.mol2 or receptor PDB
    04b footprint: footprint_per_molecule.csv (optional)

Output:
    contact_summary.csv         — same format as GNINA 05f
    contact_vs_footprint.csv    — contact distance vs energy contribution

Location: 01_src/molecular_docking/m04_dock6_analysis/contact_mapping.py
Project: molecular_docking
Module: 04d (DOCK6 analysis)
Version: 1.0 (2026-03-22)
"""

import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# PDB/MOL2 PARSING
# =============================================================================

def _read_receptor_atoms(receptor_path: str) -> List[Dict[str, Any]]:
    """
    Read receptor atom coordinates from PDB or mol2.
    Returns list of dicts: {atom_name, res_name, res_num, chain, x, y, z}
    """
    path = Path(receptor_path)
    atoms = []

    if path.suffix.lower() in (".pdb", ".ent"):
        with open(path) as f:
            for line in f:
                if not line.startswith("ATOM"):
                    continue
                try:
                    atoms.append({
                        "atom_name": line[12:16].strip(),
                        "res_name": line[17:20].strip(),
                        "res_num": int(line[22:26].strip()),
                        "chain": line[21].strip() or "A",
                        "x": float(line[30:38]),
                        "y": float(line[38:46]),
                        "z": float(line[46:54]),
                    })
                except (ValueError, IndexError):
                    continue

    elif path.suffix.lower() == ".mol2":
        with open(path) as f:
            content = f.read()

        # Parse ATOM section
        in_atoms = False
        for line in content.split("\n"):
            if line.startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                continue
            if line.startswith("@<TRIPOS>") and in_atoms:
                break
            if in_atoms and line.strip():
                parts = line.split()
                if len(parts) >= 9:
                    try:
                        # mol2 format: atom_id atom_name x y z atom_type res_id res_name charge
                        subst = parts[7] if len(parts) > 7 else ""
                        # Parse residue info from substructure name (e.g., "ALA123")
                        m = re.match(r"([A-Z]{1,4})(\d+)", subst)
                        res_name = m.group(1) if m else subst[:3]
                        res_num = int(m.group(2)) if m else 0

                        atoms.append({
                            "atom_name": parts[1],
                            "res_name": res_name,
                            "res_num": res_num,
                            "chain": "A",  # mol2 doesn't have chain directly
                            "x": float(parts[2]),
                            "y": float(parts[3]),
                            "z": float(parts[4]),
                        })
                    except (ValueError, IndexError):
                        continue
    else:
        logger.error(f"  Unsupported receptor format: {path.suffix}")

    return atoms


def _parse_ligand_coords(mol2_path: str, pose_idx: int = 0) -> np.ndarray:
    """Extract atom coordinates from a specific pose in a scored mol2."""
    with open(mol2_path) as f:
        content = f.read()

    blocks = content.split("@<TRIPOS>MOLECULE")
    if pose_idx + 1 >= len(blocks):
        return np.empty((0, 3))

    block = "@<TRIPOS>MOLECULE" + blocks[pose_idx + 1]
    coords = []
    in_atoms = False
    for line in block.split("\n"):
        if line.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            continue
        if line.startswith("@<TRIPOS>") and in_atoms:
            break
        if in_atoms and line.strip():
            parts = line.split()
            if len(parts) >= 5:
                try:
                    coords.append([float(parts[2]), float(parts[3]), float(parts[4])])
                except ValueError:
                    continue

    return np.array(coords) if coords else np.empty((0, 3))


# =============================================================================
# CONTACT COMPUTATION
# =============================================================================

def _compute_contacts(
        ligand_coords: np.ndarray,
        receptor_atoms: List[Dict[str, Any]],
        cutoff: float = 4.5,
) -> List[Dict[str, Any]]:
    """
    Find receptor residues within cutoff distance of any ligand atom.

    Returns list of dicts per residue:
        {residue_id, residue_name, residue_number, chain,
         min_distance, n_contacts, closest_atom}
    """
    if len(ligand_coords) == 0:
        return []

    cutoff_sq = cutoff ** 2
    residue_contacts = {}

    for atom in receptor_atoms:
        ax, ay, az = atom["x"], atom["y"], atom["z"]
        res_id = f"{atom['res_name']}{atom['res_num']}.{atom['chain']}"

        # Vectorized distance to all ligand atoms
        diff = ligand_coords - np.array([ax, ay, az])
        dsq = np.sum(diff ** 2, axis=1)
        min_dsq = float(dsq.min())

        if min_dsq <= cutoff_sq:
            d = np.sqrt(min_dsq)
            if res_id not in residue_contacts:
                residue_contacts[res_id] = {
                    "residue_id": res_id,
                    "residue_name": atom["res_name"],
                    "residue_number": atom["res_num"],
                    "chain": atom["chain"],
                    "min_distance": d,
                    "n_contacts": 0,
                    "closest_atom": atom["atom_name"],
                }
            residue_contacts[res_id]["n_contacts"] += 1
            if d < residue_contacts[res_id]["min_distance"]:
                residue_contacts[res_id]["min_distance"] = d
                residue_contacts[res_id]["closest_atom"] = atom["atom_name"]

    return sorted(residue_contacts.values(), key=lambda r: r["min_distance"])


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def run_contact_mapping(
        docking_dir: Union[str, Path],
        receptor_path: Union[str, Path],
        output_dir: Union[str, Path],
        footprint_csv: Optional[Union[str, Path]] = None,
        contact_cutoff: float = 4.5,
        best_pose_only: bool = True,
        all_poses_csv: Optional[Union[str, Path]] = None,
) -> Dict[str, Any]:
    """
    Map contacts for all DOCK6 docked molecules.

    Args:
        docking_dir:    Path to 01c_dock6_run output
        receptor_path:  Path to receptor PDB or mol2
        output_dir:     Output directory
        footprint_csv:  Path to footprint_per_molecule.csv (from 04b, optional)
        contact_cutoff: Distance cutoff in Angstroms
        best_pose_only: Only analyze best pose per molecule
        all_poses_csv:  Path to dock6_all_poses.csv (to identify best pose)

    Returns:
        Dict with: success, n_molecules, output paths
    """
    docking_dir = Path(docking_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("  04d DOCK6 Contact Mapping v1.0")
    logger.info("=" * 60)
    logger.info(f"  Contact cutoff: {contact_cutoff} Å")

    # --- Load receptor ---
    receptor_atoms = _read_receptor_atoms(str(receptor_path))
    logger.info(f"  Receptor atoms: {len(receptor_atoms)}")

    if not receptor_atoms:
        return {"success": False, "error": f"No atoms read from {receptor_path}"}

    # --- Load best pose info ---
    best_pose_ids = {}  # name -> pose_id (0-indexed)
    if all_poses_csv and Path(all_poses_csv).exists():
        df_poses = pd.read_csv(all_poses_csv)
        for name, grp in df_poses.groupby("Name"):
            best_idx = grp["Grid_Score"].idxmin() if "Grid_Score" in grp.columns else grp.index[0]
            # pose_id within molecule (0-indexed)
            best_pose_ids[name] = int(grp.index.get_loc(best_idx))

    # --- Load footprint data ---
    df_fp = None
    if footprint_csv and Path(footprint_csv).exists():
        df_fp = pd.read_csv(footprint_csv)
        logger.info(f"  Footprint data: {len(df_fp)} rows")

    # --- Process each molecule ---
    mol_dirs = sorted([
        d for d in docking_dir.iterdir()
        if d.is_dir() and not d.name.startswith(".")
    ])

    contact_rows = []
    cross_rows = []
    n_parsed = 0

    for d in mol_dirs:
        name = d.name
        scored_mol2 = d / f"{name}_scored.mol2"
        if not scored_mol2.exists():
            continue

        # Determine which pose to analyze
        pose_idx = best_pose_ids.get(name, 0) if best_pose_only else 0
        coords = _parse_ligand_coords(str(scored_mol2), pose_idx)
        if len(coords) == 0:
            continue

        # Compute contacts
        contacts = _compute_contacts(coords, receptor_atoms, contact_cutoff)
        n_parsed += 1

        for c in contacts:
            contact_rows.append({
                "Name": name,
                "pose_id": pose_idx,
                **c,
            })

        # Cross-reference with footprint
        if df_fp is not None:
            mol_fp = df_fp[df_fp["Name"] == name]
            contact_ids = {c["residue_id"] for c in contacts}
            fp_ids = set(mol_fp["residue_id"].unique()) if "residue_id" in mol_fp.columns else set()

            # Residues that are contacts AND have footprint energy
            for c in contacts:
                rid = c["residue_id"]
                fp_row = mol_fp[mol_fp["residue_id"] == rid] if "residue_id" in mol_fp.columns else pd.DataFrame()
                fp_vdw = float(fp_row["vdw"].iloc[0]) if len(fp_row) > 0 and "vdw" in fp_row.columns else None
                fp_es = float(fp_row["es"].iloc[0]) if len(fp_row) > 0 and "es" in fp_row.columns else None
                fp_total = float(fp_row["total"].iloc[0]) if len(fp_row) > 0 and "total" in fp_row.columns else None

                cross_rows.append({
                    "Name": name,
                    "residue_id": rid,
                    "residue_name": c["residue_name"],
                    "residue_number": c["residue_number"],
                    "chain": c["chain"],
                    "min_distance": round(c["min_distance"], 2),
                    "n_contacts": c["n_contacts"],
                    "fp_vdw": fp_vdw,
                    "fp_es": fp_es,
                    "fp_total": fp_total,
                    "is_contact": True,
                    "has_energy": fp_total is not None and fp_total < -0.1,
                    "category": _categorize_contact(c["min_distance"],
                                                     fp_total),
                })

            # Residues with footprint energy but NOT geometric contacts
            for _, fp_r in mol_fp.iterrows():
                if "residue_id" not in fp_r:
                    continue
                rid = fp_r["residue_id"]
                if rid in contact_ids:
                    continue
                if "total" in fp_r and fp_r["total"] < -0.5:
                    cross_rows.append({
                        "Name": name,
                        "residue_id": rid,
                        "residue_name": fp_r.get("residue_name", ""),
                        "residue_number": int(fp_r.get("residue_number", 0)),
                        "chain": fp_r.get("chain", ""),
                        "min_distance": None,  # not a geometric contact
                        "n_contacts": 0,
                        "fp_vdw": fp_r.get("vdw"),
                        "fp_es": fp_r.get("es"),
                        "fp_total": fp_r.get("total"),
                        "is_contact": False,
                        "has_energy": True,
                        "category": "long_range_ES",
                    })

    logger.info(f"  Parsed: {n_parsed} molecules")

    if not contact_rows:
        return {"success": False, "error": "No contacts found"}

    # --- Save contact_summary.csv ---
    df_contacts = pd.DataFrame(contact_rows)

    # Per-residue consensus: how many molecules contact each residue?
    res_consensus = []
    for rid, grp in df_contacts.groupby("residue_id"):
        n_mol = grp["Name"].nunique()
        res_consensus.append({
            "residue_id": rid,
            "residue_name": grp.iloc[0]["residue_name"],
            "residue_number": grp.iloc[0]["residue_number"],
            "chain": grp.iloc[0]["chain"],
            "n_molecules": n_mol,
            "frac_molecules": round(n_mol / n_parsed, 3),
            "mean_min_distance": round(grp["min_distance"].mean(), 2),
            "mean_n_contacts": round(grp["n_contacts"].mean(), 1),
        })

    df_consensus = pd.DataFrame(res_consensus)
    df_consensus.sort_values("n_molecules", ascending=False, inplace=True)

    summary_csv = output_dir / "contact_summary.csv"
    df_consensus.to_csv(summary_csv, index=False, encoding="utf-8")
    logger.info(f"  Saved: {summary_csv} ({len(df_consensus)} residues)")

    # --- Save contact_vs_footprint.csv ---
    if cross_rows:
        df_cross = pd.DataFrame(cross_rows)
        cross_csv = output_dir / "contact_vs_footprint.csv"
        df_cross.to_csv(cross_csv, index=False, encoding="utf-8")
        logger.info(f"  Saved: {cross_csv} ({len(df_cross)} entries)")

        # Statistics
        n_contact_energy = sum(1 for _, r in df_cross.iterrows()
                               if r["is_contact"] and r["has_energy"])
        n_contact_no_energy = sum(1 for _, r in df_cross.iterrows()
                                  if r["is_contact"] and not r["has_energy"])
        n_long_range = sum(1 for _, r in df_cross.iterrows()
                           if not r["is_contact"] and r["has_energy"])
        logger.info(f"  Contact + energy:    {n_contact_energy}")
        logger.info(f"  Contact, no energy:  {n_contact_no_energy}")
        logger.info(f"  Long-range ES:       {n_long_range}")
    else:
        cross_csv = None

    # --- Save per-molecule contacts ---
    detail_csv = output_dir / "contact_detail.csv"
    df_contacts.to_csv(detail_csv, index=False, encoding="utf-8")

    logger.info("=" * 60)

    return {
        "success": True,
        "n_molecules": n_parsed,
        "n_residues": len(df_consensus),
        "contact_summary_csv": str(summary_csv),
        "contact_vs_footprint_csv": str(cross_csv) if cross_csv else None,
        "contact_detail_csv": str(detail_csv),
        "output_dir": str(output_dir),
    }


def _categorize_contact(distance: float, energy: Optional[float]) -> str:
    """Categorize a contact by distance vs energy."""
    if energy is None:
        return "contact_only"
    if distance is not None and distance < 4.5 and energy < -0.5:
        return "contact_and_energy"
    if distance is not None and distance < 4.5 and energy >= -0.5:
        return "contact_no_energy"
    if (distance is None or distance >= 4.5) and energy < -0.5:
        return "long_range_ES"
    return "weak"

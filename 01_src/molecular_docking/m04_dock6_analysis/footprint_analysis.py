"""
Footprint Analysis - Core Module (04b)
=========================================
Per-residue vdW + ES energy decomposition via DOCK6 footprint scoring.

Pipeline:
    Step 1: Build residue mapping (mol2 sequential → PDB original numbering)
    Step 2: Parse per-residue vdW + ES from footprint TXT files
    Step 3: Cross-molecule consensus (which residues always contribute)
    Step 4: Compare each molecule's footprint vs reference (UDX)

DOCK6 6.13 footprint output format:
    - {name}_fps_scored.mol2           → summary scores per pose
    - {name}_fps_footprint_scored.txt  → per-residue tabular data
    - {name}_fps_hbond_scored.txt      → H-bond details

Residue numbering:
    ChimeraX mol2 renumbers residues sequentially (1..N).
    DOCK6 footprint inherits this sequential numbering.
    This module remaps to PDB original numbering (e.g., 141→392)
    so that footprint and contact results share the same residue IDs.

Input:
    01d_footprint_rescore/{name}/{name}_fps_footprint_scored.txt
    00b_receptor_preparation/rec_charged.mol2  (sequential numbering)
    00b_receptor_preparation/rec_noH.pdb       (PDB numbering)

Output:
    footprint_per_molecule.csv    — residue × molecule energy matrix (PDB numbering)
    residue_consensus.csv         — which residues always contribute
    vs_reference_comparison.csv   — delta vdW/ES vs reference per residue
    pharmacophore_residues.json   — residues contacted by >80% of molecules
    molecule_footprint_summary.csv — one row per molecule with totals
    residue_mapping.csv           — sequential→PDB mapping for reference

Location: 01_src/molecular_docking/m04_dock6_analysis/footprint_analysis.py
Project: molecular_docking
Module: 04b (DOCK6 analysis)
Version: 3.0 (2026-03-23) — adds sequential→PDB residue remapping

Reference: Balius et al. J Chem Inf Model 2011, 51(8):1942-56
"""

import json
import logging
import re
from pathlib import Path
from typing import Dict, List, Any, Optional, Union, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# RESIDUE MAPPING: mol2 sequential → PDB original
# =============================================================================

def build_residue_mapping(
        receptor_mol2: str,
        receptor_pdb: str,
) -> Dict[str, str]:
    """
    Build mapping from mol2 sequential numbering to PDB original numbering.

    ChimeraX mol2 SUBSTRUCTURE section lists residues sequentially (1..N).
    PDB CA atoms list residues with original numbering + chain.
    Both have the same number of residues in the same order.

    Returns:
        Dict mapping "RESseq" → "RES_pdbnum.chain"
        e.g., {"TRP141" → "TRP392.A", "HIS84" → "HIS335.A"}
    """
    # --- Parse mol2 SUBSTRUCTURE ---
    mol2_residues = []  # list of (seq_id, resname)
    with open(receptor_mol2, "r") as f:
        in_substructure = False
        for line in f:
            if "@<TRIPOS>SUBSTRUCTURE" in line:
                in_substructure = True
                continue
            if in_substructure:
                if line.startswith("@"):
                    break
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        seq_id = int(parts[0])
                        resname = parts[1]
                        mol2_residues.append((seq_id, resname))
                    except (ValueError, IndexError):
                        continue

    # --- Parse PDB CA atoms ---
    pdb_residues = []  # list of (resname, resid, chain)
    with open(receptor_pdb, "r") as f:
        for line in f:
            if (line.startswith("ATOM") or line.startswith("HETATM")) and " CA " in line:
                resname = line[17:20].strip()
                chain = line[21].strip() or "A"
                try:
                    resid = int(line[22:26].strip())
                except ValueError:
                    continue
                pdb_residues.append((resname, resid, chain))

    # --- Build mapping ---
    if len(mol2_residues) != len(pdb_residues):
        logger.warning(f"  Residue count mismatch: mol2={len(mol2_residues)}, PDB={len(pdb_residues)}")
        logger.warning("  Falling back to sequential numbering (no remapping)")
        return {}

    mapping = {}
    mismatches = 0
    for (seq_id, mol2_name), (pdb_name, pdb_resid, chain) in zip(mol2_residues, pdb_residues):
        # Verify residue names match
        if mol2_name[:3].upper() != pdb_name[:3].upper():
            mismatches += 1
            if mismatches <= 3:
                logger.warning(f"  Residue name mismatch at seq={seq_id}: "
                               f"mol2={mol2_name}, PDB={pdb_name}{pdb_resid}.{chain}")

        seq_key = f"{mol2_name}{seq_id}"
        pdb_key = f"{pdb_name}{pdb_resid}.{chain}"
        mapping[seq_key] = pdb_key

    if mismatches > 3:
        logger.warning(f"  ... and {mismatches - 3} more mismatches")
    if mismatches > len(mol2_residues) * 0.1:
        logger.error(f"  Too many mismatches ({mismatches}/{len(mol2_residues)}), disabling remapping")
        return {}

    logger.info(f"  Residue mapping: {len(mapping)} residues (mol2 sequential → PDB original)")
    return mapping


# =============================================================================
# FOOTPRINT TXT PARSER (per-residue tabular data)
# =============================================================================

def parse_footprint_txt(
        txt_path: str,
        residue_mapping: Optional[Dict[str, str]] = None,
) -> List[Dict[str, Any]]:
    """
    Parse a DOCK6 footprint TXT file (*_fps_footprint_scored.txt).

    Applies residue_mapping to convert sequential→PDB numbering if provided.

    Returns:
        List of pose dicts, each with:
            - pose_id (int)
            - fps_score, fps_vdw_energy, fps_es_energy (float or None)
            - residue_footprint: list of per-residue dicts
    """
    with open(txt_path, "r") as f:
        content = f.read()

    # Split into pose blocks by the "### Molecule:" separator
    blocks = re.split(r"#{20,}\s*\n\s*###\s+Molecule:", content)

    poses = []

    for block_idx, block in enumerate(blocks):
        # Skip preamble before first molecule
        if block_idx == 0 and "resname" not in block:
            continue

        lines = block.strip().split("\n")

        # --- Parse summary fields from ## lines ---
        fps_score = None
        fps_vdw_energy = None
        fps_es_energy = None
        fps_vdw_es_energy = None
        fps_num_hbond = None

        for line in lines:
            line_s = line.strip()
            if "Footprint_Similarity_Score:" in line_s:
                try:
                    fps_score = float(line_s.split(":")[-1].strip())
                except ValueError:
                    pass
            elif "FPS_vdw_energy:" in line_s:
                try:
                    fps_vdw_energy = float(line_s.split(":")[-1].strip())
                except ValueError:
                    pass
            elif "FPS_es_energy:" in line_s:
                try:
                    fps_es_energy = float(line_s.split(":")[-1].strip())
                except ValueError:
                    pass
            elif "FPS_vdw+es_energy:" in line_s:
                try:
                    fps_vdw_es_energy = float(line_s.split(":")[-1].strip())
                except ValueError:
                    pass
            elif "FPS_num_hbond:" in line_s:
                try:
                    fps_num_hbond = int(line_s.split(":")[-1].strip())
                except ValueError:
                    pass

        # --- Parse per-residue tabular data ---
        residue_footprint = []
        in_table = False

        for line in lines:
            line_s = line.strip()

            # Detect header row
            if line_s.startswith("resname") and "resid" in line_s:
                in_table = True
                continue

            # Detect end of table (empty line or new ## block)
            if in_table and (not line_s or line_s.startswith("#")):
                in_table = False
                continue

            if in_table:
                parts = line_s.split()
                if len(parts) >= 8:
                    try:
                        resname = parts[0]
                        resid = int(parts[1])
                        vdw_ref = float(parts[2])
                        es_ref = float(parts[3])
                        hb_ref = int(float(parts[4]))
                        vdw_pose = float(parts[5])
                        es_pose = float(parts[6])
                        hb_pose = int(float(parts[7]))

                        # Sequential key (from DOCK6 output)
                        seq_key = f"{resname}{resid}"

                        # Remap to PDB numbering if mapping available
                        if residue_mapping and seq_key in residue_mapping:
                            residue_id = residue_mapping[seq_key]
                            # Parse PDB residue_id back to components
                            m_pdb = re.match(r"([A-Z]{1,4})(\d+)\.(\w+)", residue_id)
                            if m_pdb:
                                resname_out = m_pdb.group(1)
                                resid_out = int(m_pdb.group(2))
                                chain_out = m_pdb.group(3)
                            else:
                                resname_out = resname
                                resid_out = resid
                                chain_out = "A"
                                residue_id = f"{resname}{resid}.A"
                        else:
                            resname_out = resname
                            resid_out = resid
                            chain_out = "A"
                            residue_id = f"{resname}{resid}.A"

                        total_pose = vdw_pose + es_pose
                        total_ref = vdw_ref + es_ref

                        residue_footprint.append({
                            "residue_id": residue_id,
                            "residue_name": resname_out,
                            "residue_number": resid_out,
                            "chain": chain_out,
                            "vdw": round(vdw_pose, 6),
                            "es": round(es_pose, 6),
                            "total": round(total_pose, 6),
                            "ref_vdw": round(vdw_ref, 6),
                            "ref_es": round(es_ref, 6),
                            "ref_total": round(total_ref, 6),
                            "delta_vdw": round(vdw_pose - vdw_ref, 6),
                            "delta_es": round(es_pose - es_ref, 6),
                            "delta_total": round(total_pose - total_ref, 6),
                            "hb_pose": hb_pose,
                            "hb_ref": hb_ref,
                        })
                    except (ValueError, IndexError) as e:
                        logger.debug(f"  Skipping malformed line: {line_s} ({e})")
                        continue

        pose_id = block_idx if block_idx > 0 else 0
        poses.append({
            "pose_id": pose_id,
            "fps_score": fps_score,
            "fps_vdw_energy": fps_vdw_energy,
            "fps_es_energy": fps_es_energy,
            "fps_vdw_es_energy": fps_vdw_es_energy,
            "fps_num_hbond": fps_num_hbond,
            "residue_footprint": residue_footprint,
            "n_residues": len(residue_footprint),
        })

    return poses


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def run_footprint_analysis(
        footprint_dir: Union[str, Path],
        output_dir: Union[str, Path],
        receptor_mol2: Optional[str] = None,
        receptor_pdb: Optional[str] = None,
        pharmacophore_threshold: float = 0.8,
        energy_cutoff: float = -0.5,
        best_pose_only: bool = True,
) -> Dict[str, Any]:
    """
    Analyze DOCK6 footprint re-scoring results.

    Args:
        footprint_dir:  Directory with {name}/{name}_fps_footprint_scored.txt
        output_dir:     Output directory for analysis files
        receptor_mol2:  Path to rec_charged.mol2 (sequential numbering, for mapping)
        receptor_pdb:   Path to rec_noH.pdb (PDB original numbering, for mapping)
        pharmacophore_threshold: Fraction of molecules that must contact a residue
        energy_cutoff:  Minimum total energy for a residue to count as "contributing"
        best_pose_only: If True, analyze only best pose per molecule

    Returns:
        Dict with: success, n_molecules, output paths
    """
    footprint_dir = Path(footprint_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("  04b DOCK6 Footprint Analysis v3.0")
    logger.info("=" * 60)

    # --- Build residue mapping ---
    residue_mapping = {}
    if receptor_mol2 and receptor_pdb and Path(receptor_mol2).exists() and Path(receptor_pdb).exists():
        logger.info(f"  Building residue mapping (mol2 → PDB)...")
        residue_mapping = build_residue_mapping(receptor_mol2, receptor_pdb)
        if residue_mapping:
            # Save mapping for reference
            mapping_csv = output_dir / "residue_mapping.csv"
            mapping_rows = []
            for seq_key, pdb_key in sorted(residue_mapping.items(),
                                            key=lambda x: int(re.search(r'\d+', x[0]).group())):
                mapping_rows.append({"mol2_sequential": seq_key, "pdb_original": pdb_key})
            pd.DataFrame(mapping_rows).to_csv(mapping_csv, index=False)
            logger.info(f"  Saved: {mapping_csv}")
    else:
        logger.warning("  No receptor mol2/PDB provided — using sequential numbering")
        logger.warning("  (Residue IDs will NOT match contact mapping)")

    # --- Find molecule directories ---
    mol_dirs = sorted([
        d for d in footprint_dir.iterdir()
        if d.is_dir() and not d.name.startswith(".")
    ])

    all_rows = []       # flat table: Name × residue × energy
    mol_summaries = []  # one row per molecule
    n_parsed = 0
    n_failed = 0

    for d in mol_dirs:
        name = d.name

        # --- Find footprint TXT file ---
        fps_txt = d / f"{name}_fps_footprint_scored.txt"
        if not fps_txt.exists():
            txt_candidates = list(d.glob("*_fps_footprint_scored.txt"))
            if txt_candidates:
                fps_txt = txt_candidates[0]
            else:
                logger.debug(f"  No footprint TXT found for {name}, skipping")
                n_failed += 1
                continue

        poses = parse_footprint_txt(str(fps_txt), residue_mapping=residue_mapping)
        if not poses:
            logger.warning(f"  No poses parsed from TXT: {name}")
            n_failed += 1
            continue

        # Filter out empty poses
        poses = [p for p in poses if p["residue_footprint"]]
        if not poses:
            logger.warning(f"  No residue data in any pose: {name}")
            n_failed += 1
            continue

        # Select best pose (lowest FPS_Score = most similar to reference)
        if best_pose_only and len(poses) > 1:
            scored = [p for p in poses if p["fps_score"] is not None]
            if scored:
                poses = [min(scored, key=lambda p: p["fps_score"])]
            else:
                poses = [poses[0]]

        for pose in poses:
            fps_score = pose.get("fps_score")
            fps_vdw = pose.get("fps_vdw_energy")
            fps_es = pose.get("fps_es_energy")

            for res in pose["residue_footprint"]:
                all_rows.append({
                    "Name": name,
                    "pose_id": pose["pose_id"],
                    "fps_score": fps_score,
                    **res,
                })

            # Summary per molecule
            fp = pose["residue_footprint"]
            total_vdw = sum(r["vdw"] for r in fp)
            total_es = sum(r["es"] for r in fp)
            n_contributing = sum(1 for r in fp if r["total"] < energy_cutoff)

            mol_summaries.append({
                "Name": name,
                "pose_id": pose["pose_id"],
                "fps_score": fps_score,
                "fps_vdw_energy": fps_vdw,
                "fps_es_energy": fps_es,
                "fps_vdw_es_energy": pose.get("fps_vdw_es_energy"),
                "fps_num_hbond": pose.get("fps_num_hbond"),
                "total_vdw": round(total_vdw, 3),
                "total_es": round(total_es, 3),
                "total_energy": round(total_vdw + total_es, 3),
                "n_residues_total": len(fp),
                "n_residues_contributing": n_contributing,
            })

        n_parsed += 1

    logger.info(f"  Parsed: {n_parsed} molecules, {n_failed} failed")

    if not all_rows:
        logger.warning("  No footprint data found in any TXT file!")
        logger.warning("  Possible causes:")
        logger.warning("    1. Footprint re-scoring hasn't been run (run 01d first)")
        logger.warning("    2. TXT files not generated (check write_footprints=yes)")
        return {
            "success": False,
            "error": "No footprint data found. Run 01d footprint re-scoring first.",
            "n_parsed": n_parsed,
        }

    # --- Build DataFrames ---
    df_all = pd.DataFrame(all_rows)
    df_summary = pd.DataFrame(mol_summaries)

    # --- Save footprint_per_molecule.csv ---
    fps_csv = output_dir / "footprint_per_molecule.csv"
    df_all.to_csv(fps_csv, index=False, encoding="utf-8")
    logger.info(f"  Saved: {fps_csv} ({len(df_all)} rows)")

    # --- Residue consensus ---
    residue_stats = []
    all_names = df_all["Name"].nunique()

    for res_id, grp in df_all.groupby("residue_id"):
        n_mol = grp["Name"].nunique()
        n_contributing = grp[grp["total"] < energy_cutoff]["Name"].nunique()

        residue_stats.append({
            "residue_id": res_id,
            "residue_name": grp.iloc[0]["residue_name"],
            "residue_number": grp.iloc[0]["residue_number"],
            "chain": grp.iloc[0]["chain"],
            "n_molecules_present": n_mol,
            "n_molecules_contributing": n_contributing,
            "frac_present": round(n_mol / all_names, 3) if all_names > 0 else 0,
            "frac_contributing": round(n_contributing / all_names, 3) if all_names > 0 else 0,
            "mean_vdw": round(grp["vdw"].mean(), 4),
            "mean_es": round(grp["es"].mean(), 4),
            "mean_total": round(grp["total"].mean(), 4),
            "std_total": round(grp["total"].std(), 4) if len(grp) > 1 else 0.0,
            "min_total": round(grp["total"].min(), 4),
            "max_total": round(grp["total"].max(), 4),
            "ref_vdw": round(grp["ref_vdw"].iloc[0], 4) if "ref_vdw" in grp.columns else 0.0,
            "ref_es": round(grp["ref_es"].iloc[0], 4) if "ref_es" in grp.columns else 0.0,
        })

    df_consensus = pd.DataFrame(residue_stats)
    df_consensus.sort_values("mean_total", ascending=True, inplace=True)
    df_consensus.reset_index(drop=True, inplace=True)

    consensus_csv = output_dir / "residue_consensus.csv"
    df_consensus.to_csv(consensus_csv, index=False, encoding="utf-8")
    logger.info(f"  Saved: {consensus_csv} ({len(df_consensus)} residues)")

    # --- Pharmacophore residues ---
    pharma = df_consensus[
        df_consensus["frac_contributing"] >= pharmacophore_threshold
    ].copy()

    pharma_list = []
    for _, row in pharma.iterrows():
        pharma_list.append({
            "residue_id": row["residue_id"],
            "residue_name": row["residue_name"],
            "residue_number": int(row["residue_number"]),
            "chain": row["chain"],
            "frac_contributing": row["frac_contributing"],
            "mean_total": row["mean_total"],
        })

    pharma_json = output_dir / "pharmacophore_residues.json"
    with open(pharma_json, "w") as f:
        json.dump({
            "threshold": pharmacophore_threshold,
            "energy_cutoff": energy_cutoff,
            "n_molecules": all_names,
            "n_pharmacophore_residues": len(pharma_list),
            "numbering": "PDB_original" if residue_mapping else "mol2_sequential",
            "residues": pharma_list,
        }, f, indent=2)
    logger.info(f"  Saved: {pharma_json} ({len(pharma_list)} pharmacophore residues)")

    # --- vs reference comparison ---
    ref_csv = None
    if "delta_total" in df_all.columns:
        ref_comparison = []
        for name, grp in df_all.groupby("Name"):
            sig = grp[grp["delta_total"].abs() > 0.5]
            for _, row in sig.iterrows():
                ref_comparison.append({
                    "Name": name,
                    "residue_id": row["residue_id"],
                    "vdw": row["vdw"],
                    "es": row["es"],
                    "ref_vdw": row["ref_vdw"],
                    "ref_es": row["ref_es"],
                    "delta_vdw": row["delta_vdw"],
                    "delta_es": row["delta_es"],
                    "delta_total": row["delta_total"],
                })

        if ref_comparison:
            df_ref = pd.DataFrame(ref_comparison)
            ref_csv = output_dir / "vs_reference_comparison.csv"
            df_ref.to_csv(ref_csv, index=False, encoding="utf-8")
            logger.info(f"  Saved: {ref_csv} ({len(df_ref)} significant deltas)")
        else:
            logger.info("  No significant deltas vs reference found.")

    # --- Save molecule summary ---
    summary_csv = output_dir / "molecule_footprint_summary.csv"
    df_summary.to_csv(summary_csv, index=False, encoding="utf-8")
    logger.info(f"  Saved: {summary_csv}")

    # --- Log summary ---
    logger.info("")
    logger.info(f"  Molecules analyzed:      {n_parsed}")
    logger.info(f"  Total residues tracked:  {len(df_consensus)}")
    logger.info(f"  Pharmacophore residues:  {len(pharma_list)} (>{pharmacophore_threshold * 100:.0f}%)")
    logger.info(f"  Numbering:               {'PDB original' if residue_mapping else 'mol2 sequential'}")
    if len(pharma_list) > 0:
        top5 = pharma_list[:5]
        for r in top5:
            logger.info(f"    {r['residue_id']}: {r['frac_contributing'] * 100:.0f}% "
                        f"(mean={r['mean_total']:.2f} kcal/mol)")
    logger.info("=" * 60)

    return {
        "success": True,
        "n_molecules": n_parsed,
        "n_residues": len(df_consensus),
        "n_pharmacophore": len(pharma_list),
        "numbering": "PDB_original" if residue_mapping else "mol2_sequential",
        "footprint_per_molecule_csv": str(fps_csv),
        "residue_consensus_csv": str(consensus_csv),
        "pharmacophore_json": str(pharma_json),
        "vs_reference_csv": str(ref_csv) if ref_csv else None,
        "molecule_summary_csv": str(summary_csv),
        "output_dir": str(output_dir),
    }
"""
Footprint Analysis - Core Module (04b)
=========================================
Per-residue vdW + ES energy decomposition via DOCK6 footprint scoring.

Pipeline:
    Step 1: Re-score existing poses with footprint_similarity_score_primary
            (calls dock6_runner.run_footprint_rescoring — fast, no re-docking)
    Step 2: Parse per-residue vdW + ES from re-scored mol2 headers
    Step 3: Cross-molecule consensus (which residues always contribute)
    Step 4: Compare each molecule's footprint vs reference (UDX)

DOCK6 footprint output format (in scored mol2 header):
    ##########  FPS_vdw_<RES><NUM>.<CHAIN>:  <value>
    ##########  FPS_es_<RES><NUM>.<CHAIN>:   <value>
    ##########  FPS.vdw_<RES><NUM>.<CHAIN>:  <value>  (reference)
    ##########  FPS.es_<RES><NUM>.<CHAIN>:   <value>  (reference)

The exact format varies by DOCK6 version (6.9 vs 6.13). This parser
handles both formats via regex.

Input:
    01c_dock6_run/{name}/{name}_scored.mol2  — poses to re-score
    00b receptor: rec_charged.mol2
    01d best_poses: UDX best pose as reference

Output:
    footprint_per_molecule.csv    — residue × molecule energy matrix
    residue_consensus.csv         — which residues always contribute
    vs_reference_comparison.csv   — delta vdW/ES vs reference per residue
    pharmacophore_residues.json   — residues contacted by >80% of molecules

Location: 01_src/molecular_docking/m04_dock6_analysis/footprint_analysis.py
Project: molecular_docking
Module: 04b (DOCK6 analysis)
Version: 1.0 (2026-03-22)

Reference: Balius et al. J Chem Inf Model 2011, 51(8):1942-56
"""

import json
import logging
import re
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# FOOTPRINT MOL2 PARSER
# =============================================================================

# Patterns for DOCK6 footprint fields in mol2 headers.
# DOCK6 writes per-residue energies as:
#   ##########  FPS_vdw_ALA123.A:  -2.340
#   ##########  FPS_es_ALA123.A:   -0.123
# Reference (ligand used as footprint reference):
#   ##########  FPS.vdw_ALA123.A:  -1.800
#   ##########  FPS.es_ALA123.A:   -0.090
# Some DOCK6 versions use FPS.vdw vs FPS_vdw — we handle both.

_FPS_LIGAND_VDW = re.compile(
    r"##########\s+FPS[_.]vdw[_.]([A-Z]{1,4})(\d+)\.(\w+)\s*:\s*([-\d.eE+]+)"
)
_FPS_LIGAND_ES = re.compile(
    r"##########\s+FPS[_.]es[_.]([A-Z]{1,4})(\d+)\.(\w+)\s*:\s*([-\d.eE+]+)"
)
_FPS_REF_VDW = re.compile(
    r"##########\s+FPS\.vdw[_.]([A-Z]{1,4})(\d+)\.(\w+)\s*:\s*([-\d.eE+]+)"
)
_FPS_REF_ES = re.compile(
    r"##########\s+FPS\.es[_.]([A-Z]{1,4})(\d+)\.(\w+)\s*:\s*([-\d.eE+]+)"
)

# Also match the Euclidean distance score
_FPS_EUCLIDEAN = re.compile(
    r"##########\s+FPS_Score\s*:\s*([-\d.eE+]+)"
)
# Alternative: FPS_vdw_score, FPS_es_score
_FPS_VDW_SCORE = re.compile(
    r"##########\s+FPS[_.]vdw[_.]?[Ss]core\s*:\s*([-\d.eE+]+)"
)
_FPS_ES_SCORE = re.compile(
    r"##########\s+FPS[_.]es[_.]?[Ss]core\s*:\s*([-\d.eE+]+)"
)

# Standard DOCK6 header fields
_HEADER_PATTERN = re.compile(
    r"##########\s+(\S+)\s*:\s*(.*)"
)


def parse_footprint_mol2(mol2_path: str) -> List[Dict[str, Any]]:
    """
    Parse a DOCK6 footprint-scored mol2 file.

    Returns a list of poses, each with:
        - name: molecule name
        - pose_id: 0-indexed pose number
        - header_fields: dict of standard DOCK6 header fields
        - fps_score: Euclidean footprint distance (float or None)
        - residue_footprint: list of dicts per residue:
            {residue_name, residue_number, chain, residue_id,
             vdw, es, total, ref_vdw, ref_es, ref_total,
             delta_vdw, delta_es, delta_total}
    """
    with open(mol2_path, "r") as f:
        content = f.read()

    # Split into MOLECULE blocks
    blocks = content.split("@<TRIPOS>MOLECULE")
    if len(blocks) < 2:
        logger.warning(f"  No MOLECULE blocks in {mol2_path}")
        return []

    poses = []

    for block_idx, block in enumerate(blocks[1:], 0):
        lines = block.strip().split("\n")
        mol_name = lines[0].strip() if lines else f"pose_{block_idx}"

        # Collect header fields from ## lines
        header_fields = {}
        for line in lines:
            m = _HEADER_PATTERN.match(line)
            if m:
                key, val = m.group(1), m.group(2).strip()
                try:
                    header_fields[key] = float(val)
                except ValueError:
                    header_fields[key] = val

        # Parse FPS score
        fps_score = None
        fps_vdw_score = None
        fps_es_score = None
        for line in lines:
            m = _FPS_EUCLIDEAN.match(line)
            if m:
                fps_score = float(m.group(1))
            m = _FPS_VDW_SCORE.match(line)
            if m:
                fps_vdw_score = float(m.group(1))
            m = _FPS_ES_SCORE.match(line)
            if m:
                fps_es_score = float(m.group(1))

        # Parse per-residue footprint
        ligand_vdw = {}  # {res_id: value}
        ligand_es = {}
        ref_vdw = {}
        ref_es = {}

        for line in lines:
            # Ligand vdW (try ligand-specific pattern first, excluding ref pattern)
            # We need to be careful: FPS.vdw is reference, FPS_vdw is ligand
            m = _FPS_LIGAND_VDW.match(line)
            if m and not line.strip().startswith("##########  FPS.vdw"):
                res_name, res_num, chain = m.group(1), m.group(2), m.group(3)
                res_id = f"{res_name}{res_num}.{chain}"
                ligand_vdw[res_id] = float(m.group(4))

            m = _FPS_LIGAND_ES.match(line)
            if m and not line.strip().startswith("##########  FPS.es"):
                res_name, res_num, chain = m.group(1), m.group(2), m.group(3)
                res_id = f"{res_name}{res_num}.{chain}"
                ligand_es[res_id] = float(m.group(4))

            # Reference vdW/ES
            m = _FPS_REF_VDW.match(line)
            if m:
                res_name, res_num, chain = m.group(1), m.group(2), m.group(3)
                res_id = f"{res_name}{res_num}.{chain}"
                ref_vdw[res_id] = float(m.group(4))

            m = _FPS_REF_ES.match(line)
            if m:
                res_name, res_num, chain = m.group(1), m.group(2), m.group(3)
                res_id = f"{res_name}{res_num}.{chain}"
                ref_es[res_id] = float(m.group(4))

        # Combine residue data
        all_residues = set(ligand_vdw) | set(ligand_es) | set(ref_vdw) | set(ref_es)
        residue_footprint = []

        for res_id in sorted(all_residues):
            # Parse residue_id back to components
            m_res = re.match(r"([A-Z]{1,4})(\d+)\.(\w+)", res_id)
            if not m_res:
                continue

            vdw_val = ligand_vdw.get(res_id, 0.0)
            es_val = ligand_es.get(res_id, 0.0)
            rvdw = ref_vdw.get(res_id, 0.0)
            res_val = ref_es.get(res_id, 0.0)

            residue_footprint.append({
                "residue_id": res_id,
                "residue_name": m_res.group(1),
                "residue_number": int(m_res.group(2)),
                "chain": m_res.group(3),
                "vdw": round(vdw_val, 4),
                "es": round(es_val, 4),
                "total": round(vdw_val + es_val, 4),
                "ref_vdw": round(rvdw, 4),
                "ref_es": round(res_val, 4),
                "ref_total": round(rvdw + res_val, 4),
                "delta_vdw": round(vdw_val - rvdw, 4),
                "delta_es": round(es_val - res_val, 4),
                "delta_total": round((vdw_val + es_val) - (rvdw + res_val), 4),
            })

        poses.append({
            "name": mol_name,
            "pose_id": block_idx,
            "header_fields": header_fields,
            "fps_score": fps_score,
            "fps_vdw_score": fps_vdw_score,
            "fps_es_score": fps_es_score,
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
        pharmacophore_threshold: float = 0.8,
        energy_cutoff: float = -0.5,
        best_pose_only: bool = True,
) -> Dict[str, Any]:
    """
    Analyze DOCK6 footprint re-scoring results.

    Reads footprint-scored mol2 files (from dock6_runner.run_footprint_rescoring)
    and produces per-residue energy analysis, consensus, and reference comparison.

    Args:
        footprint_dir:  Directory with {name}/{name}_footprint_scored.mol2
        output_dir:     Output directory for analysis files
        pharmacophore_threshold: Fraction of molecules that must contact a residue
                                 for it to be pharmacophoric (default: 0.8)
        energy_cutoff:  Minimum total energy for a residue to count as "contributing"
                        (kcal/mol, must be more negative than this)
        best_pose_only: If True, analyze only best pose per molecule (by FPS_Score)

    Returns:
        Dict with: success, n_molecules, output paths
    """
    footprint_dir = Path(footprint_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("  04b DOCK6 Footprint Analysis v1.0")
    logger.info("=" * 60)

    # --- Find footprint-scored mol2 files ---
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
        fps_mol2 = d / f"{name}_footprint_scored.mol2"
        if not fps_mol2.exists():
            continue

        poses = parse_footprint_mol2(str(fps_mol2))
        if not poses:
            logger.warning(f"  No poses parsed: {name}")
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
            fps_vdw = pose.get("fps_vdw_score")
            fps_es = pose.get("fps_es_score")

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
                "fps_vdw_score": fps_vdw,
                "fps_es_score": fps_es,
                "total_vdw": round(total_vdw, 3),
                "total_es": round(total_es, 3),
                "total_energy": round(total_vdw + total_es, 3),
                "n_residues_total": len(fp),
                "n_residues_contributing": n_contributing,
            })

        n_parsed += 1

    logger.info(f"  Parsed: {n_parsed} molecules, {n_failed} failed")

    if not all_rows:
        # No footprint data — may mean re-scoring hasn't been run yet
        # or DOCK6 didn't write footprint fields. Provide diagnostic info.
        logger.warning("  No footprint data found in any mol2 file!")
        logger.warning("  Possible causes:")
        logger.warning("    1. Footprint re-scoring hasn't been run (run 04b with --rescore)")
        logger.warning("    2. DOCK6 version doesn't write per-residue FPS fields")
        logger.warning("  Check: grep -i 'FPS_vdw\\|FPS_es' <mol2_file> | head -20")
        return {
            "success": False,
            "error": "No footprint data found. Run footprint re-scoring first.",
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
    # For each residue: how many molecules have it contributing?
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
            # Reference energies (same for all molecules)
            "ref_vdw": round(grp["ref_vdw"].iloc[0], 4) if "ref_vdw" in grp.columns else 0.0,
            "ref_es": round(grp["ref_es"].iloc[0], 4) if "ref_es" in grp.columns else 0.0,
        })

    df_consensus = pd.DataFrame(residue_stats)
    df_consensus.sort_values("mean_total", ascending=True, inplace=True)
    df_consensus.reset_index(drop=True, inplace=True)

    consensus_csv = output_dir / "residue_consensus.csv"
    df_consensus.to_csv(consensus_csv, index=False, encoding="utf-8")
    logger.info(f"  Saved: {consensus_csv} ({len(df_consensus)} residues)")

    # --- Pharmacophore residues (contacted by >threshold of molecules) ---
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
            "residues": pharma_list,
        }, f, indent=2)
    logger.info(f"  Saved: {pharma_json} ({len(pharma_list)} pharmacophore residues)")

    # --- vs reference comparison ---
    # Per-molecule delta from reference
    if "delta_total" in df_all.columns:
        # Pivot: for each molecule, which residues differ most from reference?
        ref_comparison = []
        for name, grp in df_all.groupby("Name"):
            sig = grp[grp["delta_total"].abs() > 0.5]  # significant deltas
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
            ref_csv = None
            logger.info("  No significant deltas vs reference found.")
    else:
        ref_csv = None

    # --- Save molecule summary ---
    summary_csv = output_dir / "molecule_footprint_summary.csv"
    df_summary.to_csv(summary_csv, index=False, encoding="utf-8")
    logger.info(f"  Saved: {summary_csv}")

    # --- Log summary ---
    logger.info("")
    logger.info(f"  Molecules analyzed:      {n_parsed}")
    logger.info(f"  Total residues tracked:  {len(df_consensus)}")
    logger.info(f"  Pharmacophore residues:  {len(pharma_list)} (>{pharmacophore_threshold*100:.0f}%)")
    if len(pharma_list) > 0:
        top5 = pharma_list[:5]
        for r in top5:
            logger.info(f"    {r['residue_id']}: {r['frac_contributing']*100:.0f}% "
                        f"(mean={r['mean_total']:.2f} kcal/mol)")
    logger.info("=" * 60)

    return {
        "success": True,
        "n_molecules": n_parsed,
        "n_residues": len(df_consensus),
        "n_pharmacophore": len(pharma_list),
        "footprint_per_molecule_csv": str(fps_csv),
        "residue_consensus_csv": str(consensus_csv),
        "pharmacophore_json": str(pharma_json),
        "vs_reference_csv": str(ref_csv) if ref_csv else None,
        "molecule_summary_csv": str(summary_csv),
        "output_dir": str(output_dir),
    }

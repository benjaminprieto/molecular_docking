#!/usr/bin/env python3
"""
GNINA Runner - Core Module (02b)
==================================
Executes GNINA docking for each prepared ligand.

GNINA = Vina + CNN scoring in a single run.
Produces three scores per pose:
  - vina_affinity: Empirical Vina score (kcal/mol, more negative = better)
  - cnn_score: Probability of correct pose (0-1, higher = better)
  - cnn_affinity: CNN-predicted affinity (kcal/mol)

Supports:
  - dock mode: full docking (find best pose)
  - rescore mode: score existing poses without moving them
  - autobox: define binding box from reference ligand

GNINA accepts PDB + SDF/mol2 directly — no PDBQT conversion needed.

Location: 01_src/molecular_docking/m02_gnina/gnina_runner.py

Project: molecular_docking
Module: 02b (core)
Version: 1.0
"""

import json
import logging
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class GninaResult:
    """Container for one molecule's GNINA docking result."""
    name: str
    success: bool
    vina_affinity: Optional[float] = None
    cnn_score: Optional[float] = None
    cnn_affinity: Optional[float] = None
    n_poses: int = 0
    output_file: Optional[str] = None
    ligand_file: Optional[str] = None
    ligand_format: Optional[str] = None
    error: Optional[str] = None
    runtime: float = 0.0
    binding_box_source: str = "json"


# =============================================================================
# GNINA DETECTION
# =============================================================================

def find_gnina(gnina_path: Optional[str] = None) -> Optional[str]:
    """Find GNINA executable."""
    paths_to_try = []
    if gnina_path:
        paths_to_try.append(gnina_path)

    paths_to_try.extend([
        "gnina",
        "/usr/local/bin/gnina",
        str(Path.home() / "bin" / "gnina"),
    ])

    for path in paths_to_try:
        try:
            result = subprocess.run(
                [path, "--version"], capture_output=True, text=True, timeout=10
            )
            if result.returncode == 0:
                return path
        except (subprocess.SubprocessError, FileNotFoundError):
            continue
    return None


def check_gnina_available(gnina_path: Optional[str] = None) -> Tuple[bool, str]:
    """Check if GNINA is available and return version info."""
    path = find_gnina(gnina_path)
    if path:
        try:
            result = subprocess.run(
                [path, "--version"], capture_output=True, text=True, timeout=10
            )
            version = (result.stdout.strip() or result.stderr.strip())[:120]
            return True, f"GNINA found at {path}: {version}"
        except Exception:
            return True, f"GNINA found at {path}"
    return False, "GNINA not found"


# =============================================================================
# GNINA OUTPUT PARSING
# =============================================================================

def parse_gnina_output(stdout: str) -> List[Dict[str, float]]:
    """
    Parse GNINA stdout to extract pose scores.

    GNINA output format:
      mode |   affinity | intramol | CNN-score | CNN-affinity
         1       -8.320    -0.123      0.847        -7.92
         2       -7.952    -0.098      0.721        -7.45
         ...

    Returns:
        List of dicts: [{mode, vina_affinity, cnn_score, cnn_affinity}, ...]
    """
    poses = []

    for line in stdout.split('\n'):
        stripped = line.strip()
        if not stripped or stripped.startswith('#') or 'mode' in stripped.lower() or '---' in stripped:
            continue

        parts = stripped.split()
        if len(parts) >= 4:
            try:
                mode = int(parts[0])
                vina_aff = float(parts[1])
                # GNINA output varies: sometimes 4 cols (mode, affinity, cnn_score, cnn_aff)
                # sometimes 5 cols (mode, affinity, intramol, cnn_score, cnn_aff)
                if len(parts) >= 5:
                    cnn_sc = float(parts[-2])
                    cnn_aff = float(parts[-1])
                else:
                    cnn_sc = float(parts[2])
                    cnn_aff = float(parts[3])

                poses.append({
                    "mode": mode,
                    "vina_affinity": vina_aff,
                    "cnn_score": cnn_sc,
                    "cnn_affinity": cnn_aff,
                })
            except (ValueError, IndexError):
                continue

    return poses


# =============================================================================
# LIGAND FILE SELECTION
# =============================================================================

def select_ligand_file(
        mol_input: Dict[str, Any],
        ligand_format: str = "auto",
) -> Tuple[Optional[str], str]:
    """
    Select ligand file based on format preference.

    GNINA can read both mol2 and SDF. Default: SDF from 00c.
    """
    sdf_file = mol_input.get('ligand_sdf')
    mol2_file = mol_input.get('ligand_mol2')

    sdf_exists = sdf_file and Path(sdf_file).exists()
    mol2_exists = mol2_file and Path(mol2_file).exists()

    if ligand_format == "mol2":
        return (mol2_file, "mol2") if mol2_exists else (None, "")
    elif ligand_format == "sdf":
        return (sdf_file, "sdf") if sdf_exists else (None, "")
    else:  # auto: prefer sdf (from 00c), fallback to mol2
        if sdf_exists:
            return sdf_file, "sdf"
        elif mol2_exists:
            return mol2_file, "mol2"
        return None, ""


# =============================================================================
# SINGLE MOLECULE DOCKING
# =============================================================================

def run_gnina_single(
        protein_pdb: str,
        ligand_file: str,
        output_file: str,
        # Binding box
        center_x: float = 0.0,
        center_y: float = 0.0,
        center_z: float = 0.0,
        size_x: float = 25.0,
        size_y: float = 25.0,
        size_z: float = 25.0,
        # Autobox
        autobox_ligand: Optional[str] = None,
        autobox_add: float = 6.0,
        use_autobox: bool = False,
        # GNINA params
        mode: str = "dock",
        cnn_scoring: str = "rescore",
        exhaustiveness: int = 8,
        num_modes: int = 9,
        seed: int = 42,
        # Execution
        gnina_path: str = "gnina",
        timeout: int = 600,
        atom_term_data: bool = True,
) -> GninaResult:
    """Run GNINA for a single ligand."""
    name = Path(ligand_file).stem
    fmt = Path(ligand_file).suffix.lstrip('.')
    start = time.time()

    # Validate inputs
    if not Path(protein_pdb).exists():
        return GninaResult(name=name, success=False,
                           error=f"Protein not found: {protein_pdb}")
    if not Path(ligand_file).exists():
        return GninaResult(name=name, success=False,
                           error=f"Ligand not found: {ligand_file}",
                           ligand_file=ligand_file, ligand_format=fmt)

    Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    # Resolve to absolute paths
    protein_pdb = str(Path(protein_pdb).resolve())
    ligand_file = str(Path(ligand_file).resolve())
    output_file = str(Path(output_file).resolve())

    # Build command
    cmd = [
        gnina_path,
        "--receptor", protein_pdb,
        "--ligand", ligand_file,
        "--out", output_file,
    ]

    # Binding box
    box_source = "json"
    if use_autobox and autobox_ligand and Path(autobox_ligand).exists():
        cmd.extend([
            "--autobox_ligand", str(Path(autobox_ligand).resolve()),
            "--autobox_add", str(autobox_add),
        ])
        box_source = "autobox"
    else:
        cmd.extend([
            "--center_x", f"{center_x:.3f}",
            "--center_y", f"{center_y:.3f}",
            "--center_z", f"{center_z:.3f}",
            "--size_x", f"{size_x:.1f}",
            "--size_y", f"{size_y:.1f}",
            "--size_z", f"{size_z:.1f}",
        ])

    # Mode
    if mode == "rescore":
        cmd.append("--score_only")
    else:
        cmd.extend([
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", str(num_modes),
        ])

    # CNN scoring
    cmd.extend(["--cnn_scoring", cnn_scoring])

    # Seed
    cmd.extend(["--seed", str(seed)])

    # Per-atom interaction terms embedded in output SDF
    if atom_term_data:
        cmd.append("--atom_term_data")

    # Execute
    try:
        logger.debug(f"  CMD: {' '.join(cmd)}")
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout,
        )
        elapsed = time.time() - start

        if result.returncode != 0:
            err = (result.stderr or result.stdout)[:300].strip()
            return GninaResult(
                name=name, success=False,
                error=f"GNINA exit {result.returncode}: {err}",
                ligand_file=ligand_file, ligand_format=fmt,
                runtime=elapsed, binding_box_source=box_source,
            )

        # Parse output
        poses = parse_gnina_output(result.stdout)
        if poses:
            best = poses[0]
            return GninaResult(
                name=name, success=True,
                vina_affinity=best["vina_affinity"],
                cnn_score=best["cnn_score"],
                cnn_affinity=best["cnn_affinity"],
                n_poses=len(poses),
                output_file=output_file,
                ligand_file=ligand_file, ligand_format=fmt,
                runtime=elapsed, binding_box_source=box_source,
            )
        else:
            return GninaResult(
                name=name, success=False,
                error="No poses parsed from GNINA output",
                ligand_file=ligand_file, ligand_format=fmt,
                runtime=elapsed, binding_box_source=box_source,
            )

    except subprocess.TimeoutExpired:
        return GninaResult(
            name=name, success=False, error=f"Timeout ({timeout}s)",
            ligand_file=ligand_file, ligand_format=fmt,
            runtime=time.time() - start, binding_box_source=box_source,
        )
    except Exception as e:
        return GninaResult(
            name=name, success=False, error=str(e),
            ligand_file=ligand_file, ligand_format=fmt,
            runtime=time.time() - start, binding_box_source=box_source,
        )


# =============================================================================
# BATCH DOCKING
# =============================================================================

def run_gnina_batch(
        inputs: List[Dict[str, Any]],
        output_dir: str,
        mode: str = "dock",
        cnn_scoring: str = "rescore",
        exhaustiveness: int = 8,
        num_modes: int = 9,
        seed: int = 42,
        n_workers: int = 2,
        timeout: int = 600,
        atom_term_data: bool = True,
        gnina_path: str = "gnina",
        ligand_format: str = "auto",
        # Autobox
        use_autobox: bool = False,
        autobox_ligand: Optional[str] = None,
        autobox_add: float = 6.0,
) -> List[GninaResult]:
    """Run GNINA on multiple molecules in parallel."""
    poses_dir = Path(output_dir) / "poses"
    poses_dir.mkdir(parents=True, exist_ok=True)

    results = []

    logger.info(f"Running GNINA {mode} on {len(inputs)} molecules "
                f"(workers={n_workers}, exh={exhaustiveness}, cnn={cnn_scoring})")

    def process_molecule(mol_input: Dict[str, Any]) -> GninaResult:
        name = mol_input['name']

        # Select ligand
        ligand_file, fmt = select_ligand_file(mol_input, ligand_format)
        if not ligand_file:
            return GninaResult(name=name, success=False,
                               error=f"No {ligand_format} file available")

        protein = mol_input.get('protein_pdb', '')
        output_file = str(poses_dir / f"{name}_gnina.sdf")

        # Binding box
        box = mol_input.get('binding_box', {}) or {}
        cx = box.get('center_x', 0.0)
        cy = box.get('center_y', 0.0)
        cz = box.get('center_z', 0.0)
        sx = box.get('size_x', 25.0)
        sy = box.get('size_y', 25.0)
        sz = box.get('size_z', 25.0)

        return run_gnina_single(
            protein_pdb=protein,
            ligand_file=ligand_file,
            output_file=output_file,
            center_x=cx, center_y=cy, center_z=cz,
            size_x=sx, size_y=sy, size_z=sz,
            autobox_ligand=autobox_ligand,
            autobox_add=autobox_add,
            use_autobox=use_autobox,
            mode=mode, cnn_scoring=cnn_scoring,
            exhaustiveness=exhaustiveness, num_modes=num_modes,
            seed=seed, gnina_path=gnina_path, timeout=timeout,
            atom_term_data=atom_term_data,
        )

    # Parallel execution
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(process_molecule, inp): inp for inp in inputs}

        for i, future in enumerate(as_completed(futures), 1):
            result = future.result()
            results.append(result)

            if result.success:
                score_str = (f"Vina={result.vina_affinity:.1f} "
                             f"CNN={result.cnn_score:.3f} "
                             f"CNN_aff={result.cnn_affinity:.1f}")
                logger.info(f"  [{i}/{len(inputs)}] ✓ {result.name} "
                            f"[{result.ligand_format}] {score_str} "
                            f"({result.n_poses} poses, {result.runtime:.1f}s)")
            else:
                err_short = result.error[:60] if result.error else "Unknown"
                logger.info(f"  [{i}/{len(inputs)}] ✗ {result.name}: {err_short}")

    return results


# =============================================================================
# LOAD INPUTS
# =============================================================================

def load_gnina_inputs_json(inputs_json: str) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """
    Load gnina_inputs.json from 02a.

    Returns:
        Tuple of (molecule_list, top_level_config)
    """
    with open(inputs_json, 'r') as f:
        data = json.load(f)

    if isinstance(data, list):
        return data, {}

    molecules = data.get('molecules', [])
    top_box = data.get('binding_box')
    if top_box:
        for mol in molecules:
            if not mol.get('binding_box'):
                mol['binding_box'] = top_box

    config = {
        "autobox_ligand": data.get("autobox_ligand"),
        "use_autobox": data.get("use_autobox", False),
    }

    return molecules, config


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def run_gnina_docking(
        inputs_json: str,
        output_dir: str,
        # GNINA params
        mode: str = "dock",
        cnn_scoring: str = "rescore",
        exhaustiveness: int = 8,
        num_modes: int = 9,
        seed: int = 42,
        # Execution
        n_workers: int = 2,
        timeout: int = 600,
        atom_term_data: bool = True,
        gnina_path: Optional[str] = None,
        ligand_format: str = "auto",
        # Autobox
        use_autobox: bool = False,
        autobox_ligand: Optional[str] = None,
        autobox_add: float = 6.0,
) -> Dict[str, Any]:
    """Run the complete GNINA docking pipeline."""
    logger.info("=" * 60)
    logger.info(f"GNINA DOCKING (02b) v1.0 — mode: {mode}")
    logger.info("=" * 60)

    # Resolve GNINA
    gnina = find_gnina(gnina_path)
    if not gnina:
        logger.error("GNINA not found!")
        return {"success": False, "error": "GNINA not found"}

    ok, version_msg = check_gnina_available(gnina)
    logger.info(f"  {version_msg}")

    # Load inputs
    logger.info(f"Loading inputs from {inputs_json}")
    inputs, top_config = load_gnina_inputs_json(inputs_json)
    logger.info(f"  Loaded {len(inputs)} molecules")

    # Resolve autobox from inputs JSON if not overridden
    if not use_autobox and top_config.get("use_autobox"):
        use_autobox = True
    if not autobox_ligand and top_config.get("autobox_ligand"):
        autobox_ligand = top_config["autobox_ligand"]

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Run docking
    start_time = time.time()
    results = run_gnina_batch(
        inputs=inputs,
        output_dir=output_dir,
        mode=mode, cnn_scoring=cnn_scoring,
        exhaustiveness=exhaustiveness, num_modes=num_modes,
        seed=seed, n_workers=n_workers, timeout=timeout,
        atom_term_data=atom_term_data,
        gnina_path=gnina, ligand_format=ligand_format,
        use_autobox=use_autobox, autobox_ligand=autobox_ligand,
        autobox_add=autobox_add,
    )
    total_time = time.time() - start_time

    # Convert to DataFrame
    results_data = []
    for r in results:
        results_data.append({
            "name": r.name,
            "success": r.success,
            "vina_affinity": r.vina_affinity,
            "cnn_score": r.cnn_score,
            "cnn_affinity": r.cnn_affinity,
            "n_poses": r.n_poses,
            "output_file": r.output_file,
            "ligand_file": r.ligand_file,
            "ligand_format": r.ligand_format,
            "error": r.error,
            "runtime": r.runtime,
            "binding_box_source": r.binding_box_source,
        })

    df = pd.DataFrame(results_data)

    # Merge with metadata
    inputs_df = pd.DataFrame(inputs)
    if 'metadata' in inputs_df.columns:
        meta_cols = pd.json_normalize(inputs_df['metadata'])
        meta_cols['name'] = inputs_df['name']
        df = df.merge(meta_cols, on='name', how='left')

    # Save CSV
    results_csv = output_path / "gnina_results.csv"
    df.to_csv(results_csv, index=False)
    logger.info(f"Saved: {results_csv}")

    # Save JSON
    results_json_file = output_path / "gnina_results.json"
    with open(results_json_file, 'w') as f:
        json.dump(results_data, f, indent=2, default=str)

    # Summary text
    summary_file = output_path / "gnina_summary.txt"
    _write_summary(summary_file, results, total_time, mode, cnn_scoring, exhaustiveness, num_modes)

    # Log summary
    n_success = int(df['success'].sum())
    logger.info("")
    logger.info("=" * 60)
    logger.info("GNINA DOCKING COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Mode:             {mode}")
    logger.info(f"CNN scoring:      {cnn_scoring}")
    logger.info(f"Total molecules:  {len(df)}")
    logger.info(f"Successful:       {n_success} ({100 * n_success / max(len(df), 1):.1f}%)")
    logger.info(f"Total time:       {total_time:.1f}s ({total_time / 60:.1f} min)")

    if n_success > 0:
        valid = df[df['success']]
        logger.info("")
        logger.info("Score statistics (best pose per molecule):")
        if valid['vina_affinity'].notna().any():
            logger.info(f"  Vina affinity:   {valid['vina_affinity'].min():.2f} to "
                        f"{valid['vina_affinity'].max():.2f} kcal/mol")
        if valid['cnn_score'].notna().any():
            logger.info(f"  CNN score:       {valid['cnn_score'].min():.3f} to "
                        f"{valid['cnn_score'].max():.3f}")
        if valid['cnn_affinity'].notna().any():
            logger.info(f"  CNN affinity:    {valid['cnn_affinity'].min():.2f} to "
                        f"{valid['cnn_affinity'].max():.2f} kcal/mol")

    return {
        "success": True,
        "n_molecules": len(df),
        "n_success": n_success,
        "total_time": total_time,
        "results_csv": str(results_csv),
        "results_json": str(results_json_file),
        "summary_txt": str(summary_file),
        "poses_dir": str(output_path / "poses"),
        "dataframe": df,
    }


# =============================================================================
# SUMMARY WRITER
# =============================================================================

def _write_summary(path, results, total_time, mode, cnn_scoring, exhaustiveness, num_modes):
    n_ok = sum(1 for r in results if r.success)
    n_fail = sum(1 for r in results if not r.success)
    w = 78

    lines = [
        "=" * w,
        "02b GNINA DOCKING — SUMMARY",
        "=" * w,
        "",
        f"Mode:              {mode}",
        f"CNN scoring:       {cnn_scoring}",
        f"Exhaustiveness:    {exhaustiveness}",
        f"Num modes:         {num_modes}",
        f"Total molecules:   {len(results)}",
        f"Successful:        {n_ok}",
        f"Failed:            {n_fail}",
        f"Total time:        {total_time:.0f}s ({total_time / 60:.1f} min)",
        f"Avg per molecule:  {total_time / max(n_ok, 1):.1f}s",
        "",
        "-" * w,
        f"{'Name':<30} {'Status':>7} {'Vina':>8} {'CNN_Sc':>7} {'CNN_Aff':>8} {'Poses':>5} {'Time':>7}",
        "-" * w,
    ]

    sorted_results = sorted(
        results,
        key=lambda r: r.cnn_score if r.cnn_score is not None else -999.0,
        reverse=True
    )

    for r in sorted_results:
        status = "OK" if r.success else "FAIL"
        vina = f"{r.vina_affinity:.2f}" if r.vina_affinity is not None else "—"
        cnn = f"{r.cnn_score:.3f}" if r.cnn_score is not None else "—"
        cnn_a = f"{r.cnn_affinity:.2f}" if r.cnn_affinity is not None else "—"
        poses = str(r.n_poses) if r.success else "—"
        lines.append(
            f"{r.name:<30} {status:>7} {vina:>8} {cnn:>7} {cnn_a:>8} {poses:>5} {r.runtime:>7.1f}"
        )

    if n_fail > 0:
        lines.extend(["", "FAILURES:"])
        for r in results:
            if not r.success:
                lines.append(f"  {r.name}: {r.error or 'unknown'}")

    lines.extend(["", "=" * w])

    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


if __name__ == '__main__':
    print("GNINA Runner - Core Module (02b) v1.0")
    print("Use 02b_gnina_runner.py CLI for execution")
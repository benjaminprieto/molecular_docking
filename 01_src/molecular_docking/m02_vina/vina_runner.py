#!/usr/bin/env python3
"""
Vina Runner - Core Module (02b)
=================================
Executes AutoDock Vina or Vina-GPU docking for each prepared ligand.

Supports:
  - AutoDock Vina 1.2.x (CPU, standard)
  - Vina-GPU 2.1 (GPU-accelerated)
  - Autobox with reference ligand
  - Parallel execution with ThreadPoolExecutor

Input:
  - vina_inputs.json from 02a (receptor PDBQT, ligand PDBQTs, binding box)

Output per molecule:
  - {name}/{name}_vina.pdbqt    Docked poses
  - {name}/vina.log             Vina stdout

Global output:
  - vina_results.csv            Results table
  - vina_results.json           Results in JSON format
  - vina_summary.txt            Human-readable summary

Location: 01_src/molecular_docking/m02_vina/vina_runner.py

Project: molecular_docking
Module: 02b (core)
Version: 1.0
"""

import json
import logging
import os
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
class VinaResult:
    """Container for one molecule's Vina docking result."""
    name: str
    success: bool
    affinity: Optional[float] = None        # Best pose (kcal/mol)
    rmsd_lb: Optional[float] = None         # RMSD lower bound (best pose)
    rmsd_ub: Optional[float] = None         # RMSD upper bound (best pose)
    n_poses: int = 0
    all_affinities: Optional[List[float]] = None  # All poses
    output_file: Optional[str] = None       # Output PDBQT
    ligand_file: Optional[str] = None       # Input PDBQT used
    log_file: Optional[str] = None
    error: Optional[str] = None
    runtime: float = 0.0
    engine: str = "vina"                    # "vina" or "vina-gpu"


# =============================================================================
# VINA DETECTION & VALIDATION
# =============================================================================

def find_vina(
        vina_path: Optional[str] = None,
        engine: str = "vina",
) -> Optional[str]:
    """
    Find Vina executable.

    Args:
        vina_path: Explicit path (highest priority).
        engine:    "vina" or "vina-gpu".

    Returns:
        Path to executable, or None.
    """
    paths_to_try = []

    if vina_path:
        paths_to_try.append(vina_path)

    if engine == "vina-gpu":
        paths_to_try.extend([
            "vina-gpu",
            "Vina-GPU",
            "/usr/local/bin/vina-gpu",
            str(Path.home() / "bin" / "vina-gpu"),
            "vina_gpu",
            str(Path.home() / "bin" / "Vina-GPU"),
        ])
    else:
        paths_to_try.extend([
            "vina",
            "/usr/local/bin/vina",
            str(Path.home() / "bin" / "vina"),
            "/opt/vina/bin/vina",
        ])

    for path in paths_to_try:
        try:
            result = subprocess.run(
                [path, "--version"],
                capture_output=True, text=True, timeout=10
            )
            if result.returncode == 0:
                return path
        except (subprocess.SubprocessError, FileNotFoundError):
            continue

    return None


def check_vina_available(
        vina_path: Optional[str] = None,
        engine: str = "vina",
) -> Tuple[bool, str]:
    """Check if Vina is available and return version info."""
    path = find_vina(vina_path, engine)
    if path:
        try:
            result = subprocess.run(
                [path, "--version"],
                capture_output=True, text=True, timeout=10
            )
            version = (result.stdout.strip() or result.stderr.strip())[:120]
            return True, f"{engine} found at {path}: {version}"
        except Exception:
            return True, f"{engine} found at {path}"
    return False, f"{engine} not found"


# =============================================================================
# VINA OUTPUT PARSING
# =============================================================================

def parse_vina_output(log_content: str) -> List[Dict[str, float]]:
    """
    Parse Vina stdout/log to extract pose scores.

    Vina output format:
      mode |   affinity | dist from best mode
           | (kcal/mol) | rmsd l.b.| rmsd u.b.
      -----+------------+----------+----------
         1       -8.320      0.000      0.000
         2       -7.952      1.834      2.456
         ...

    Returns:
        List of dicts: [{mode, affinity, rmsd_lb, rmsd_ub}, ...]
    """
    poses = []
    in_results = False

    for line in log_content.split('\n'):
        stripped = line.strip()

        # Detect results table
        if stripped.startswith("-----+"):
            in_results = True
            continue

        if in_results and stripped:
            parts = stripped.split()
            if len(parts) >= 4:
                try:
                    mode = int(parts[0])
                    affinity = float(parts[1])
                    rmsd_lb = float(parts[2])
                    rmsd_ub = float(parts[3])
                    poses.append({
                        "mode": mode,
                        "affinity": affinity,
                        "rmsd_lb": rmsd_lb,
                        "rmsd_ub": rmsd_ub,
                    })
                except (ValueError, IndexError):
                    # End of results table
                    if poses:
                        break

    return poses


# =============================================================================
# SINGLE MOLECULE DOCKING
# =============================================================================

def run_vina_single(
        receptor_pdbqt: str,
        ligand_pdbqt: str,
        output_pdbqt: str,
        log_file: str,
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
        # Search params
        exhaustiveness: int = 8,
        num_modes: int = 9,
        energy_range: float = 3.0,
        seed: int = 42,
        # Execution
        vina_path: str = "vina",
        timeout: int = 600,
        cpu: int = 1,
        # Vina-GPU specific
        engine: str = "vina",
        gpu_batch_size: int = 0,
) -> VinaResult:
    """
    Run Vina/Vina-GPU for a single ligand.

    Args:
        receptor_pdbqt: Receptor PDBQT file.
        ligand_pdbqt:   Ligand PDBQT file.
        output_pdbqt:   Output PDBQT with docked poses.
        log_file:       Log file path.
        center_x/y/z:   Binding box center (used if no autobox_ligand).
        size_x/y/z:     Binding box size (used if no autobox_ligand).
        autobox_ligand:  Reference ligand for autobox (overrides center/size).
        autobox_add:     Padding for autobox.
        exhaustiveness:  Search thoroughness.
        num_modes:       Max poses to generate.
        energy_range:    Max energy diff from best pose (kcal/mol).
        seed:            Random seed for reproducibility.
        vina_path:       Path to Vina executable.
        timeout:         Max seconds per molecule.
        cpu:             CPU threads per molecule.
        engine:          "vina" or "vina-gpu".
        gpu_batch_size:  Batch size for Vina-GPU (0 = auto).

    Returns:
        VinaResult with scores and file paths.
    """
    name = Path(ligand_pdbqt).stem.replace("_vina", "")
    start = time.time()

    # Validate inputs
    if not Path(receptor_pdbqt).exists():
        return VinaResult(
            name=name, success=False,
            error=f"Receptor not found: {receptor_pdbqt}",
            engine=engine,
        )
    if not Path(ligand_pdbqt).exists():
        return VinaResult(
            name=name, success=False,
            error=f"Ligand not found: {ligand_pdbqt}",
            ligand_file=ligand_pdbqt,
            engine=engine,
        )

    # Ensure output directory exists
    Path(output_pdbqt).parent.mkdir(parents=True, exist_ok=True)

    # Build command
    cmd = [
        vina_path,
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--out", str(output_pdbqt),
    ]

    # Binding box: autobox or explicit coordinates
    if autobox_ligand and Path(autobox_ligand).exists():
        cmd.extend([
            "--autobox_ligand", str(autobox_ligand),
            "--autobox_add", str(autobox_add),
        ])
    else:
        cmd.extend([
            "--center_x", f"{center_x:.3f}",
            "--center_y", f"{center_y:.3f}",
            "--center_z", f"{center_z:.3f}",
            "--size_x", f"{size_x:.1f}",
            "--size_y", f"{size_y:.1f}",
            "--size_z", f"{size_z:.1f}",
        ])

    # Search parameters
    cmd.extend([
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--energy_range", str(energy_range),
        "--seed", str(seed),
        "--cpu", str(cpu),
    ])

    # Vina-GPU specific
    if engine == "vina-gpu" and gpu_batch_size > 0:
        cmd.extend(["--batch_size", str(gpu_batch_size)])

    # Execute
    try:
        logger.debug(f"  CMD: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=str(Path(output_pdbqt).parent),
        )

        elapsed = time.time() - start

        # Save log
        log_content = result.stdout + "\n" + result.stderr
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        with open(log_file, 'w') as f:
            f.write(log_content)

        # Parse results
        if result.returncode == 0 and Path(output_pdbqt).exists():
            poses = parse_vina_output(log_content)
            if poses:
                best = poses[0]
                return VinaResult(
                    name=name,
                    success=True,
                    affinity=best["affinity"],
                    rmsd_lb=best["rmsd_lb"],
                    rmsd_ub=best["rmsd_ub"],
                    n_poses=len(poses),
                    all_affinities=[p["affinity"] for p in poses],
                    output_file=str(output_pdbqt),
                    ligand_file=ligand_pdbqt,
                    log_file=log_file,
                    runtime=elapsed,
                    engine=engine,
                )
            else:
                return VinaResult(
                    name=name, success=False,
                    error="No poses parsed from Vina output",
                    ligand_file=ligand_pdbqt,
                    log_file=log_file, runtime=elapsed,
                    engine=engine,
                )
        else:
            err = (result.stderr or result.stdout)[:300].strip()
            return VinaResult(
                name=name, success=False,
                error=f"Vina returned {result.returncode}: {err}",
                ligand_file=ligand_pdbqt,
                log_file=log_file, runtime=elapsed,
                engine=engine,
            )

    except subprocess.TimeoutExpired:
        elapsed = time.time() - start
        return VinaResult(
            name=name, success=False,
            error=f"Timeout ({timeout}s)",
            ligand_file=ligand_pdbqt,
            runtime=elapsed, engine=engine,
        )
    except Exception as e:
        elapsed = time.time() - start
        return VinaResult(
            name=name, success=False,
            error=str(e),
            ligand_file=ligand_pdbqt,
            runtime=elapsed, engine=engine,
        )


# =============================================================================
# BATCH DOCKING
# =============================================================================

def run_vina_batch(
        inputs: List[Dict[str, Any]],
        output_dir: str,
        # Search params
        exhaustiveness: int = 8,
        num_modes: int = 9,
        energy_range: float = 3.0,
        seed: int = 42,
        # Execution
        n_workers: int = 4,
        cpu_per_worker: int = 1,
        timeout: int = 600,
        vina_path: str = "vina",
        # Engine
        engine: str = "vina",
        gpu_batch_size: int = 0,
        # Autobox
        autobox_ligand: Optional[str] = None,
        autobox_add: float = 6.0,
        use_autobox: bool = False,
) -> List[VinaResult]:
    """
    Run Vina for a batch of molecules in parallel.

    Args:
        inputs:          List of molecule dicts from vina_inputs.json.
        output_dir:      Base output directory.
        exhaustiveness:  Search thoroughness.
        num_modes:       Max poses per molecule.
        energy_range:    Max energy difference from best.
        seed:            Random seed.
        n_workers:       Parallel workers.
        cpu_per_worker:  CPU threads per Vina process.
        timeout:         Seconds per molecule.
        vina_path:       Path to Vina executable.
        engine:          "vina" or "vina-gpu".
        gpu_batch_size:  Vina-GPU batch size.
        autobox_ligand:  Reference ligand for autobox.
        autobox_add:     Autobox padding.
        use_autobox:     If True, use autobox_ligand instead of coordinates.

    Returns:
        List of VinaResult objects.
    """
    output_path = Path(output_dir)
    poses_dir = output_path / "poses"
    poses_dir.mkdir(parents=True, exist_ok=True)

    results = []

    logger.info(f"Running {engine} on {len(inputs)} molecules "
                f"(workers={n_workers}, exh={exhaustiveness}, modes={num_modes})")

    def process_molecule(mol_input: Dict[str, Any]) -> VinaResult:
        name = mol_input['name']
        mol_out = poses_dir / name
        mol_out.mkdir(parents=True, exist_ok=True)

        out_pdbqt = mol_out / f"{name}_vina.pdbqt"
        log = mol_out / "vina.log"

        # Resolve binding box
        box = mol_input.get('binding_box', {})
        cx = box.get('center_x', 0.0)
        cy = box.get('center_y', 0.0)
        cz = box.get('center_z', 0.0)
        sx = box.get('size_x', 25.0)
        sy = box.get('size_y', 25.0)
        sz = box.get('size_z', 25.0)

        # Autobox override
        ab_ligand = None
        if use_autobox and autobox_ligand:
            ab_ligand = autobox_ligand

        return run_vina_single(
            receptor_pdbqt=mol_input['receptor_pdbqt'],
            ligand_pdbqt=mol_input['ligand_pdbqt'],
            output_pdbqt=str(out_pdbqt),
            log_file=str(log),
            center_x=cx, center_y=cy, center_z=cz,
            size_x=sx, size_y=sy, size_z=sz,
            autobox_ligand=ab_ligand,
            autobox_add=autobox_add,
            exhaustiveness=exhaustiveness,
            num_modes=num_modes,
            energy_range=energy_range,
            seed=seed,
            vina_path=vina_path,
            timeout=timeout,
            cpu=cpu_per_worker,
            engine=engine,
            gpu_batch_size=gpu_batch_size,
        )

    # Parallel execution
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(process_molecule, inp): inp for inp in inputs}

        for i, future in enumerate(as_completed(futures), 1):
            result = future.result()
            results.append(result)

            if result.success:
                score = f"ΔG={result.affinity:.2f}" if result.affinity else ""
                logger.info(f"  [{i}/{len(inputs)}] ✓ {result.name} "
                            f"{score} ({result.n_poses} poses, {result.runtime:.1f}s)")
            else:
                err_short = result.error[:60] if result.error else "Unknown"
                logger.info(f"  [{i}/{len(inputs)}] ✗ {result.name}: {err_short}")

    return results


# =============================================================================
# LOAD INPUTS
# =============================================================================

def load_vina_inputs_json(inputs_json: str) -> List[Dict[str, Any]]:
    """
    Load vina_inputs.json from 02a.

    Handles:
      - v1.0 dict format with 'molecules' key
      - Propagates top-level binding_box to molecules that lack one

    Returns:
        List of molecule dicts ready for run_vina_batch.
    """
    with open(inputs_json, 'r') as f:
        data = json.load(f)

    if isinstance(data, list):
        return data

    molecules = data.get('molecules', [])
    top_box = data.get('binding_box')
    if top_box:
        for mol in molecules:
            if not mol.get('binding_box'):
                mol['binding_box'] = top_box

    return molecules


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def run_vina_docking(
        inputs_json: str,
        output_dir: str,
        # Search
        exhaustiveness: int = 8,
        num_modes: int = 9,
        energy_range: float = 3.0,
        seed: int = 42,
        # Execution
        n_workers: int = 4,
        cpu_per_worker: int = 1,
        timeout: int = 600,
        vina_path: Optional[str] = None,
        # Engine
        engine: str = "vina",
        gpu_batch_size: int = 0,
        # Autobox
        use_autobox: bool = False,
        autobox_ligand: Optional[str] = None,
        autobox_add: float = 6.0,
) -> Dict[str, Any]:
    """
    Run the complete Vina docking pipeline.

    Args:
        inputs_json:     Path to vina_inputs.json from 02a.
        output_dir:      Output directory.
        exhaustiveness:  Search exhaustiveness.
        num_modes:       Max number of poses per molecule.
        energy_range:    Max energy difference (kcal/mol).
        seed:            Random seed.
        n_workers:       Parallel workers.
        cpu_per_worker:  CPU threads per Vina call.
        timeout:         Seconds per molecule.
        vina_path:       Path to Vina/Vina-GPU executable.
        engine:          "vina" or "vina-gpu".
        gpu_batch_size:  Vina-GPU batch size (0=auto).
        use_autobox:     Use autobox with reference ligand.
        autobox_ligand:  Reference ligand for autobox.
        autobox_add:     Autobox padding.

    Returns:
        Dict with results, stats, and file paths.
    """
    logger.info("=" * 60)
    logger.info(f"VINA DOCKING (02b) v1.0 — engine: {engine}")
    logger.info("=" * 60)

    # Resolve Vina executable
    vina = find_vina(vina_path, engine)
    if not vina:
        logger.error(f"{engine} not found!")
        return {"success": False, "error": f"{engine} not found"}

    ok, version_msg = check_vina_available(vina, engine)
    logger.info(f"  {version_msg}")

    # Load inputs
    logger.info(f"Loading inputs from {inputs_json}")
    inputs = load_vina_inputs_json(inputs_json)
    logger.info(f"  Loaded {len(inputs)} molecules")

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Run docking
    start_time = time.time()
    results = run_vina_batch(
        inputs=inputs,
        output_dir=output_dir,
        exhaustiveness=exhaustiveness,
        num_modes=num_modes,
        energy_range=energy_range,
        seed=seed,
        n_workers=n_workers,
        cpu_per_worker=cpu_per_worker,
        timeout=timeout,
        vina_path=vina,
        engine=engine,
        gpu_batch_size=gpu_batch_size,
        autobox_ligand=autobox_ligand,
        autobox_add=autobox_add,
        use_autobox=use_autobox,
    )
    total_time = time.time() - start_time

    # Convert to DataFrame
    results_data = []
    for r in results:
        results_data.append({
            "name": r.name,
            "success": r.success,
            "affinity": r.affinity,
            "rmsd_lb": r.rmsd_lb,
            "rmsd_ub": r.rmsd_ub,
            "n_poses": r.n_poses,
            "output_file": r.output_file,
            "ligand_file": r.ligand_file,
            "log_file": r.log_file,
            "error": r.error,
            "runtime": r.runtime,
            "engine": r.engine,
        })

    df = pd.DataFrame(results_data)

    # Merge with original metadata
    inputs_df = pd.DataFrame(inputs)
    if 'metadata' in inputs_df.columns:
        meta_cols = pd.json_normalize(inputs_df['metadata'])
        meta_cols['name'] = inputs_df['name']
        df = df.merge(meta_cols, on='name', how='left')

    # Save CSV
    results_csv = output_path / "vina_results.csv"
    df.to_csv(results_csv, index=False)
    logger.info(f"Saved: {results_csv}")

    # Save JSON
    results_json_file = output_path / "vina_results.json"
    with open(results_json_file, 'w') as f:
        json.dump(results_data, f, indent=2, default=str)

    # Summary text
    summary_file = output_path / "vina_summary.txt"
    _write_summary(summary_file, results, total_time, engine, exhaustiveness, num_modes)

    # Log summary
    n_success = int(df['success'].sum())
    logger.info("")
    logger.info("=" * 60)
    logger.info("VINA DOCKING COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Engine:           {engine}")
    logger.info(f"Total molecules:  {len(df)}")
    logger.info(f"Successful:       {n_success} ({100 * n_success / max(len(df), 1):.1f}%)")
    logger.info(f"Total time:       {total_time:.1f}s ({total_time / 60:.1f} min)")

    if n_success > 0:
        valid = df[df['success']]
        logger.info("")
        logger.info("Score statistics (best pose per molecule):")
        if valid['affinity'].notna().any():
            logger.info(f"  Affinity range:  {valid['affinity'].min():.2f} to "
                        f"{valid['affinity'].max():.2f} kcal/mol")
            logger.info(f"  Mean affinity:   {valid['affinity'].mean():.2f} kcal/mol")
        logger.info(f"  Avg poses:       {valid['n_poses'].mean():.1f}")
        logger.info(f"  Avg runtime:     {valid['runtime'].mean():.1f}s/mol")

    return {
        "success": True,
        "n_molecules": len(df),
        "n_success": n_success,
        "total_time": total_time,
        "engine": engine,
        "results_csv": str(results_csv),
        "results_json": str(results_json_file),
        "summary_txt": str(summary_file),
        "poses_dir": str(output_path / "poses"),
        "dataframe": df,
    }


# =============================================================================
# SUMMARY WRITER
# =============================================================================

def _write_summary(
        path: Path,
        results: List[VinaResult],
        total_time: float,
        engine: str,
        exhaustiveness: int,
        num_modes: int,
):
    """Write human-readable summary file."""
    n_ok = sum(1 for r in results if r.success)
    n_fail = sum(1 for r in results if not r.success)
    w = 70

    lines = [
        "=" * w,
        "02b VINA DOCKING — SUMMARY",
        "=" * w,
        "",
        f"Engine:            {engine}",
        f"Exhaustiveness:    {exhaustiveness}",
        f"Num modes:         {num_modes}",
        f"Total molecules:   {len(results)}",
        f"Successful:        {n_ok}",
        f"Failed:            {n_fail}",
        f"Total time:        {total_time:.0f}s ({total_time / 60:.1f} min)",
        f"Avg per molecule:  {total_time / max(n_ok, 1):.1f}s",
        "",
        "-" * w,
        f"{'Name':<30} {'Status':>8} {'Affinity':>10} {'Poses':>6} {'Time(s)':>8}",
        "-" * w,
    ]

    # Sort by affinity (best first)
    sorted_results = sorted(
        results,
        key=lambda r: r.affinity if r.affinity is not None else 999.0
    )

    for r in sorted_results:
        status = "OK" if r.success else "FAILED"
        aff = f"{r.affinity:.2f}" if r.affinity is not None else "—"
        poses = str(r.n_poses) if r.success else "—"
        lines.append(
            f"{r.name:<30} {status:>8} {aff:>10} {poses:>6} {r.runtime:>8.1f}"
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
    print("Vina Runner - Core Module (02b) v1.0")
    print("Use 02b_vina_runner.py CLI for execution")
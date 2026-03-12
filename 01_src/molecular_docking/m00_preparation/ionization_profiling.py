"""
Ionization Profiling - Core Module (00c)
==========================================
Wrapper around ionprofile to generate pH-dependent protonation states.

Reads unique_molecules.csv from 00a, calls ionprofile to:
  1. Protonate each molecule at each pH in the gradient
  2. Calculate formal charges at each pH
  3. Generate individual SDF files with 3D coords + explicit H

The user reviews outputs before proceeding to 00d (antechamber).

Location: 01_src/molecular_docking/m00_preparation/ionization_profiling.py
Project: molecular_docking
Module: 00c (core)
Version: 2.0
"""

import logging
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import pandas as pd

logger = logging.getLogger(__name__)


def _prepare_ionprofile_input(
        molecules_csv: Union[str, Path],
        work_dir: Path,
) -> Path:
    """
    Convert unique_molecules.csv into a format ionprofile can read.

    ionprofile expects: mol_id, smiles columns.
    00a outputs: Name, SMILES columns.

    Returns:
        Path to prepared CSV.
    """
    df = pd.read_csv(molecules_csv)

    # Normalize column names for ionprofile
    col_map = {}
    if "Name" in df.columns:
        col_map["Name"] = "mol_id"
    if "SMILES" in df.columns:
        col_map["SMILES"] = "smiles"

    # Try alternates
    if "mol_id" not in col_map.values():
        for col in ["name", "NAME", "ID", "Molecule_Name"]:
            if col in df.columns:
                col_map[col] = "mol_id"
                break
    if "smiles" not in col_map.values():
        for col in ["Smile", "smi", "canonical_smiles"]:
            if col in df.columns:
                col_map[col] = "smiles"
                break

    df = df.rename(columns=col_map)

    if "mol_id" not in df.columns or "smiles" not in df.columns:
        raise ValueError(
            f"Cannot find mol_id/smiles columns. "
            f"Available: {list(df.columns)}"
        )

    # Write minimal CSV for ionprofile
    work_dir.mkdir(parents=True, exist_ok=True)
    out_path = work_dir / "ionprofile_input.csv"
    df[["mol_id", "smiles"]].to_csv(out_path, index=False)

    logger.info(f"  Prepared {len(df)} molecules for ionprofile")
    return out_path


def run_ionization_profiling(
        molecules_csv: Union[str, Path],
        output_dir: Union[str, Path],
        docking_ph: Union[float, List[float]] = 7.2,
        ph_range: Optional[float] = None,
        ph_step: float = 0.1,
        engine: str = "openbabel",
        precision: float = 0.5,
        output_formats: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Run ionization profiling via ionprofile.

    Args:
        molecules_csv:  Path to unique_molecules.csv from 00a.
        output_dir:     Directory for outputs.
        docking_ph:     Target pH(s). Float or list of floats.
                        If list, generates SDF for each pH.
        ph_range:       If set, extends gradient ±ph_range around docking_ph.
                        E.g., docking_ph=7.2, ph_range=0.5 → 7.7-6.7
        ph_step:        Step between pH values (default 0.1).
        engine:         Ionprofile engine: "openbabel", "dimorphite", "qupkake".
        precision:      Engine precision parameter.
        output_formats: List of output formats. Default: ["csv", "excel", "sdf"].

    Returns:
        Dict with: n_molecules, ph_values, output_dir, output_files,
                   structures_dir, ionprofile_result
    """
    # --- Validate ionprofile is available ---
    try:
        from ionprofile import run_profiling
        from ionprofile.profiling.ionizer import check_dependencies
    except ImportError:
        return {
            "success": False,
            "error": "ionprofile not installed. "
                     "Install with: pip install git+https://github.com/benjaminprieto/ionization.git",
        }

    molecules_csv = Path(molecules_csv)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not molecules_csv.exists():
        return {
            "success": False,
            "error": f"Input CSV not found: {molecules_csv}. Run 00a first.",
        }

    if output_formats is None:
        output_formats = ["csv", "excel", "sdf"]

    # --- Determine pH gradient ---
    if isinstance(docking_ph, list):
        ph_max = max(docking_ph)
        ph_min = min(docking_ph)
        if ph_max == ph_min:
            ph_max += 0.1  # ionprofile needs a range
    else:
        if ph_range:
            ph_max = round(docking_ph + ph_range, 2)
            ph_min = round(docking_ph - ph_range, 2)
        else:
            # Single pH: create minimal gradient around it
            ph_max = docking_ph
            ph_min = docking_ph
            ph_step = 0.1

    # --- Check dependencies ---
    deps = check_dependencies(engine)
    if not deps.get("fully_operational"):
        logger.warning(f"Engine '{engine}' not fully operational: {deps}")

    # --- Prepare input ---
    work_dir = output_dir / "_work"
    input_csv = _prepare_ionprofile_input(molecules_csv, work_dir)

    # --- Run ionprofile ---
    logger.info("=" * 60)
    logger.info("  Ionization Profiling (via ionprofile)")
    logger.info("=" * 60)
    logger.info(f"  Input:      {molecules_csv.name}")
    logger.info(f"  Engine:     {engine}")
    logger.info(f"  pH:         {ph_min} -> {ph_max} (step {ph_step})")
    logger.info(f"  Formats:    {output_formats}")
    logger.info(f"  Output:     {output_dir}")

    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")

    try:
        result = run_profiling(
            input_path=str(input_csv),
            output_dir=str(output_dir),
            ph_max=ph_max,
            ph_min=ph_min,
            ph_step=ph_step,
            precision=precision,
            engine=engine,
            output_formats=output_formats,
            run_id=run_id,
            smiles_column="smiles",
            id_column="mol_id",
        )
    except Exception as e:
        logger.error(f"ionprofile failed: {e}")
        return {"success": False, "error": str(e)}

    # --- Locate SDF structures ---
    structures_dir = output_dir / run_id / "structures"
    ph_folders = {}
    if structures_dir.exists():
        for ph_dir in sorted(structures_dir.iterdir()):
            if ph_dir.is_dir() and ph_dir.name.startswith("pH"):
                n_sdf = len(list(ph_dir.glob("*.sdf")))
                if n_sdf > 0:
                    ph_folders[ph_dir.name] = {
                        "path": str(ph_dir),
                        "n_molecules": n_sdf,
                    }

    # --- Summary ---
    n_molecules = result.get("n_molecules", 0)

    logger.info("")
    logger.info("=" * 60)
    logger.info("  IONIZATION PROFILING COMPLETE")
    logger.info("=" * 60)
    logger.info(f"  Molecules:  {n_molecules}")
    logger.info(f"  pH values:  {len(result.get('ph_values', []))}")

    for fmt, fpath in result.get("output_files", {}).items():
        if isinstance(fpath, str):
            logger.info(f"  {fmt.upper():10s}: {fpath}")
        elif isinstance(fpath, dict):
            for k, v in fpath.items():
                logger.info(f"  {fmt.upper()} {k}: {v}")

    if ph_folders:
        logger.info("")
        logger.info("  SDF structures (for 00d antechamber):")
        for ph_label, info in ph_folders.items():
            logger.info(f"    {ph_label}/: {info['n_molecules']} molecules")
    else:
        logger.warning("  No SDF structures generated. "
                       "Check if 'sdf' is in output_formats.")

    logger.info("=" * 60)

    # --- Clean work dir ---
    shutil.rmtree(work_dir, ignore_errors=True)

    return {
        "success": True,
        "n_molecules": n_molecules,
        "ph_values": result.get("ph_values", []),
        "run_id": run_id,
        "output_dir": str(output_dir / run_id),
        "output_files": result.get("output_files", {}),
        "structures_dir": str(structures_dir) if structures_dir.exists() else None,
        "ph_folders": ph_folders,
    }

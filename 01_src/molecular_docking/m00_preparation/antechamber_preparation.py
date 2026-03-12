"""
Antechamber Preparation - Core Module (00d)
=============================================
Converts protonated SDF files into DOCK6-ready mol2 files.

Takes individual SDF files from 00c (ionization profiling) and runs
antechamber to assign AM1-BCC charges and GAFF2 atom types.

Pipeline per molecule:
    SDF (3D + protonated) → antechamber → mol2 (AM1-BCC, GAFF2)
    If antechamber fails → obabel fallback (Gasteiger charges)

Location: 01_src/molecular_docking/m00_preparation/antechamber_preparation.py
Project: molecular_docking
Module: 00d (core)
Version: 1.0
"""

import logging
import shutil
import subprocess
import tempfile
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import pandas as pd

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


# =============================================================================
# ANTECHAMBER
# =============================================================================

def run_antechamber(
        input_sdf: Union[str, Path],
        output_mol2: Union[str, Path],
        charge_method: str = "bcc",
        atom_type: str = "gaff2",
        net_charge: int = 0,
        timeout: int = 300,
) -> bool:
    """
    Run antechamber to convert SDF -> mol2 with partial charges.

    Args:
        input_sdf: Path to input SDF (must have 3D coords)
        output_mol2: Path for output mol2
        charge_method: "bcc" (AM1-BCC) or "gas" (Gasteiger)
        atom_type: "gaff2" or "gaff"
        net_charge: Net formal charge of the molecule
        timeout: Max seconds for antechamber

    Returns:
        True if successful
    """
    input_sdf = Path(input_sdf)
    output_mol2 = Path(output_mol2)
    output_mol2.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Copy input SDF to temp dir
        tmp_sdf = tmpdir / "input.sdf"
        shutil.copy2(input_sdf, tmp_sdf)

        tmp_mol2 = tmpdir / "output.mol2"

        cmd = [
            "antechamber",
            "-i", str(tmp_sdf),
            "-fi", "sdf",
            "-o", str(tmp_mol2),
            "-fo", "mol2",
            "-c", charge_method,
            "-at", atom_type,
            "-nc", str(net_charge),
            "-pf", "y",  # clean intermediate files
            "-dr", "n",  # no interactive prompts
        ]

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
                cwd=str(tmpdir),
            )

            if result.returncode != 0:
                logger.debug(f"    antechamber failed (rc={result.returncode})")
                if result.stderr:
                    logger.debug(f"    stderr: {result.stderr[:300]}")
                return False

            if not tmp_mol2.exists() or tmp_mol2.stat().st_size == 0:
                logger.debug("    antechamber produced no output")
                return False

            # Copy result to final location
            shutil.copy2(tmp_mol2, output_mol2)
            return True

        except FileNotFoundError:
            logger.error("    antechamber not found in PATH")
            return False
        except subprocess.TimeoutExpired:
            logger.warning(f"    antechamber timed out ({timeout}s)")
            return False
        except Exception as e:
            logger.debug(f"    antechamber error: {e}")
            return False


def mol2_from_obabel(
        input_sdf: Union[str, Path],
        output_mol2: Union[str, Path],
) -> bool:
    """
    Fallback: convert SDF -> mol2 using OpenBabel with Gasteiger charges.
    """
    try:
        output_mol2 = Path(output_mol2)
        output_mol2.parent.mkdir(parents=True, exist_ok=True)

        cmd = [
            "obabel", str(input_sdf),
            "-O", str(output_mol2),
            "--partialcharge", "gasteiger",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        return (
            result.returncode == 0
            and output_mol2.exists()
            and output_mol2.stat().st_size > 0
        )
    except Exception as e:
        logger.debug(f"    obabel fallback failed: {e}")
        return False


# =============================================================================
# MOL2 VALIDATION
# =============================================================================

def validate_mol2(mol2_path: Union[str, Path]) -> Dict[str, Any]:
    """
    Validate a mol2 file for DOCK6 compatibility.
    """
    result = {
        "exists": False,
        "has_atoms": False,
        "has_bonds": False,
        "has_charges": False,
        "n_atoms": 0,
        "n_bonds": 0,
        "valid": False,
    }

    path = Path(mol2_path)
    if not path.exists() or path.stat().st_size == 0:
        return result
    result["exists"] = True

    text = path.read_text()
    result["has_atoms"] = "@<TRIPOS>ATOM" in text
    result["has_bonds"] = "@<TRIPOS>BOND" in text

    if result["has_atoms"]:
        in_atom = False
        charges = []
        for line in text.split("\n"):
            if "@<TRIPOS>ATOM" in line:
                in_atom = True
                continue
            if line.startswith("@<TRIPOS>") and in_atom:
                break
            if in_atom and line.strip():
                parts = line.split()
                if len(parts) >= 9:
                    try:
                        charges.append(float(parts[8]))
                    except (ValueError, IndexError):
                        pass

        result["n_atoms"] = len(charges)
        result["has_charges"] = (
            len(charges) > 0
            and not all(abs(c) < 1e-10 for c in charges)
        )

    if result["has_bonds"]:
        in_bond = False
        n_bonds = 0
        for line in text.split("\n"):
            if "@<TRIPOS>BOND" in line:
                in_bond = True
                continue
            if line.startswith("@<TRIPOS>") and in_bond:
                break
            if in_bond and line.strip():
                n_bonds += 1
        result["n_bonds"] = n_bonds

    result["valid"] = (
        result["has_atoms"]
        and result["has_bonds"]
        and result["has_charges"]
        and result["n_atoms"] > 0
    )

    return result


# =============================================================================
# CHARGE DETECTION
# =============================================================================

def detect_formal_charge(sdf_path: Union[str, Path]) -> int:
    """
    Detect formal charge from SDF file.

    Tries:
        1. SDF property 'formal_charge'
        2. RDKit formal charge from structure
        3. Default 0
    """
    path = Path(sdf_path)

    if RDKIT_AVAILABLE:
        try:
            supplier = Chem.SDMolSupplier(str(path), removeHs=False)
            mol = next(iter(supplier), None)
            if mol is not None:
                # Check property first
                try:
                    return int(mol.GetProp("formal_charge"))
                except (KeyError, ValueError):
                    pass
                # Calculate from structure
                return Chem.GetFormalCharge(mol)
        except Exception:
            pass

    # Fallback: parse SDF text for formal_charge property
    try:
        text = path.read_text()
        for i, line in enumerate(text.split("\n")):
            if ">  <formal_charge>" in line or "> <formal_charge>" in line:
                lines = text.split("\n")
                if i + 1 < len(lines):
                    return int(lines[i + 1].strip())
    except Exception:
        pass

    return 0


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def run_antechamber_preparation(
        sdf_dir: Union[str, Path],
        output_dir: Union[str, Path],
        charge_method: str = "bcc",
        atom_type: str = "gaff2",
        obabel_fallback: bool = True,
        timeout_per_molecule: int = 300,
) -> Dict[str, Any]:
    """
    Run antechamber on all SDF files in a directory.

    Args:
        sdf_dir:              Directory with individual SDF files (from 00c).
        output_dir:           Directory for output mol2 files.
        charge_method:        "bcc" (AM1-BCC) or "gas" (Gasteiger).
        atom_type:            "gaff2" or "gaff".
        obabel_fallback:      Try obabel if antechamber fails.
        timeout_per_molecule: Max seconds per molecule for antechamber.

    Returns:
        Dict with: n_total, n_ok, n_failed, results, mol2_dir
    """
    sdf_dir = Path(sdf_dir)
    output_dir = Path(output_dir)
    mol2_dir = output_dir / "mol2"
    mol2_dir.mkdir(parents=True, exist_ok=True)

    if not sdf_dir.exists():
        return {
            "success": False,
            "error": f"SDF directory not found: {sdf_dir}. Run 00c first.",
        }

    sdf_files = sorted(sdf_dir.glob("*.sdf"))
    if not sdf_files:
        return {
            "success": False,
            "error": f"No SDF files found in {sdf_dir}",
        }

    logger.info("=" * 60)
    logger.info("  Antechamber Preparation")
    logger.info("=" * 60)
    logger.info(f"  SDF dir:    {sdf_dir} ({len(sdf_files)} files)")
    logger.info(f"  Output:     {mol2_dir}")
    logger.info(f"  Charges:    {charge_method} ({atom_type})")
    logger.info(f"  Timeout:    {timeout_per_molecule}s per molecule")
    logger.info(f"  Fallback:   {'obabel' if obabel_fallback else 'none'}")

    results = []
    n_total = len(sdf_files)

    for idx, sdf_path in enumerate(sdf_files):
        name = sdf_path.stem
        logger.info(f"  [{idx + 1}/{n_total}] {name}")

        mol2_path = mol2_dir / f"{name}.mol2"

        # Detect formal charge
        net_charge = detect_formal_charge(sdf_path)

        # Try antechamber
        success = run_antechamber(
            input_sdf=sdf_path,
            output_mol2=mol2_path,
            charge_method=charge_method,
            atom_type=atom_type,
            net_charge=net_charge,
            timeout=timeout_per_molecule,
        )
        method = f"antechamber_{charge_method}" if success else None

        # Fallback to obabel
        if not success and obabel_fallback:
            logger.info(f"    antechamber failed, trying obabel (Gasteiger)")
            success = mol2_from_obabel(sdf_path, mol2_path)
            if success:
                method = "obabel_gasteiger"

        # Validate
        if success:
            validation = validate_mol2(mol2_path)
            res = {
                "Name": name,
                "status": "OK" if validation["valid"] else "INVALID",
                "method": method,
                "net_charge": net_charge,
                "mol2_valid": validation["valid"],
                "mol2_n_atoms": validation["n_atoms"],
                "mol2_path": str(mol2_path),
                "sdf_path": str(sdf_path),
                "error": None if validation["valid"] else "mol2 validation failed",
            }
            if validation["valid"]:
                logger.info(f"    -> OK ({method}, "
                            f"{validation['n_atoms']} atoms, Q={net_charge})")
            else:
                logger.warning(f"    -> INVALID mol2 ({method})")
        else:
            res = {
                "Name": name,
                "status": "FAILED",
                "method": None,
                "net_charge": net_charge,
                "mol2_valid": False,
                "mol2_n_atoms": 0,
                "mol2_path": None,
                "sdf_path": str(sdf_path),
                "error": "antechamber + obabel both failed",
            }
            logger.warning(f"    -> FAILED (all methods)")

        results.append(res)

    # --- Summary ---
    n_ok = sum(1 for r in results if r["status"] == "OK")
    n_failed = n_total - n_ok
    methods = Counter(r["method"] for r in results if r["method"])

    logger.info("")
    logger.info("=" * 60)
    logger.info(f"  RESULTS: {n_ok}/{n_total} molecules prepared")
    if n_failed > 0:
        logger.info(f"  FAILED:  {n_failed} molecules")
    logger.info(f"  mol2:    {mol2_dir}/")
    logger.info("=" * 60)

    # --- Save summary CSV ---
    df_results = pd.DataFrame(results)
    csv_path = output_dir / "antechamber_summary.csv"
    df_results.to_csv(csv_path, index=False)

    # --- Save summary TXT ---
    w = 70
    lines = [
        "=" * w,
        "00d ANTECHAMBER PREPARATION - SUMMARY",
        "=" * w,
        "",
        f"Date:              {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"SDF input:         {sdf_dir}",
        f"Molecules:         {n_total}",
        f"Successful:        {n_ok}",
        f"Failed:            {n_failed}",
        f"Charge method:     {charge_method} ({atom_type})",
        "",
    ]

    if methods:
        lines.append("Methods used:")
        for m, count in sorted(methods.items()):
            lines.append(f"  {m}: {count}")
        lines.append("")

    lines.extend([
        "-" * w,
        f"{'Name':<30} {'Status':>8} {'Method':<25} {'Q':>4} {'Atoms':>6}",
        "-" * w,
    ])

    for r in results:
        method_str = r.get("method") or r.get("error", "—")
        if len(method_str) > 23:
            method_str = method_str[:20] + "..."
        lines.append(
            f"{r['Name']:<30} {r['status']:>8} {method_str:<25} "
            f"{r['net_charge']:>4} {r['mol2_n_atoms']:>6}"
        )

    if n_failed > 0:
        lines.extend(["", "FAILURES:"])
        for r in results:
            if r["status"] != "OK":
                lines.append(f"  {r['Name']}: {r.get('error', 'unknown')}")

    lines.extend(["", "=" * w])

    txt_path = output_dir / "antechamber_summary.txt"
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    return {
        "success": True,
        "n_total": n_total,
        "n_ok": n_ok,
        "n_failed": n_failed,
        "results": results,
        "mol2_dir": str(mol2_dir),
        "summary_csv": str(csv_path),
        "summary_txt": str(txt_path),
    }
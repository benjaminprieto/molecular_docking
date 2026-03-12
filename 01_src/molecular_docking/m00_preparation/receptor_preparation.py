"""
Receptor Preparation - Core Module (00b)
==========================================
Prepara el receptor para DOCK6: limpieza + protonacion al pH +
generacion de mol2 con cargas parciales.

DOCK6 requiere:
  - rec_charged.mol2  (Sybyl atom types + partial charges)
  - rec_noH.pdb       (sin H, para DMS surface en 01a)

Strategies (sin Chimera):
  A. pdb2pqr  → pH-aware via PROPKA, AMBER charges, inyectadas en mol2
  B. reduce   → AmberTools reduce + tleap ff14SB, cargas extraidas de prmtop
  C. obabel   → simple, Gasteiger charges (menos preciso pero funcional)

Todas generan mol2 con Sybyl atom types (via obabel) porque DOCK6 usa
vdw_AMBER_parm99.defn que mapea Sybyl types a parametros VdW.

Pipeline:
    1. Limpiar PDB (agua, alt conf, HETATM, cadenas)
    2. Protonar al pH de docking (herramienta segun strategy)
    3. Generar mol2 con Sybyl types (obabel)
    4. Inyectar cargas parciales (AMBER o Gasteiger)
    5. Generar rec_noH.pdb (strip H del protonado)
    6. Validar outputs
    7. Generar protonation_report.json

Location: 01_src/molecular_docking/m00_preparation/receptor_preparation.py
Project: molecular_docking
Module: 00b (core)
Version: 1.0
"""

import json
import logging
import re
import shutil
import subprocess
import tempfile
from collections import OrderedDict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple, Union

logger = logging.getLogger(__name__)


# =============================================================================
# PDB CLEANING
# =============================================================================

def clean_pdb(
        input_pdb: str,
        output_pdb: str,
        remove_water: bool = True,
        remove_hetatm: bool = True,
        remove_alt_conformations: bool = True,
        keep_chains: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Clean a PDB file for docking preparation.

    Operations:
      - Remove water molecules (HOH, WAT, TIP, TIP3)
      - Remove HETATM records (ligands, ions, cofactors)
      - Keep only first alternate conformation (altLoc A or ' ')
      - Filter by chain ID
      - Preserve CONECT records (for disulfide bonds)

    Args:
        input_pdb: Path to input PDB
        output_pdb: Path for cleaned PDB
        remove_water: Remove water molecules
        remove_hetatm: Remove HETATM records
        remove_alt_conformations: Keep only first alt conformation
        keep_chains: Chain IDs to keep (None = all)

    Returns:
        Dict with cleaning statistics
    """
    water_residues = {"HOH", "WAT", "TIP", "TIP3", "SOL"}
    stats = {
        "atoms_input": 0,
        "atoms_output": 0,
        "waters_removed": 0,
        "hetatm_removed": 0,
        "alt_conf_removed": 0,
        "chains_removed": 0,
    }

    output_lines = []

    with open(input_pdb) as f:
        for line in f:
            is_atom = line.startswith("ATOM")
            is_hetatm = line.startswith("HETATM")

            if is_atom or is_hetatm:
                stats["atoms_input"] += 1
                res_name = line[17:20].strip()
                chain_id = line[21].strip()
                alt_loc = line[16].strip()

                # Chain filter
                if keep_chains and chain_id not in keep_chains:
                    stats["chains_removed"] += 1
                    continue

                # Water filter
                if remove_water and res_name in water_residues:
                    stats["waters_removed"] += 1
                    continue

                # HETATM filter
                if remove_hetatm and is_hetatm:
                    stats["hetatm_removed"] += 1
                    continue

                # Alt conformation filter (keep '' or 'A')
                if remove_alt_conformations and alt_loc and alt_loc != "A":
                    stats["alt_conf_removed"] += 1
                    continue

                # Clear alt loc indicator for kept atoms
                if remove_alt_conformations and alt_loc == "A":
                    line = line[:16] + " " + line[17:]

                stats["atoms_output"] += 1
                output_lines.append(line)

            elif line.startswith("TER") or line.startswith("END"):
                output_lines.append(line)
            elif line.startswith("CONECT"):
                output_lines.append(line)
            elif line.startswith(("HEADER", "TITLE", "REMARK", "CRYST")):
                output_lines.append(line)

    # Write cleaned PDB
    Path(output_pdb).parent.mkdir(parents=True, exist_ok=True)
    with open(output_pdb, "w") as f:
        f.writelines(output_lines)

    logger.info(f"  Cleaned PDB: {stats['atoms_input']} -> {stats['atoms_output']} atoms")
    if stats["waters_removed"]:
        logger.info(f"    Removed {stats['waters_removed']} water atoms")
    if stats["hetatm_removed"]:
        logger.info(f"    Removed {stats['hetatm_removed']} HETATM atoms")
    if stats["alt_conf_removed"]:
        logger.info(f"    Removed {stats['alt_conf_removed']} alt conformations")

    return stats


def strip_hydrogens(input_pdb: str, output_pdb: str) -> int:
    """
    Remove hydrogen atoms from PDB and save as rec_noH.pdb.
    Needed for DMS surface generation in 01a.

    Returns:
        Number of hydrogen atoms removed
    """
    n_removed = 0
    output_lines = []

    with open(input_pdb) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16].strip()
                element = line[76:78].strip() if len(line) > 76 else ""

                # Detect hydrogen by element column or atom name
                is_hydrogen = (
                    element == "H"
                    or (not element and (atom_name.startswith("H") or atom_name in ("1H", "2H", "3H")))
                )
                if is_hydrogen:
                    n_removed += 1
                    continue

            output_lines.append(line)

    Path(output_pdb).parent.mkdir(parents=True, exist_ok=True)
    with open(output_pdb, "w") as f:
        f.writelines(output_lines)

    logger.info(f"  Stripped {n_removed} hydrogen atoms -> {Path(output_pdb).name}")
    return n_removed


# =============================================================================
# PQR CHARGE PARSING
# =============================================================================

def parse_pqr_charges(pqr_path: str) -> Dict[str, float]:
    """
    Parse a PQR file and extract per-atom partial charges.

    PQR format: like PDB but columns 55-62 = charge, 63-70 = radius

    Returns:
        Dict mapping "chainID:resNum:atomName" -> charge
    """
    charges = {}
    with open(pqr_path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            try:
                atom_name = line[12:16].strip()
                chain_id = line[21].strip() or "A"
                res_num = line[22:26].strip()
                # PQR: charge is after coordinates
                parts = line[30:].split()
                # x, y, z, charge, radius
                if len(parts) >= 5:
                    charge = float(parts[3])
                    key = f"{chain_id}:{res_num}:{atom_name}"
                    charges[key] = charge
            except (ValueError, IndexError):
                continue

    return charges


def inject_charges_into_mol2(
        mol2_path: str,
        charges: Dict[str, float],
        output_path: str,
) -> Tuple[int, int]:
    """
    Replace partial charges in a mol2 file with values from PQR.

    Matches atoms by residue number + atom name.

    Returns:
        (n_matched, n_total) atoms
    """
    lines = Path(mol2_path).read_text().split("\n")
    output_lines = []
    in_atom_section = False
    n_matched = 0
    n_total = 0

    for line in lines:
        if "@<TRIPOS>ATOM" in line:
            in_atom_section = True
            output_lines.append(line)
            continue
        if line.startswith("@<TRIPOS>") and in_atom_section:
            in_atom_section = False
            output_lines.append(line)
            continue

        if in_atom_section and line.strip():
            n_total += 1
            parts = line.split()
            if len(parts) >= 9:
                atom_name = parts[1]
                # parts[7] = subst_name like "ALA123" or "123" or "1"
                subst_name = parts[7] if len(parts) > 7 else ""
                # Extract residue number from subst_name
                res_match = re.search(r"(\d+)", subst_name)
                res_num = res_match.group(1) if res_match else "0"

                # Try matching with different chain IDs
                matched = False
                for chain in ["A", "B", "C", "D", ""]:
                    key = f"{chain}:{res_num}:{atom_name}"
                    if key in charges:
                        # Replace charge (last column before newline)
                        new_charge = f"{charges[key]:>10.4f}"
                        # Reconstruct line preserving alignment
                        # mol2 ATOM format: id name x y z type resid resname charge
                        parts[8] = f"{charges[key]:.4f}"
                        line = (
                            f"{parts[0]:>7s} {parts[1]:<8s} "
                            f"{float(parts[2]):>10.4f}{float(parts[3]):>10.4f}{float(parts[4]):>10.4f} "
                            f"{parts[5]:<8s} {parts[6]:>5s} {parts[7]:<8s} {charges[key]:>10.4f}"
                        )
                        n_matched += 1
                        matched = True
                        break

        output_lines.append(line)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write("\n".join(output_lines))

    return n_matched, n_total


# =============================================================================
# STRATEGY A: PDB2PQR (pH-aware, recommended)
# =============================================================================

def prepare_with_pdb2pqr(
        clean_pdb: str,
        output_dir: str,
        docking_ph: float = 7.2,
        force_field: str = "AMBER",
) -> Dict[str, Any]:
    """
    Protonate receptor using PDB2PQR + PROPKA and generate DOCK6 mol2.

    Pipeline:
      1. pdb2pqr --ff=AMBER --titration-state-method=propka --with-ph=X
      2. obabel: protonated PDB -> mol2 (Sybyl types)
      3. Inject AMBER charges from PQR into mol2

    This is the most accurate method for pH-specific docking because
    PROPKA predicts per-residue pKa values.
    """
    output_dir = Path(output_dir)
    pqr_path = output_dir / "receptor.pqr"
    protonated_pdb = output_dir / "receptor_protonated.pdb"

    # --- Step 1: PDB2PQR ---
    logger.info(f"  PDB2PQR: protonating at pH {docking_ph} (ff={force_field})")

    cmd = [
        "pdb2pqr",
        "--ff", force_field,
        "--titration-state-method", "propka",
        "--with-ph", str(docking_ph),
        "--keep-chain",
        "--pdb-output", str(protonated_pdb),
        "--log-level", "WARNING",
        str(clean_pdb),
        str(pqr_path),
    ]

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=300,
        )
        if result.returncode != 0:
            logger.error(f"    PDB2PQR failed (rc={result.returncode})")
            if result.stderr:
                logger.error(f"    {result.stderr[:500]}")
            return {"success": False, "error": "PDB2PQR failed", "tool": "pdb2pqr"}

        if not pqr_path.exists():
            return {"success": False, "error": "PDB2PQR produced no PQR output", "tool": "pdb2pqr"}

    except FileNotFoundError:
        return {"success": False, "error": "pdb2pqr not found in PATH", "tool": "pdb2pqr"}
    except subprocess.TimeoutExpired:
        return {"success": False, "error": "PDB2PQR timed out", "tool": "pdb2pqr"}

    logger.info(f"    PQR: {pqr_path.name}")
    if protonated_pdb.exists():
        logger.info(f"    Protonated PDB: {protonated_pdb.name}")

    # --- Step 2: Generate mol2 with Sybyl types via obabel ---
    # Use the protonated PDB if available, else PQR
    source_for_mol2 = str(protonated_pdb) if protonated_pdb.exists() else str(clean_pdb)
    mol2_gasteiger = output_dir / "rec_sybyl_gasteiger.mol2"

    logger.info("  obabel: PDB -> mol2 (Sybyl atom types)")
    obabel_cmd = [
        "obabel", source_for_mol2,
        "-O", str(mol2_gasteiger),
        "--partialcharge", "gasteiger",
    ]
    try:
        ob_result = subprocess.run(obabel_cmd, capture_output=True, text=True, timeout=120)
        if ob_result.returncode != 0 or not mol2_gasteiger.exists():
            return {"success": False, "error": "obabel mol2 conversion failed", "tool": "pdb2pqr"}
    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        return {"success": False, "error": f"obabel failed: {e}", "tool": "pdb2pqr"}

    # --- Step 3: Inject AMBER charges from PQR ---
    logger.info("  Injecting AMBER charges from PQR into mol2...")
    charges = parse_pqr_charges(str(pqr_path))

    rec_charged = output_dir / "rec_charged.mol2"
    n_matched, n_total = inject_charges_into_mol2(
        str(mol2_gasteiger), charges, str(rec_charged),
    )

    match_pct = (n_matched / n_total * 100) if n_total > 0 else 0
    logger.info(f"    Charges injected: {n_matched}/{n_total} atoms ({match_pct:.1f}%)")

    if match_pct < 50:
        logger.warning(f"    Low charge match rate ({match_pct:.1f}%). "
                       "Residue numbering may differ between PQR and mol2.")

    # --- Parse PROPKA results ---
    propka_report = _parse_propka_log(output_dir)

    return {
        "success": True,
        "tool": "pdb2pqr",
        "rec_charged_mol2": str(rec_charged),
        "pqr_path": str(pqr_path),
        "protonated_pdb": str(protonated_pdb) if protonated_pdb.exists() else None,
        "charge_match_rate": round(match_pct, 1),
        "n_pqr_charges": len(charges),
        "propka_report": propka_report,
    }


def _parse_propka_log(output_dir: Path) -> Dict[str, Any]:
    """Try to parse PROPKA pKa predictions from PDB2PQR output."""
    report = {"titratable_residues": []}

    # PDB2PQR may write propka output to various locations
    propka_files = list(output_dir.glob("*.propka")) + list(output_dir.glob("*.pka"))

    for propka_file in propka_files:
        try:
            text = propka_file.read_text()
            for line in text.split("\n"):
                # Look for lines like: "ASP  48 A   3.45"
                match = re.match(r"\s*(ASP|GLU|HIS|LYS|CYS|TYR)\s+(\d+)\s+(\S+)\s+([\d.]+)", line)
                if match:
                    report["titratable_residues"].append({
                        "residue": match.group(1),
                        "number": int(match.group(2)),
                        "chain": match.group(3),
                        "pKa": float(match.group(4)),
                    })
        except Exception:
            pass

    return report


# =============================================================================
# STRATEGY B: reduce + tleap
# =============================================================================

def prepare_with_reduce(
        clean_pdb: str,
        output_dir: str,
        docking_ph: float = 7.2,
) -> Dict[str, Any]:
    """
    Protonate receptor using AmberTools reduce + tleap.

    Pipeline:
      1. reduce -build: add hydrogens
      2. tleap: load with ff14SB, savemol2
      3. obabel: convert tleap mol2 (AMBER types) -> mol2 (Sybyl types)
         preserving coordinates and charges

    Note: reduce is less pH-aware than PDB2PQR/PROPKA.
    HIS protonation is determined by local environment, not pH.
    """
    output_dir = Path(output_dir)
    reduced_pdb = output_dir / "receptor_reduced.pdb"

    # --- Step 1: reduce ---
    logger.info("  reduce: adding hydrogens")
    try:
        result = subprocess.run(
            ["reduce", "-build", "-nuclear", str(clean_pdb)],
            capture_output=True, text=True, timeout=120,
        )
        # reduce writes to stdout
        if result.stdout:
            reduced_pdb.write_text(result.stdout)
            logger.info(f"    -> {reduced_pdb.name}")
        else:
            return {"success": False, "error": "reduce produced no output", "tool": "reduce"}
    except FileNotFoundError:
        return {"success": False, "error": "reduce not found in PATH", "tool": "reduce"}
    except subprocess.TimeoutExpired:
        return {"success": False, "error": "reduce timed out", "tool": "reduce"}

    # --- Step 2: tleap -> prmtop (for charges) + mol2 ---
    logger.info("  tleap: assigning ff14SB charges")
    tleap_mol2 = output_dir / "rec_tleap.mol2"
    prmtop = output_dir / "receptor.prmtop"
    inpcrd = output_dir / "receptor.inpcrd"

    tleap_script = f"""\
source leaprc.protein.ff14SB
rec = loadpdb {reduced_pdb}
savemol2 rec {tleap_mol2} 1
saveamberparm rec {prmtop} {inpcrd}
quit
"""
    tleap_in = output_dir / "tleap.in"
    tleap_in.write_text(tleap_script)

    try:
        result = subprocess.run(
            ["tleap", "-f", str(tleap_in)],
            capture_output=True, text=True, timeout=120,
            cwd=str(output_dir),
        )
        # tleap may show warnings but still succeed
        if not tleap_mol2.exists() or tleap_mol2.stat().st_size == 0:
            logger.warning("    tleap savemol2 failed, trying obabel fallback")
    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        logger.warning(f"    tleap failed: {e}")

    # --- Step 3: Generate Sybyl-typed mol2 via obabel ---
    # tleap's savemol2 outputs AMBER atom types, not Sybyl.
    # We use obabel on the reduced PDB to get proper Sybyl types.
    rec_charged = output_dir / "rec_charged.mol2"

    logger.info("  obabel: generating mol2 with Sybyl types")
    obabel_cmd = [
        "obabel", str(reduced_pdb),
        "-O", str(rec_charged),
        "--partialcharge", "gasteiger",
    ]

    try:
        ob_result = subprocess.run(obabel_cmd, capture_output=True, text=True, timeout=120)
        if ob_result.returncode != 0 or not rec_charged.exists():
            return {"success": False, "error": "obabel conversion failed", "tool": "reduce"}
    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        return {"success": False, "error": f"obabel failed: {e}", "tool": "reduce"}

    # If tleap produced charges, inject them
    if tleap_mol2.exists() and tleap_mol2.stat().st_size > 0:
        logger.info("  Extracting AMBER charges from tleap mol2...")
        tleap_charges = _extract_charges_from_mol2(str(tleap_mol2))
        if tleap_charges:
            n_matched, n_total = inject_charges_into_mol2(
                str(rec_charged), tleap_charges, str(rec_charged),
            )
            match_pct = (n_matched / n_total * 100) if n_total > 0 else 0
            logger.info(f"    AMBER charges injected: {n_matched}/{n_total} ({match_pct:.1f}%)")

    return {
        "success": True,
        "tool": "reduce",
        "rec_charged_mol2": str(rec_charged),
        "reduced_pdb": str(reduced_pdb),
        "tleap_mol2": str(tleap_mol2) if tleap_mol2.exists() else None,
        "prmtop": str(prmtop) if prmtop.exists() else None,
    }


def _extract_charges_from_mol2(mol2_path: str) -> Dict[str, float]:
    """Extract charges from a mol2 file, keyed by residue:atom_name."""
    charges = {}
    in_atom = False

    for line in Path(mol2_path).read_text().split("\n"):
        if "@<TRIPOS>ATOM" in line:
            in_atom = True
            continue
        if line.startswith("@<TRIPOS>") and in_atom:
            break
        if in_atom and line.strip():
            parts = line.split()
            if len(parts) >= 9:
                atom_name = parts[1]
                subst_name = parts[7] if len(parts) > 7 else ""
                res_match = re.search(r"(\d+)", subst_name)
                res_num = res_match.group(1) if res_match else "0"
                try:
                    charge = float(parts[8])
                    for chain in ["A", "B", "C", "D", ""]:
                        charges[f"{chain}:{res_num}:{atom_name}"] = charge
                except (ValueError, IndexError):
                    pass
    return charges


# =============================================================================
# STRATEGY C: obabel (simple)
# =============================================================================

def prepare_with_obabel(
        clean_pdb: str,
        output_dir: str,
        docking_ph: float = 7.2,
) -> Dict[str, Any]:
    """
    Protonate receptor using OpenBabel.

    Simple approach: obabel -p pH adds hydrogens, then generates mol2
    with Gasteiger charges and Sybyl atom types.

    Less accurate than PDB2PQR for pH-specific protonation,
    but always works and produces valid DOCK6 input.
    """
    output_dir = Path(output_dir)
    rec_charged = output_dir / "rec_charged.mol2"

    logger.info(f"  obabel: protonating at pH {docking_ph} + Gasteiger charges")

    cmd = [
        "obabel", str(clean_pdb),
        "-O", str(rec_charged),
        "-p", str(docking_ph),
        "--partialcharge", "gasteiger",
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0 or not rec_charged.exists() or rec_charged.stat().st_size == 0:
            return {"success": False, "error": "obabel failed", "tool": "obabel"}
    except FileNotFoundError:
        return {"success": False, "error": "obabel not found in PATH", "tool": "obabel"}
    except subprocess.TimeoutExpired:
        return {"success": False, "error": "obabel timed out", "tool": "obabel"}

    logger.info(f"    -> {rec_charged.name} ({rec_charged.stat().st_size} bytes)")

    return {
        "success": True,
        "tool": "obabel",
        "rec_charged_mol2": str(rec_charged),
    }


# =============================================================================
# MOL2 VALIDATION
# =============================================================================

def validate_prepared_mol2(mol2_path: Union[str, Path]) -> Dict[str, Any]:
    """
    Validate a receptor mol2 for DOCK6 compatibility.

    Checks:
      - File exists and non-empty
      - Has @<TRIPOS>ATOM and @<TRIPOS>BOND sections
      - Reasonable atom count (>100 for a protein)
      - Charges are present (not all zero)
      - Sybyl atom types present (C.3, N.am, O.2, etc.)
    """
    result = {
        "exists": False,
        "has_atoms": False,
        "has_bonds": False,
        "has_charges": False,
        "has_sybyl_types": False,
        "n_atoms": 0,
        "n_residues": 0,
        "total_charge": 0.0,
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
        atom_types = set()
        residue_ids = set()

        for line in text.split("\n"):
            if "@<TRIPOS>ATOM" in line:
                in_atom = True
                continue
            if line.startswith("@<TRIPOS>") and in_atom:
                break
            if in_atom and line.strip():
                parts = line.split()
                if len(parts) >= 9:
                    atom_types.add(parts[5])
                    try:
                        charges.append(float(parts[8]))
                    except (ValueError, IndexError):
                        pass
                    if len(parts) > 6:
                        residue_ids.add(parts[6])

        result["n_atoms"] = len(charges)
        result["n_residues"] = len(residue_ids)
        result["has_charges"] = (
            len(charges) > 0
            and not all(abs(c) < 1e-10 for c in charges)
        )
        result["total_charge"] = round(sum(charges), 2) if charges else 0.0

        # Check for Sybyl types (contain dots: C.3, N.am, O.2, S.3, etc.)
        sybyl_types = [t for t in atom_types if "." in t]
        result["has_sybyl_types"] = len(sybyl_types) > 0

    result["valid"] = (
        result["has_atoms"]
        and result["has_bonds"]
        and result["has_charges"]
        and result["n_atoms"] > 50  # A protein should have many atoms
    )

    return result


# =============================================================================
# MAIN PIPELINE FUNCTION
# =============================================================================

def run_receptor_preparation(
        receptor_pdb: Union[str, Path],
        output_dir: Union[str, Path],
        docking_ph: float = 7.2,
        protonation_tool: str = "pdb2pqr",
        force_field: str = "AMBER",
        chain: Optional[str] = None,
        remove_water: bool = True,
        remove_hetatm: bool = True,
        remove_alt_conformations: bool = True,
) -> Dict[str, Any]:
    """
    Run the complete receptor preparation pipeline.

    Args:
        receptor_pdb: Path to input PDB file
        output_dir: Directory for output files
        docking_ph: pH for protonation
        protonation_tool: "pdb2pqr" | "reduce" | "obabel"
        force_field: "AMBER" | "CHARMM" | "PARSE" (for PDB2PQR)
        chain: Chain ID to keep (None = all)
        remove_water: Remove water molecules
        remove_hetatm: Remove HETATM records
        remove_alt_conformations: Keep only first alt conformation

    Returns:
        Dict with output paths, validation, and protonation report
    """
    receptor_pdb = Path(receptor_pdb)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not receptor_pdb.exists():
        return {"success": False, "error": f"Receptor PDB not found: {receptor_pdb}"}

    logger.info("=" * 60)
    logger.info("  Receptor Preparation Pipeline")
    logger.info("=" * 60)
    logger.info(f"  Input:   {receptor_pdb.name}")
    logger.info(f"  pH:      {docking_ph}")
    logger.info(f"  Tool:    {protonation_tool}")
    logger.info(f"  Output:  {output_dir}")

    report = {
        "start_time": datetime.now().isoformat(),
        "input_pdb": str(receptor_pdb),
        "docking_ph": docking_ph,
        "protonation_tool": protonation_tool,
    }

    # --- Step 1: Clean PDB ---
    logger.info("\nStep 1: Cleaning PDB")
    clean_path = output_dir / "receptor_clean.pdb"
    keep_chains = [chain] if chain else None

    clean_stats = clean_pdb(
        str(receptor_pdb), str(clean_path),
        remove_water=remove_water,
        remove_hetatm=remove_hetatm,
        remove_alt_conformations=remove_alt_conformations,
        keep_chains=keep_chains,
    )
    report["clean_stats"] = clean_stats

    # --- Step 2: Protonate + generate mol2 ---
    logger.info(f"\nStep 2: Protonation ({protonation_tool})")

    strategies = {
        "pdb2pqr": lambda: prepare_with_pdb2pqr(
            str(clean_path), str(output_dir), docking_ph, force_field,
        ),
        "reduce": lambda: prepare_with_reduce(
            str(clean_path), str(output_dir), docking_ph,
        ),
        "obabel": lambda: prepare_with_obabel(
            str(clean_path), str(output_dir), docking_ph,
        ),
    }

    if protonation_tool not in strategies:
        return {"success": False,
                "error": f"Unknown protonation tool: {protonation_tool}. "
                         f"Options: {list(strategies.keys())}"}

    prep_result = strategies[protonation_tool]()

    if not prep_result.get("success"):
        # Try fallback to obabel
        if protonation_tool != "obabel":
            logger.warning(f"  {protonation_tool} failed, trying obabel fallback...")
            prep_result = prepare_with_obabel(
                str(clean_path), str(output_dir), docking_ph,
            )

    if not prep_result.get("success"):
        report["error"] = prep_result.get("error", "All protonation methods failed")
        return {"success": False, "report": report, **prep_result}

    report["protonation_result"] = {
        k: v for k, v in prep_result.items() if k != "success"
    }

    rec_charged = prep_result["rec_charged_mol2"]

    # --- Step 3: Generate rec_noH.pdb ---
    logger.info("\nStep 3: Generating rec_noH.pdb")
    rec_noH = output_dir / "rec_noH.pdb"

    # Source: protonated PDB if available, else clean PDB
    protonated_pdb = prep_result.get("protonated_pdb") or prep_result.get("reduced_pdb")
    source_pdb = protonated_pdb if protonated_pdb and Path(protonated_pdb).exists() else str(clean_path)
    strip_hydrogens(source_pdb, str(rec_noH))

    # --- Step 4: Validate mol2 ---
    logger.info("\nStep 4: Validating rec_charged.mol2")
    validation = validate_prepared_mol2(rec_charged)
    report["validation"] = validation

    if validation["valid"]:
        logger.info(f"  VALID: {validation['n_atoms']} atoms, "
                     f"{validation['n_residues']} residues, "
                     f"total charge: {validation['total_charge']}")
        if validation["has_sybyl_types"]:
            logger.info("  Sybyl atom types: present")
        else:
            logger.warning("  WARNING: Sybyl atom types NOT detected "
                          "(DOCK6 may not score correctly)")
    else:
        logger.warning(f"  WARNING: mol2 validation issues: {validation}")

    # --- Save report ---
    report["end_time"] = datetime.now().isoformat()
    report["success"] = True

    report_path = output_dir / "protonation_report.json"
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2, default=str)

    # --- Summary TXT ---
    summary_path = output_dir / "protonation_summary.txt"
    w = 70
    lines = [
        "=" * w,
        "00b RECEPTOR PREPARATION - SUMMARY",
        "=" * w,
        "",
        f"Date:              {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"Input:             {receptor_pdb.name}",
        f"pH:                {docking_ph}",
        f"Tool:              {prep_result.get('tool', protonation_tool)}",
        "",
        f"Atoms (input):     {clean_stats['atoms_input']}",
        f"Atoms (cleaned):   {clean_stats['atoms_output']}",
        f"Waters removed:    {clean_stats['waters_removed']}",
        f"HETATM removed:    {clean_stats['hetatm_removed']}",
        "",
        f"mol2 atoms:        {validation.get('n_atoms', 0)}",
        f"mol2 residues:     {validation.get('n_residues', 0)}",
        f"Total charge:      {validation.get('total_charge', 0):.2f}",
        f"Sybyl types:       {'yes' if validation.get('has_sybyl_types') else 'NO'}",
        f"Valid:             {'YES' if validation.get('valid') else 'NO'}",
        "",
    ]

    # PROPKA results if available
    propka = prep_result.get("propka_report", {})
    titratable = propka.get("titratable_residues", [])
    if titratable:
        lines.extend([
            "-" * w,
            "PROPKA pKa Predictions (titratable residues)",
            "-" * w,
        ])
        # Show residues where protonation state differs from standard
        for res in titratable:
            marker = ""
            if res["residue"] == "HIS" and res["pKa"] > docking_ph:
                marker = " <- PROTONATED (HIP)"
            elif res["residue"] in ("ASP", "GLU") and res["pKa"] > docking_ph:
                marker = " <- NEUTRAL"
            elif res["residue"] == "LYS" and res["pKa"] < docking_ph:
                marker = " <- NEUTRAL"
            elif res["residue"] == "CYS" and res["pKa"] < docking_ph:
                marker = " <- DEPROTONATED"

            lines.append(
                f"  {res['residue']}{res['number']:>4d} {res['chain']}: "
                f"pKa = {res['pKa']:.2f}{marker}"
            )
        lines.append("")

    lines.extend([
        "=" * w,
        "",
        f"rec_charged.mol2:  {rec_charged}",
        f"rec_noH.pdb:       {rec_noH}",
    ])

    summary_path.write_text("\n".join(lines))
    logger.info(f"  Summary: {summary_path}")

    logger.info("")
    logger.info(f"{'=' * 60}")
    logger.info(f"  Receptor prepared: {Path(rec_charged).name}")
    logger.info(f"  rec_noH.pdb:       {rec_noH.name}")
    logger.info(f"{'=' * 60}")

    return {
        "success": True,
        "rec_charged_mol2": rec_charged,
        "rec_noH_pdb": str(rec_noH),
        "report": report,
        "report_path": str(report_path),
        "summary_path": str(summary_path),
        "validation": validation,
    }
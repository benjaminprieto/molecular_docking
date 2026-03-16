#!/usr/bin/env python3
"""
Vina Preparation - Core Module (02a)
======================================
Prepares input files for AutoDock Vina docking.

Converts receptor PDB → PDBQT and ligand SDF → PDBQT using Meeko,
calculates the binding box, and generates vina_inputs.json for 02b.

FORMAT CONVERSION:
  Receptor: mk_prepare_receptor.py (Meeko) or prepare_receptor (ADFRsuite)
  Ligands:  mk_prepare_ligand.py (Meeko) — SDF from 00c → PDBQT

BINDING BOX:
  Priority 1: Reference ligand coordinates (± padding)
  Priority 2: Binding site center from 00e (± fixed size)
  Priority 3: User-specified coordinates from campaign_config

Location: 01_src/molecular_docking/m02_vina/vina_preparation.py

Project: molecular_docking
Module: 02a (core)
Version: 1.0
"""

import logging
import json
import subprocess
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Union
from dataclasses import dataclass, field, asdict

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


# =============================================================================
# DEPENDENCY CHECKS
# =============================================================================

def _check_tool(cmd: List[str], timeout: int = 10) -> Tuple[bool, str]:
    """Check if a command-line tool is available."""
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        if result.returncode == 0:
            version = (result.stdout.strip() or result.stderr.strip())[:120]
            return True, version
    except (subprocess.SubprocessError, FileNotFoundError):
        pass
    return False, ""


def check_meeko() -> Tuple[bool, str]:
    """Check if Meeko mk_prepare_ligand.py is available."""
    ok, ver = _check_tool(["mk_prepare_ligand.py", "--help"])
    if ok:
        return True, f"Meeko mk_prepare_ligand.py found"
    # Try as python module
    ok2, ver2 = _check_tool(["python", "-m", "meeko.cli.mk_prepare_ligand", "--help"])
    if ok2:
        return True, f"Meeko found (python -m)"
    return False, "Meeko mk_prepare_ligand.py not found"


def check_meeko_receptor() -> Tuple[bool, str]:
    """Check if Meeko mk_prepare_receptor.py is available."""
    ok, ver = _check_tool(["mk_prepare_receptor.py", "--help"])
    if ok:
        return True, "mk_prepare_receptor.py found"
    return False, "mk_prepare_receptor.py not found"


def check_prepare_receptor() -> Tuple[bool, str]:
    """Check if ADFRsuite prepare_receptor is available."""
    ok, ver = _check_tool(["prepare_receptor", "--help"])
    if ok:
        return True, "prepare_receptor (ADFRsuite) found"
    return False, "prepare_receptor not found"


def check_openbabel() -> Tuple[bool, str]:
    """Check if OpenBabel is available."""
    try:
        result = subprocess.run(
            ['obabel', '-V'], capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            version = result.stdout.strip() or result.stderr.strip()
            return True, f"OpenBabel found: {version}"
    except (subprocess.SubprocessError, FileNotFoundError):
        pass
    return False, "OpenBabel not found"


# Lazy checks — evaluated at import time
MEEKO_AVAILABLE, MEEKO_MSG = check_meeko()
MEEKO_RECEPTOR_AVAILABLE, MEEKO_RECEPTOR_MSG = check_meeko_receptor()
ADFR_AVAILABLE, ADFR_MSG = check_prepare_receptor()
OPENBABEL_AVAILABLE, OPENBABEL_MSG = check_openbabel()


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class VinaBindingBox:
    """3D binding box for Vina docking."""
    center_x: float
    center_y: float
    center_z: float
    size_x: float = 25.0
    size_y: float = 25.0
    size_z: float = 25.0
    source: str = "unknown"
    padding: float = 6.0

    @property
    def center(self) -> Tuple[float, float, float]:
        return (self.center_x, self.center_y, self.center_z)

    @property
    def size(self) -> Tuple[float, float, float]:
        return (self.size_x, self.size_y, self.size_z)

    @property
    def volume(self) -> float:
        return self.size_x * self.size_y * self.size_z

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    def to_vina_args(self) -> List[str]:
        """Generate Vina command-line arguments."""
        return [
            "--center_x", f"{self.center_x:.3f}",
            "--center_y", f"{self.center_y:.3f}",
            "--center_z", f"{self.center_z:.3f}",
            "--size_x", f"{self.size_x:.1f}",
            "--size_y", f"{self.size_y:.1f}",
            "--size_z", f"{self.size_z:.1f}",
        ]


@dataclass
class VinaInput:
    """Container for one molecule's Vina input files."""
    name: str
    receptor_pdbqt: str
    ligand_pdbqt: str
    binding_box: Optional[VinaBindingBox] = None
    ligand_source_sdf: Optional[str] = None
    smiles: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def _to_native_type(value):
    """Convert numpy/pandas types to native Python types."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    if hasattr(value, 'item'):          # numpy scalar
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    return value


# =============================================================================
# RECEPTOR PREPARATION: PDB → PDBQT
# =============================================================================

def prepare_receptor_pdbqt(
        pdb_path: str,
        output_pdbqt: str,
        remove_water: bool = True,
        remove_hetatm: bool = True,
) -> Tuple[bool, str]:
    """
    Convert receptor PDB to PDBQT.

    Tries in order:
      1. mk_prepare_receptor.py (Meeko >= 0.5)
      2. prepare_receptor (ADFRsuite)
      3. obabel fallback

    Args:
        pdb_path:      Input PDB file.
        output_pdbqt:  Output PDBQT file.
        remove_water:  Remove water molecules.
        remove_hetatm: Remove heteroatoms.

    Returns:
        Tuple of (success, message).
    """
    Path(output_pdbqt).parent.mkdir(parents=True, exist_ok=True)

    # --- Option 1: Meeko mk_prepare_receptor.py ---
    if MEEKO_RECEPTOR_AVAILABLE:
        try:
            cmd = [
                "mk_prepare_receptor.py",
                "-i", str(pdb_path),
                "-o", str(output_pdbqt),
            ]
            logger.info(f"  Receptor PDBQT via Meeko: {' '.join(cmd)}")
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=120
            )
            if result.returncode == 0 and Path(output_pdbqt).exists():
                size_kb = Path(output_pdbqt).stat().st_size / 1024
                return True, f"Meeko mk_prepare_receptor.py ({size_kb:.0f} KB)"
            else:
                logger.warning(f"  Meeko receptor failed: {result.stderr[:200]}")
        except Exception as e:
            logger.warning(f"  Meeko receptor error: {e}")

    # --- Option 2: ADFRsuite prepare_receptor ---
    if ADFR_AVAILABLE:
        try:
            cmd = ["prepare_receptor", "-r", str(pdb_path), "-o", str(output_pdbqt)]
            if remove_water:
                cmd.append("-U")
                cmd.append("waters")
            logger.info(f"  Receptor PDBQT via ADFRsuite: {' '.join(cmd)}")
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=120
            )
            if result.returncode == 0 and Path(output_pdbqt).exists():
                size_kb = Path(output_pdbqt).stat().st_size / 1024
                return True, f"ADFRsuite prepare_receptor ({size_kb:.0f} KB)"
            else:
                logger.warning(f"  ADFRsuite failed: {result.stderr[:200]}")
        except Exception as e:
            logger.warning(f"  ADFRsuite error: {e}")

    # --- Option 3: OpenBabel fallback ---
    if OPENBABEL_AVAILABLE:
        try:
            cmd = [
                "obabel", str(pdb_path), "-O", str(output_pdbqt),
                "-xr",  # receptor mode (no torsions)
            ]
            if remove_water:
                cmd.append("-d")  # delete hydrogens first, re-add
            logger.info(f"  Receptor PDBQT via OpenBabel: {' '.join(cmd)}")
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=120
            )
            if result.returncode == 0 and Path(output_pdbqt).exists():
                size_kb = Path(output_pdbqt).stat().st_size / 1024
                return True, f"OpenBabel fallback ({size_kb:.0f} KB)"
            else:
                logger.warning(f"  OpenBabel failed: {result.stderr[:200]}")
        except Exception as e:
            logger.warning(f"  OpenBabel error: {e}")

    return False, "No tool available for receptor PDBQT conversion"


# =============================================================================
# LIGAND PREPARATION: SDF → PDBQT (Meeko)
# =============================================================================

def prepare_ligand_pdbqt_meeko(
        sdf_path: str,
        output_pdbqt: str,
        add_hydrogens: bool = False,
        pH: Optional[float] = None,
) -> Tuple[bool, str]:
    """
    Convert ligand SDF to PDBQT using Meeko mk_prepare_ligand.py.

    Meeko handles:
      - Torsion tree assignment
      - Atom typing for AutoDock
      - Charge assignment (Gasteiger)

    Args:
        sdf_path:      Input SDF file (3D, protonated from 00c).
        output_pdbqt:  Output PDBQT file.
        add_hydrogens: Add hydrogens (False if already protonated by 00c).
        pH:            Protonate at this pH (None if already done).

    Returns:
        Tuple of (success, error_message).
    """
    Path(output_pdbqt).parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "mk_prepare_ligand.py",
        "-i", str(sdf_path),
        "-o", str(output_pdbqt),
    ]
    if add_hydrogens:
        cmd.append("--add_h")
    if pH is not None:
        cmd.extend(["--pH", str(pH)])

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=120
        )
        if result.returncode == 0 and Path(output_pdbqt).exists():
            return True, ""
        else:
            err = (result.stderr or result.stdout)[:300].strip()
            return False, f"Meeko error: {err}"
    except subprocess.TimeoutExpired:
        return False, "Meeko timeout (120s)"
    except Exception as e:
        return False, f"Meeko exception: {e}"


def prepare_ligand_pdbqt_obabel(
        sdf_path: str,
        output_pdbqt: str,
) -> Tuple[bool, str]:
    """Fallback: convert SDF → PDBQT with OpenBabel."""
    if not OPENBABEL_AVAILABLE:
        return False, "OpenBabel not available"

    try:
        cmd = ["obabel", str(sdf_path), "-O", str(output_pdbqt), "--partialcharge", "gasteiger"]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if result.returncode == 0 and Path(output_pdbqt).exists():
            return True, ""
        return False, f"obabel error: {(result.stderr or '')[:200]}"
    except Exception as e:
        return False, str(e)


def batch_prepare_ligands_pdbqt(
        sdf_dir: str,
        output_dir: str,
        molecule_names: List[str],
        add_hydrogens: bool = False,
        pH: Optional[float] = None,
) -> Dict[str, Dict[str, Any]]:
    """
    Batch convert SDF files to PDBQT using Meeko.

    Args:
        sdf_dir:        Directory with individual SDF files from 00c.
        output_dir:     Output directory for PDBQT files.
        molecule_names: List of molecule names to process.
        add_hydrogens:  Add hydrogens (typically False — 00c already protonated).
        pH:             Protonate at this pH (None if already done).

    Returns:
        Dict: name → {pdbqt: path, success: bool, error: str}
    """
    sdf_path = Path(sdf_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    results = {}
    n_ok = 0
    n_fail = 0

    logger.info(f"Converting {len(molecule_names)} ligands: SDF → PDBQT (Meeko)...")

    for name in molecule_names:
        sdf_file = sdf_path / f"{name}.sdf"

        if not sdf_file.exists():
            results[name] = {
                "pdbqt": None, "success": False,
                "error": f"SDF not found: {sdf_file}"
            }
            n_fail += 1
            continue

        pdbqt_file = output_path / f"{name}.pdbqt"
        ok, err = prepare_ligand_pdbqt_meeko(
            str(sdf_file), str(pdbqt_file),
            add_hydrogens=add_hydrogens, pH=pH,
        )

        if not ok:
            # Try obabel fallback
            logger.debug(f"  Meeko failed for {name}, trying obabel: {err}")
            ok, err = prepare_ligand_pdbqt_obabel(str(sdf_file), str(pdbqt_file))

        if ok:
            results[name] = {"pdbqt": str(pdbqt_file), "success": True, "error": ""}
            n_ok += 1
        else:
            results[name] = {"pdbqt": None, "success": False, "error": err}
            n_fail += 1
            logger.warning(f"  ✗ {name}: {err}")

    logger.info(f"  Converted: {n_ok}, Failed: {n_fail}")
    return results


# =============================================================================
# BINDING BOX CALCULATION
# =============================================================================

def calculate_binding_box_from_mol2(
        mol2_path: str,
        padding: float = 6.0,
) -> Optional[VinaBindingBox]:
    """
    Calculate binding box from a reference ligand mol2 file.

    Reads ATOM/HETATM coordinates, computes bounding box + padding.
    """
    coords = _extract_coords_mol2(mol2_path)
    if not coords:
        return None

    arr = np.array(coords)
    min_xyz = arr.min(axis=0)
    max_xyz = arr.max(axis=0)
    center = (min_xyz + max_xyz) / 2.0
    size = (max_xyz - min_xyz) + 2 * padding

    # Enforce minimum size
    size = np.maximum(size, 10.0)

    return VinaBindingBox(
        center_x=float(center[0]),
        center_y=float(center[1]),
        center_z=float(center[2]),
        size_x=float(size[0]),
        size_y=float(size[1]),
        size_z=float(size[2]),
        source="reference_ligand",
        padding=padding,
    )


def calculate_binding_box_from_pdb(
        pdb_path: str,
        padding: float = 6.0,
) -> Optional[VinaBindingBox]:
    """Calculate binding box from a PDB/PDBQT file (e.g., reference ligand)."""
    coords = _extract_coords_pdb(pdb_path)
    if not coords:
        return None

    arr = np.array(coords)
    min_xyz = arr.min(axis=0)
    max_xyz = arr.max(axis=0)
    center = (min_xyz + max_xyz) / 2.0
    size = (max_xyz - min_xyz) + 2 * padding
    size = np.maximum(size, 10.0)

    return VinaBindingBox(
        center_x=float(center[0]),
        center_y=float(center[1]),
        center_z=float(center[2]),
        size_x=float(size[0]),
        size_y=float(size[1]),
        size_z=float(size[2]),
        source="reference_ligand_pdb",
        padding=padding,
    )


def calculate_binding_box_from_sdf(
        sdf_path: str,
        padding: float = 6.0,
) -> Optional[VinaBindingBox]:
    """Calculate binding box from an SDF file."""
    coords = _extract_coords_sdf(sdf_path)
    if not coords:
        return None

    arr = np.array(coords)
    min_xyz = arr.min(axis=0)
    max_xyz = arr.max(axis=0)
    center = (min_xyz + max_xyz) / 2.0
    size = (max_xyz - min_xyz) + 2 * padding
    size = np.maximum(size, 10.0)

    return VinaBindingBox(
        center_x=float(center[0]),
        center_y=float(center[1]),
        center_z=float(center[2]),
        size_x=float(size[0]),
        size_y=float(size[1]),
        size_z=float(size[2]),
        source="reference_ligand_sdf",
        padding=padding,
    )


def make_binding_box_from_center(
        center: Tuple[float, float, float],
        size: Tuple[float, float, float] = (25.0, 25.0, 25.0),
        source: str = "user_coordinates",
) -> VinaBindingBox:
    """Create binding box from explicit center + size."""
    return VinaBindingBox(
        center_x=center[0], center_y=center[1], center_z=center[2],
        size_x=size[0], size_y=size[1], size_z=size[2],
        source=source, padding=0.0,
    )


# =============================================================================
# COORDINATE EXTRACTION HELPERS
# =============================================================================

def _extract_coords_mol2(mol2_path: str) -> List[Tuple[float, float, float]]:
    """Extract atom coordinates from mol2 file."""
    coords = []
    in_atoms = False
    try:
        with open(mol2_path, 'r') as f:
            for line in f:
                stripped = line.strip()
                if stripped.startswith("@<TRIPOS>ATOM"):
                    in_atoms = True
                    continue
                if stripped.startswith("@<TRIPOS>") and in_atoms:
                    break
                if in_atoms and stripped:
                    parts = stripped.split()
                    if len(parts) >= 5:
                        try:
                            x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                            coords.append((x, y, z))
                        except ValueError:
                            pass
    except Exception as e:
        logger.debug(f"Error reading mol2 {mol2_path}: {e}")
    return coords


def _extract_coords_pdb(pdb_path: str) -> List[Tuple[float, float, float]]:
    """Extract atom coordinates from PDB/PDBQT file."""
    coords = []
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append((x, y, z))
                    except (ValueError, IndexError):
                        pass
    except Exception as e:
        logger.debug(f"Error reading PDB {pdb_path}: {e}")
    return coords


def _extract_coords_sdf(sdf_path: str) -> List[Tuple[float, float, float]]:
    """Extract atom coordinates from SDF file (first molecule only)."""
    coords = []
    try:
        with open(sdf_path, 'r') as f:
            lines = f.readlines()
        if len(lines) < 4:
            return coords
        # Counts line is line 3 (0-indexed)
        counts = lines[3].strip().split()
        n_atoms = int(counts[0])
        for i in range(4, 4 + n_atoms):
            parts = lines[i].strip().split()
            if len(parts) >= 3:
                x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                coords.append((x, y, z))
    except Exception as e:
        logger.debug(f"Error reading SDF {sdf_path}: {e}")
    return coords


# =============================================================================
# REFERENCE LIGAND RESOLUTION
# =============================================================================

def find_reference_ligand(
        campaign_config: Dict[str, Any],
        campaign_dir: Path,
) -> Optional[str]:
    """
    Find the reference ligand for binding box calculation.

    Searches in order:
      1. campaign_config → grids.binding_site.reference_mol2
      2. campaign_config → grids.binding_site.reference_sdf
      3. molecules/ directory for common names (UDX, ligand, reference)

    Returns:
        Absolute path to reference ligand, or None.
    """
    grids = campaign_config.get("grids", {})
    bs = grids.get("binding_site", {})

    # Priority 1: reference_mol2 from campaign config
    ref_mol2 = bs.get("reference_mol2")
    if ref_mol2:
        ref_path = Path(ref_mol2) if Path(ref_mol2).is_absolute() else campaign_dir / ref_mol2
        if ref_path.exists():
            logger.info(f"  Reference ligand (mol2): {ref_path}")
            return str(ref_path)

    # Priority 2: reference_sdf
    ref_sdf = bs.get("reference_sdf")
    if ref_sdf:
        ref_path = Path(ref_sdf) if Path(ref_sdf).is_absolute() else campaign_dir / ref_sdf
        if ref_path.exists():
            logger.info(f"  Reference ligand (sdf): {ref_path}")
            return str(ref_path)

    # Priority 3: common names in molecules dir
    mol_dir = campaign_dir / campaign_config.get("molecules", {}).get("input_file", "molecules")
    if mol_dir.is_dir():
        for pattern in ["UDX.*", "reference.*", "ligand.*", "native.*"]:
            matches = list(mol_dir.glob(pattern))
            if matches:
                logger.info(f"  Reference ligand (auto-detected): {matches[0]}")
                return str(matches[0])

    return None


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def run_vina_preparation(
        receptor_pdb: str,
        sdf_dir: str,
        output_dir: str,
        molecule_names: List[str],
        # Binding box
        reference_ligand: Optional[str] = None,
        binding_box_padding: float = 6.0,
        box_center: Optional[Tuple[float, float, float]] = None,
        box_size: Optional[Tuple[float, float, float]] = None,
        # Receptor
        remove_water: bool = True,
        remove_hetatm: bool = True,
        # Ligands
        add_hydrogens: bool = False,
        pH: Optional[float] = None,
        # Reference control
        include_reference: bool = True,
        # Metadata
        molecules_csv: Optional[str] = None,
        name_column: str = "Name",
        smiles_column: str = "SMILES_mol2",
) -> Dict[str, Any]:
    """
    Run the complete Vina preparation pipeline.

    Pipeline:
      1. Convert receptor PDB → PDBQT
      2. Convert ligand SDFs → PDBQTs (Meeko)
      3. Calculate binding box
      4. Generate vina_inputs.json for 02b

    Args:
        receptor_pdb:       Input receptor PDB (from 00b rec_noH.pdb).
        sdf_dir:            Directory with protonated SDF files (from 00c).
        output_dir:         Output directory.
        molecule_names:     List of molecule names to process.
        reference_ligand:   Path to reference ligand for binding box.
        binding_box_padding: Padding around reference ligand (Å).
        box_center:         Explicit box center (x, y, z).
        box_size:           Explicit box size (sx, sy, sz).
        remove_water:       Remove water from receptor.
        remove_hetatm:      Remove heteroatoms from receptor.
        add_hydrogens:      Add H to ligands (False if 00c handled it).
        pH:                 Protonate ligands at this pH.
        include_reference:  Add reference ligand as control molecule.
        molecules_csv:      CSV with molecule metadata (for enriching vina_inputs).
        name_column:        Column name for molecule names in CSV.
        smiles_column:      Column name for SMILES in CSV.

    Returns:
        Dict with summary, paths, and statistics.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    pdbqt_dir = output_path / "ligands_pdbqt"
    pdbqt_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("VINA PREPARATION (02a) v1.0")
    logger.info("=" * 60)

    # Load molecule metadata if available
    df = None
    if molecules_csv and Path(molecules_csv).exists():
        df = pd.read_csv(molecules_csv)
        logger.info(f"Loaded metadata: {len(df)} molecules from {molecules_csv}")

    # =========================================================================
    # STEP 1: RECEPTOR PDB → PDBQT
    # =========================================================================
    logger.info("")
    logger.info("─── Step 1: Receptor PDBQT ───")
    receptor_pdbqt = output_path / "receptor.pdbqt"
    ok, msg = prepare_receptor_pdbqt(
        receptor_pdb, str(receptor_pdbqt),
        remove_water=remove_water,
        remove_hetatm=remove_hetatm,
    )
    if not ok:
        logger.error(f"Receptor PDBQT conversion failed: {msg}")
        return {"success": False, "error": msg}
    logger.info(f"  ✓ {msg}")

    # =========================================================================
    # STEP 2: LIGAND SDFs → PDBQTs
    # =========================================================================
    logger.info("")
    logger.info("─── Step 2: Ligand PDBQTs (Meeko) ───")
    ligand_results = batch_prepare_ligands_pdbqt(
        sdf_dir=sdf_dir,
        output_dir=str(pdbqt_dir),
        molecule_names=molecule_names,
        add_hydrogens=add_hydrogens,
        pH=pH,
    )

    # =========================================================================
    # STEP 3: BINDING BOX
    # =========================================================================
    logger.info("")
    logger.info("─── Step 3: Binding Box ───")
    binding_box = None

    # Priority 1: Explicit coordinates
    if box_center is not None:
        bsize = box_size or (25.0, 25.0, 25.0)
        binding_box = make_binding_box_from_center(box_center, bsize, source="user_coordinates")
        logger.info(f"  Using explicit coordinates: center={box_center}")

    # Priority 2: Reference ligand
    if binding_box is None and reference_ligand and Path(reference_ligand).exists():
        ext = Path(reference_ligand).suffix.lower()
        if ext == ".mol2":
            binding_box = calculate_binding_box_from_mol2(reference_ligand, binding_box_padding)
        elif ext == ".sdf":
            binding_box = calculate_binding_box_from_sdf(reference_ligand, binding_box_padding)
        elif ext in (".pdb", ".pdbqt"):
            binding_box = calculate_binding_box_from_pdb(reference_ligand, binding_box_padding)

        if binding_box:
            binding_box.source = "reference_ligand"
            logger.info(f"  From reference ligand: {Path(reference_ligand).name}")

    if binding_box:
        logger.info(f"  Center: ({binding_box.center_x:.2f}, "
                     f"{binding_box.center_y:.2f}, {binding_box.center_z:.2f})")
        logger.info(f"  Size:   ({binding_box.size_x:.1f}, "
                     f"{binding_box.size_y:.1f}, {binding_box.size_z:.1f})")
        logger.info(f"  Volume: {binding_box.volume:,.0f} Å³")
        logger.info(f"  Source: {binding_box.source}")
    else:
        logger.error("  Could not determine binding box!")
        return {"success": False, "error": "No binding box could be calculated"}

    # =========================================================================
    # STEP 4: REFERENCE LIGAND AS CONTROL
    # =========================================================================
    ref_pdbqt = None
    if include_reference and reference_ligand and Path(reference_ligand).exists():
        logger.info("")
        logger.info("─── Step 4: Reference Ligand Control ───")
        ref_name = Path(reference_ligand).stem.replace('_aligned', '').replace(
            '_interactions', '').upper()
        ref_pdbqt_path = pdbqt_dir / f"{ref_name}.pdbqt"

        # Convert reference to PDBQT
        ext = Path(reference_ligand).suffix.lower()
        if ext == ".mol2":
            # mol2 → SDF → PDBQT (Meeko needs SDF input)
            temp_sdf = output_path / f"_temp_{ref_name}.sdf"
            if OPENBABEL_AVAILABLE:
                subprocess.run(
                    ["obabel", reference_ligand, "-O", str(temp_sdf)],
                    capture_output=True, timeout=60
                )
                if temp_sdf.exists():
                    ok_ref, err_ref = prepare_ligand_pdbqt_meeko(
                        str(temp_sdf), str(ref_pdbqt_path)
                    )
                    if not ok_ref:
                        ok_ref, err_ref = prepare_ligand_pdbqt_obabel(
                            str(temp_sdf), str(ref_pdbqt_path)
                        )
                    temp_sdf.unlink(missing_ok=True)
                    if ok_ref:
                        ref_pdbqt = str(ref_pdbqt_path)
                        logger.info(f"  ✓ Reference {ref_name} → {ref_pdbqt_path.name}")
        elif ext == ".sdf":
            ok_ref, err_ref = prepare_ligand_pdbqt_meeko(
                reference_ligand, str(ref_pdbqt_path)
            )
            if ok_ref:
                ref_pdbqt = str(ref_pdbqt_path)
                logger.info(f"  ✓ Reference {ref_name} → {ref_pdbqt_path.name}")

    # =========================================================================
    # STEP 5: GENERATE VINA INPUTS JSON
    # =========================================================================
    logger.info("")
    logger.info("─── Step 5: Vina Inputs JSON ───")

    vina_inputs = []
    for name in molecule_names:
        lr = ligand_results.get(name, {})
        if not lr.get("success"):
            continue

        # Extract metadata from CSV
        metadata = {}
        smiles_val = None
        if df is not None and name_column in df.columns:
            rows = df[df[name_column] == name]
            if not rows.empty:
                row = rows.iloc[0]
                if smiles_column in df.columns:
                    val = row.get(smiles_column)
                    smiles_val = str(val) if pd.notna(val) else None
                for col in ['Grid_Score', 'MW', 'LogP', 'HBD', 'HBA']:
                    if col in df.columns:
                        metadata[col] = _to_native_type(row.get(col))

        vina_inputs.append(VinaInput(
            name=name,
            receptor_pdbqt=str(receptor_pdbqt),
            ligand_pdbqt=lr["pdbqt"],
            binding_box=binding_box,
            ligand_source_sdf=str(Path(sdf_dir) / f"{name}.sdf"),
            smiles=smiles_val,
            metadata=metadata,
        ))

    # Add reference as control
    if ref_pdbqt and include_reference:
        ref_name = Path(ref_pdbqt).stem
        vina_inputs.append(VinaInput(
            name=ref_name,
            receptor_pdbqt=str(receptor_pdbqt),
            ligand_pdbqt=ref_pdbqt,
            binding_box=binding_box,
            metadata={"Strategy": "reference_control", "is_reference": True},
        ))

    # Serialize to JSON
    inputs_file = output_path / "vina_inputs.json"
    inputs_data = {
        "version": "1.0",
        "engine": "vina",
        "receptor_pdbqt": str(receptor_pdbqt),
        "binding_box": binding_box.to_dict(),
        "molecules": [
            {
                "name": vi.name,
                "receptor_pdbqt": vi.receptor_pdbqt,
                "ligand_pdbqt": vi.ligand_pdbqt,
                "ligand_source_sdf": vi.ligand_source_sdf,
                "smiles": vi.smiles,
                "binding_box": vi.binding_box.to_dict() if vi.binding_box else None,
                "metadata": vi.metadata,
            }
            for vi in vina_inputs
        ],
    }
    with open(inputs_file, 'w') as f:
        json.dump(inputs_data, f, indent=2)
    logger.info(f"  Saved: {inputs_file}")

    # Save binding box JSON separately
    box_file = output_path / "binding_box.json"
    with open(box_file, 'w') as f:
        json.dump(binding_box.to_dict(), f, indent=2)

    # =========================================================================
    # SUMMARY
    # =========================================================================
    n_prepared = sum(1 for r in ligand_results.values() if r["success"])
    n_failed = sum(1 for r in ligand_results.values() if not r["success"])

    logger.info("")
    logger.info("=" * 60)
    logger.info("VINA PREPARATION COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Receptor PDBQT:     {receptor_pdbqt}")
    logger.info(f"Ligands prepared:   {n_prepared}/{len(molecule_names)}")
    logger.info(f"Ligands failed:     {n_failed}")
    logger.info(f"Vina inputs JSON:   {inputs_file}")
    logger.info(f"Binding box source: {binding_box.source}")

    return {
        "success": True,
        "n_molecules": len(molecule_names),
        "n_prepared": n_prepared,
        "n_failed": n_failed,
        "receptor_pdbqt": str(receptor_pdbqt),
        "pdbqt_dir": str(pdbqt_dir),
        "binding_box": binding_box.to_dict(),
        "inputs_json": str(inputs_file),
        "ligand_results": ligand_results,
        "vina_inputs": vina_inputs,
    }


# =============================================================================
# UTILITY
# =============================================================================

def load_vina_inputs(inputs_json: str) -> List[Dict[str, Any]]:
    """Load vina_inputs.json (used by 02b)."""
    with open(inputs_json, 'r') as f:
        data = json.load(f)

    if isinstance(data, dict):
        molecules = data.get('molecules', [])
        top_box = data.get('binding_box')
        if top_box:
            for mol in molecules:
                if not mol.get('binding_box'):
                    mol['binding_box'] = top_box
        return molecules
    elif isinstance(data, list):
        return data
    else:
        raise ValueError(f"Unexpected format in {inputs_json}")


if __name__ == '__main__':
    print("Vina Preparation - Core Module (02a) v1.0")
    print("Use 02a_vina_preparation.py CLI for execution")
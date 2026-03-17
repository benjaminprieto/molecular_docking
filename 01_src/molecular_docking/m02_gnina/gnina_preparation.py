#!/usr/bin/env python3
"""
GNINA Preparation - Core Module (02a)
=======================================
Prepares input files for GNINA docking.

GNINA accepts PDB + SDF/mol2 directly — no format conversion needed.
This module:
  1. Cleans receptor PDB (optional)
  2. Calculates binding box from reference ligand
  3. Generates gnina_inputs.json for 02b

Much simpler than Vina preparation (no PDBQT/Meeko dependency).

Location: 01_src/molecular_docking/m02_gnina/gnina_preparation.py

Project: molecular_docking
Module: 02a (core)
Version: 1.0
"""

import logging
import json
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Union
from dataclasses import dataclass, field, asdict

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class GninaBindingBox:
    """3D binding box for GNINA docking."""
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

    def to_gnina_args(self) -> List[str]:
        """Generate GNINA command-line arguments."""
        return [
            "--center_x", f"{self.center_x:.3f}",
            "--center_y", f"{self.center_y:.3f}",
            "--center_z", f"{self.center_z:.3f}",
            "--size_x", f"{self.size_x:.1f}",
            "--size_y", f"{self.size_y:.1f}",
            "--size_z", f"{self.size_z:.1f}",
        ]


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def _to_native_type(value):
    """Convert numpy/pandas types to native Python types."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    if hasattr(value, 'item'):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    return value


# =============================================================================
# BINDING BOX CALCULATION
# =============================================================================

def calculate_binding_box_from_mol2(
        mol2_path: str,
        padding: float = 6.0,
) -> Optional[GninaBindingBox]:
    """Calculate binding box from a reference ligand mol2 file."""
    coords = _extract_coords_mol2(mol2_path)
    if not coords:
        return None
    return _box_from_coords(coords, padding, "reference_ligand_mol2")


def calculate_binding_box_from_sdf(
        sdf_path: str,
        padding: float = 6.0,
) -> Optional[GninaBindingBox]:
    """Calculate binding box from an SDF file."""
    coords = _extract_coords_sdf(sdf_path)
    if not coords:
        return None
    return _box_from_coords(coords, padding, "reference_ligand_sdf")


def calculate_binding_box_from_pdb(
        pdb_path: str,
        padding: float = 6.0,
) -> Optional[GninaBindingBox]:
    """Calculate binding box from a PDB file."""
    coords = _extract_coords_pdb(pdb_path)
    if not coords:
        return None
    return _box_from_coords(coords, padding, "reference_ligand_pdb")


def _box_from_coords(
        coords: List[Tuple[float, float, float]],
        padding: float,
        source: str,
) -> GninaBindingBox:
    """Build a binding box from coordinate list."""
    arr = np.array(coords)
    min_xyz = arr.min(axis=0)
    max_xyz = arr.max(axis=0)
    center = (min_xyz + max_xyz) / 2.0
    size = (max_xyz - min_xyz) + 2 * padding
    size = np.maximum(size, 10.0)

    return GninaBindingBox(
        center_x=float(center[0]),
        center_y=float(center[1]),
        center_z=float(center[2]),
        size_x=float(size[0]),
        size_y=float(size[1]),
        size_z=float(size[2]),
        source=source,
        padding=padding,
    )


# =============================================================================
# COORDINATE EXTRACTION
# =============================================================================

def _extract_coords_mol2(mol2_path: str) -> List[Tuple[float, float, float]]:
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
                            coords.append((float(parts[2]), float(parts[3]), float(parts[4])))
                        except ValueError:
                            pass
    except Exception as e:
        logger.debug(f"Error reading mol2 {mol2_path}: {e}")
    return coords


def _extract_coords_sdf(sdf_path: str) -> List[Tuple[float, float, float]]:
    coords = []
    try:
        with open(sdf_path, 'r') as f:
            lines = f.readlines()
        if len(lines) < 4:
            return coords
        counts = lines[3].strip().split()
        n_atoms = int(counts[0])
        for i in range(4, 4 + n_atoms):
            parts = lines[i].strip().split()
            if len(parts) >= 3:
                coords.append((float(parts[0]), float(parts[1]), float(parts[2])))
    except Exception as e:
        logger.debug(f"Error reading SDF {sdf_path}: {e}")
    return coords


def _extract_coords_pdb(pdb_path: str) -> List[Tuple[float, float, float]]:
    coords = []
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                    except (ValueError, IndexError):
                        pass
    except Exception as e:
        logger.debug(f"Error reading PDB {pdb_path}: {e}")
    return coords


# =============================================================================
# REFERENCE LIGAND RESOLUTION
# =============================================================================

def find_reference_ligand(
        campaign_config: Dict[str, Any],
        campaign_dir: Path,
) -> Optional[str]:
    """Find the reference ligand for binding box calculation."""
    grids = campaign_config.get("grids", {})
    bs = grids.get("binding_site", {})

    for key in ["reference_mol2", "reference_sdf"]:
        ref = bs.get(key)
        if ref:
            ref_path = Path(ref) if Path(ref).is_absolute() else campaign_dir / ref
            if ref_path.exists():
                logger.info(f"  Reference ligand ({key}): {ref_path}")
                return str(ref_path)

    # Auto-detect in molecules dir
    mol_dir = campaign_dir / campaign_config.get("molecules", {}).get("input_file", "molecules")
    if mol_dir.is_dir():
        for pattern in ["UDX.*", "reference.*", "ligand.*", "native.*"]:
            matches = list(mol_dir.glob(pattern))
            if matches:
                logger.info(f"  Reference ligand (auto-detected): {matches[0]}")
                return str(matches[0])

    return None


# =============================================================================
# PROTEIN CLEANING
# =============================================================================

def clean_protein_pdb(
        input_pdb: str,
        output_pdb: str,
        remove_water: bool = True,
        remove_hetatm: bool = True,
) -> int:
    """Clean PDB for GNINA (remove water/HETATM, keep ATOM records)."""
    Path(output_pdb).parent.mkdir(parents=True, exist_ok=True)
    n_kept = 0

    with open(input_pdb, 'r') as fin, open(output_pdb, 'w') as fout:
        for line in fin:
            if line.startswith("ATOM"):
                fout.write(line)
                n_kept += 1
            elif line.startswith("HETATM"):
                if not remove_hetatm:
                    fout.write(line)
                    n_kept += 1
            elif line.startswith(("TER", "END", "MODEL", "ENDMDL", "CRYST1")):
                fout.write(line)

    return n_kept


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def run_gnina_preparation(
        receptor_pdb: str,
        sdf_dir: str,
        output_dir: str,
        molecule_names: List[str],
        # Binding box
        reference_ligand: Optional[str] = None,
        binding_box_padding: float = 6.0,
        autobox_ligand: Optional[str] = None,
        use_autobox: bool = False,
        box_center: Optional[Tuple[float, float, float]] = None,
        box_size: Optional[Tuple[float, float, float]] = None,
        # Receptor
        clean_protein: bool = True,
        remove_water: bool = True,
        remove_hetatm: bool = True,
        # Ligand format
        ligand_format: str = "sdf",
        # Reference control
        include_reference: bool = True,
        # Metadata
        molecules_csv: Optional[str] = None,
        name_column: str = "Name",
        smiles_column: str = "SMILES_mol2",
) -> Dict[str, Any]:
    """
    Run the complete GNINA preparation pipeline.

    GNINA accepts PDB + SDF/mol2 directly — no format conversion needed.

    Pipeline:
      1. Clean receptor PDB (optional)
      2. Calculate binding box from reference ligand
      3. Locate ligand SDF files from 00c
      4. Generate gnina_inputs.json for 02b

    Returns:
        Dict with summary, paths, and statistics.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("GNINA PREPARATION (02a) v1.0")
    logger.info("=" * 60)

    # Load molecule metadata
    df = None
    if molecules_csv and Path(molecules_csv).exists():
        df = pd.read_csv(molecules_csv)
        logger.info(f"Loaded metadata: {len(df)} molecules from {molecules_csv}")

    # =========================================================================
    # STEP 1: RECEPTOR PDB
    # =========================================================================
    logger.info("")
    logger.info("─── Step 1: Receptor PDB ───")
    protein_pdb = output_path / "protein_clean.pdb"

    if clean_protein:
        n_atoms = clean_protein_pdb(
            receptor_pdb, str(protein_pdb),
            remove_water=remove_water, remove_hetatm=remove_hetatm,
        )
        logger.info(f"  ✓ Cleaned PDB: {n_atoms} atoms → {protein_pdb.name}")
    else:
        shutil.copy2(receptor_pdb, protein_pdb)
        logger.info(f"  ✓ Copied PDB as-is: {protein_pdb.name}")

    # =========================================================================
    # STEP 2: BINDING BOX
    # =========================================================================
    logger.info("")
    logger.info("─── Step 2: Binding Box ───")
    binding_box = None

    # Priority 1: Explicit coordinates
    if box_center is not None:
        bsize = box_size or (25.0, 25.0, 25.0)
        binding_box = GninaBindingBox(
            center_x=box_center[0], center_y=box_center[1], center_z=box_center[2],
            size_x=bsize[0], size_y=bsize[1], size_z=bsize[2],
            source="user_coordinates", padding=0.0,
        )
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
        logger.warning("  Could not determine binding box from coordinates or reference ligand")
        logger.info("  GNINA will require --autobox_ligand at runtime (02b)")

    # =========================================================================
    # STEP 3: LOCATE LIGAND FILES
    # =========================================================================
    logger.info("")
    logger.info("─── Step 3: Ligand Files ───")

    sdf_path = Path(sdf_dir)
    ligand_files = {}
    n_found = 0
    n_missing = 0

    for name in molecule_names:
        sdf_file = sdf_path / f"{name}.sdf"
        if sdf_file.exists():
            ligand_files[name] = str(sdf_file.resolve())
            n_found += 1
        else:
            logger.warning(f"  ✗ SDF not found: {name}.sdf")
            n_missing += 1

    logger.info(f"  Found: {n_found}/{len(molecule_names)} SDF files")
    if n_missing:
        logger.warning(f"  Missing: {n_missing}")

    # =========================================================================
    # STEP 4: REFERENCE LIGAND AS CONTROL
    # =========================================================================
    ref_ligand_control = None
    if include_reference and reference_ligand and Path(reference_ligand).exists():
        logger.info("")
        logger.info("─── Step 4: Reference Ligand Control ───")
        ref_name = Path(reference_ligand).stem.replace('_aligned', '').replace(
            '_interactions', '').upper()
        ref_ligand_control = str(Path(reference_ligand).resolve())
        logger.info(f"  ✓ Reference control: {ref_name} ({Path(reference_ligand).name})")

    # =========================================================================
    # STEP 5: GENERATE GNINA INPUTS JSON
    # =========================================================================
    logger.info("")
    logger.info("─── Step 5: GNINA Inputs JSON ───")

    protein_abs = str(protein_pdb.resolve())
    gnina_inputs = []

    for name in molecule_names:
        if name not in ligand_files:
            continue

        # Extract metadata
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

        gnina_inputs.append({
            "name": name,
            "protein_pdb": protein_abs,
            "ligand_sdf": ligand_files[name],
            "ligand_format": "sdf",
            "smiles": smiles_val,
            "binding_box": binding_box.to_dict() if binding_box else None,
            "metadata": metadata,
        })

    # Add reference as control
    if ref_ligand_control and include_reference:
        ref_name = Path(ref_ligand_control).stem.replace('_aligned', '').replace(
            '_interactions', '').upper()
        ref_format = Path(ref_ligand_control).suffix.lstrip('.').lower()
        gnina_inputs.append({
            "name": ref_name,
            "protein_pdb": protein_abs,
            "ligand_sdf": ref_ligand_control if ref_format == "sdf" else None,
            "ligand_mol2": ref_ligand_control if ref_format == "mol2" else None,
            "ligand_format": ref_format,
            "smiles": None,
            "binding_box": binding_box.to_dict() if binding_box else None,
            "metadata": {"Strategy": "reference_control", "is_reference": True},
        })

    # Serialize
    inputs_data = {
        "version": "1.0",
        "engine": "gnina",
        "protein_pdb": protein_abs,
        "binding_box": binding_box.to_dict() if binding_box else None,
        "autobox_ligand": str(Path(autobox_ligand).resolve()) if autobox_ligand and Path(autobox_ligand).exists() else (
            str(Path(reference_ligand).resolve()) if use_autobox and reference_ligand and Path(reference_ligand).exists() else None
        ),
        "use_autobox": use_autobox,
        "molecules": gnina_inputs,
    }

    inputs_file = output_path / "gnina_inputs.json"
    with open(inputs_file, 'w') as f:
        json.dump(inputs_data, f, indent=2)
    logger.info(f"  Saved: {inputs_file}")

    # Save binding box JSON
    if binding_box:
        box_file = output_path / "binding_box.json"
        with open(box_file, 'w') as f:
            json.dump(binding_box.to_dict(), f, indent=2)

    # =========================================================================
    # SUMMARY
    # =========================================================================
    logger.info("")
    logger.info("=" * 60)
    logger.info("GNINA PREPARATION COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Protein PDB:        {protein_pdb}")
    logger.info(f"Ligands found:      {n_found}/{len(molecule_names)}")
    logger.info(f"GNINA inputs JSON:  {inputs_file}")
    if binding_box:
        logger.info(f"Binding box source: {binding_box.source}")
    logger.info(f"Autobox:            {use_autobox}")

    return {
        "success": True,
        "n_molecules": len(molecule_names),
        "n_prepared": len(gnina_inputs),
        "n_found": n_found,
        "n_missing": n_missing,
        "protein_pdb": str(protein_pdb),
        "binding_box": binding_box.to_dict() if binding_box else None,
        "inputs_json": str(inputs_file),
        "ligand_files": ligand_files,
    }


if __name__ == '__main__':
    print("GNINA Preparation - Core Module (02a) v1.0")
    print("Use 02a_gnina_preparation.py CLI for execution")

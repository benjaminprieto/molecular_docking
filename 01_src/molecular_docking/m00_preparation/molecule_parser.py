"""
Molecule Parser - Core Module (00a)
=====================================
Parsea moleculas de cualquier formato y produce tabla normalizada + SDF limpio.

Formatos soportados:
  - SDF (uno o varios archivos, con propiedades opcionales como <rmsd>)
  - CSV/XLSX con columna SMILES + columna Name
  - Directorio de mol2 individuales
  - SMILES plain text (una por linea: SMILES name)
  - Directorio con mezcla de formatos

Input:
  - Archivo unico (SDF, CSV, XLSX, SMILES txt)
  - Directorio (escanea todos los formatos soportados dentro)
  - Lista de archivos (en campaign_config)

Output:
  - unique_molecules.csv       (Name, SMILES, Heavy_Atoms, MW, LogP, ...)
  - unique_molecules.sdf       (SDF limpio, 1 entrada por molecula)
  - parse_summary.txt

Location: 01_src/molecular_docking/m00_preparation/molecule_parser.py
Project: molecular_docking
Module: 00a (core)
Version: 1.0
"""

import json
import logging
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# =============================================================================
# OPTIONAL DEPENDENCIES
# =============================================================================

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, AllChem, rdMolDescriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logger.warning("RDKit not available — molecule parsing will be limited")


# =============================================================================
# FORMAT DETECTION
# =============================================================================

SUPPORTED_EXTENSIONS = {
    ".sdf": "sdf",
    ".sd":  "sdf",
    ".csv": "csv",
    ".tsv": "csv",
    ".xlsx": "xlsx",
    ".xls": "xlsx",
    ".smi": "smiles",
    ".smiles": "smiles",
    ".txt": "smiles",
    ".mol2": "mol2",
}


def detect_input_format(input_path: Union[str, Path]) -> str:
    """
    Detect input format from file extension or directory contents.

    Returns:
        "sdf" | "csv" | "xlsx" | "smiles" | "mol2" | "mol2_dir" | "directory"
    """
    path = Path(input_path)

    if path.is_dir():
        # Check what's inside
        files = list(path.iterdir())
        extensions = {f.suffix.lower() for f in files if f.is_file()}
        mol2_files = [f for f in files if f.suffix.lower() == ".mol2"]

        if mol2_files and not extensions.intersection({".sdf", ".csv", ".xlsx"}):
            return "mol2_dir"
        return "directory"

    suffix = path.suffix.lower()
    if suffix in SUPPORTED_EXTENSIONS:
        return SUPPORTED_EXTENSIONS[suffix]

    raise ValueError(f"Unsupported file format: {suffix} ({path.name})")


def discover_files_in_directory(dir_path: Path) -> Dict[str, List[Path]]:
    """
    Scan a directory and group files by format.

    Returns:
        {"sdf": [path1, path2], "csv": [path3], "mol2": [path4, ...]}
    """
    groups: Dict[str, List[Path]] = defaultdict(list)

    for f in sorted(dir_path.iterdir()):
        if not f.is_file():
            continue
        suffix = f.suffix.lower()
        if suffix in SUPPORTED_EXTENSIONS:
            fmt = SUPPORTED_EXTENSIONS[suffix]
            groups[fmt].append(f)

    return dict(groups)


# =============================================================================
# SDF PARSING
# =============================================================================

def parse_sdf(
        sdf_path: str,
        name_source: str = "header_first_word",
        name_property: str = "_Name",
) -> List[Dict[str, Any]]:
    """
    Parse an SDF file into a list of molecule dicts.

    Handles Phermit-style headers (multiple IDs separated by spaces)
    and optional properties like <rmsd>.

    Args:
        sdf_path: Path to SDF file
        name_source: "header_first_word" | "header_full" | "property"
        name_property: SDF property for name (if name_source="property")

    Returns:
        List of dicts with: name, mol, smiles, heavy_atoms, rmsd, all_ids, source_file
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit required for SDF parsing")

    # Silence RDKit warnings
    rd_logger = RDLogger.logger()
    rd_logger.setLevel(RDLogger.ERROR)

    path = Path(sdf_path)
    logger.info(f"  Parsing SDF: {path.name}")

    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    molecules = []
    n_failed = 0

    for idx, mol in enumerate(supplier):
        if mol is None:
            n_failed += 1
            continue

        # --- Extract name ---
        header = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""

        if name_source == "header_first_word":
            name = header.split()[0] if header else f"MOL_{idx:04d}"
        elif name_source == "header_full":
            name = header if header else f"MOL_{idx:04d}"
        elif name_source == "property":
            name = mol.GetProp(name_property).strip() if mol.HasProp(name_property) else f"MOL_{idx:04d}"
        else:
            name = header.split()[0] if header else f"MOL_{idx:04d}"

        # --- Extract RMSD (optional, from Phermit) ---
        rmsd = None
        if mol.HasProp("rmsd"):
            try:
                rmsd = float(mol.GetProp("rmsd"))
            except ValueError:
                pass

        # --- Canonical SMILES ---
        mol_noH = Chem.RemoveHs(mol)
        smiles = Chem.MolToSmiles(mol_noH)
        if not smiles:
            n_failed += 1
            continue

        # --- All SDF properties ---
        extra_props = {}
        for prop_name in mol.GetPropsAsDict():
            if prop_name not in ("_Name", "rmsd"):
                extra_props[prop_name] = mol.GetProp(prop_name)

        molecules.append({
            "name": name,
            "mol": mol,
            "mol_noH": mol_noH,
            "smiles": smiles,
            "heavy_atoms": mol_noH.GetNumHeavyAtoms(),
            "rmsd": rmsd,
            "all_ids": header,
            "source_file": path.name,
            "source_format": "sdf",
            "source_sdf_path": str(path.resolve()),
            "entry_idx": idx,
            "extra_props": extra_props,
        })

    logger.info(f"    Parsed: {len(molecules)} molecules, {n_failed} failed")
    return molecules


# =============================================================================
# CSV / XLSX PARSING
# =============================================================================

def parse_tabular(
        file_path: str,
        name_column: Optional[str] = None,
        smiles_column: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """
    Parse a CSV/XLSX/TSV file with SMILES column.

    Auto-detects name and SMILES columns if not specified.
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit required for SMILES validation")

    path = Path(file_path)
    logger.info(f"  Parsing tabular: {path.name}")

    if path.suffix.lower() in (".xlsx", ".xls"):
        df = pd.read_excel(path)
    elif path.suffix.lower() == ".tsv":
        df = pd.read_csv(path, sep="\t")
    else:
        df = pd.read_csv(path)

    # Auto-detect columns
    if smiles_column is None:
        smiles_candidates = ["SMILES", "Smile", "smiles", "Smiles", "smi",
                             "canonical_smiles", "CANONICAL_SMILES"]
        for c in smiles_candidates:
            if c in df.columns:
                smiles_column = c
                break
        if smiles_column is None:
            raise ValueError(
                f"SMILES column not found. Available: {list(df.columns)[:10]}. "
                f"Set smiles_column explicitly."
            )

    if name_column is None:
        name_candidates = ["Name", "name", "NAME", "Molecule_Name", "ID",
                           "id", "compound_id", "CompoundID"]
        for c in name_candidates:
            if c in df.columns:
                name_column = c
                break
        if name_column is None:
            # Use index as name
            name_column = None

    logger.info(f"    Name column: {name_column or '(index)'}")
    logger.info(f"    SMILES column: {smiles_column}")

    rd_logger = RDLogger.logger()
    rd_logger.setLevel(RDLogger.ERROR)

    molecules = []
    n_failed = 0

    for idx, row in df.iterrows():
        smiles_raw = str(row[smiles_column]).strip()
        if not smiles_raw or smiles_raw == "nan":
            n_failed += 1
            continue

        mol = Chem.MolFromSmiles(smiles_raw)
        if mol is None:
            n_failed += 1
            logger.debug(f"    Failed to parse SMILES at row {idx}: {smiles_raw[:50]}")
            continue

        name = str(row[name_column]).strip() if name_column else f"MOL_{idx:04d}"
        smiles = Chem.MolToSmiles(mol)

        molecules.append({
            "name": name,
            "mol": Chem.AddHs(mol),
            "mol_noH": mol,
            "smiles": smiles,
            "heavy_atoms": mol.GetNumHeavyAtoms(),
            "rmsd": None,
            "all_ids": name,
            "source_file": path.name,
            "source_format": "tabular",
            "entry_idx": idx,
            "extra_props": {},
        })

    logger.info(f"    Parsed: {len(molecules)} molecules, {n_failed} failed")
    return molecules


# =============================================================================
# MOL2 PARSING (with OpenBabel fallback for DOCK6/AMBER atom types)
# =============================================================================

def _mol2_via_obabel(mol2_path: str) -> Optional[Dict[str, Any]]:
    """
    Parse a mol2 file using OpenBabel as fallback.

    DOCK6 mol2 files use AMBER atom types (os, c3, n, etc.) that RDKit
    cannot handle. OpenBabel reads them fine and can produce SMILES.

    Returns:
        Dict with name, smiles, heavy_atoms, or None if failed
    """
    import subprocess

    try:
        # Get SMILES via obabel
        result = subprocess.run(
            ["obabel", mol2_path, "-osmi", "-e"],
            capture_output=True, text=True, timeout=30,
        )
        if result.returncode != 0 or not result.stdout.strip():
            return None

        smiles_line = result.stdout.strip().split("\n")[0]
        parts = smiles_line.split()
        if not parts:
            return None

        smiles_raw = parts[0]
        obabel_name = parts[1] if len(parts) > 1 else Path(mol2_path).stem

        # Canonicalize with RDKit if available
        if RDKIT_AVAILABLE:
            mol = Chem.MolFromSmiles(smiles_raw)
            if mol is None:
                return None
            smiles = Chem.MolToSmiles(mol)
            heavy_atoms = mol.GetNumHeavyAtoms()
            mol_with_H = Chem.AddHs(mol)
        else:
            smiles = smiles_raw
            heavy_atoms = sum(1 for c in smiles if c.isupper())  # rough estimate
            mol_with_H = None
            mol = None

        return {
            "smiles": smiles,
            "heavy_atoms": heavy_atoms,
            "mol": mol_with_H,
            "mol_noH": mol,
            "obabel_name": obabel_name,
        }

    except (subprocess.TimeoutExpired, FileNotFoundError, Exception) as e:
        logger.debug(f"    obabel fallback failed for {mol2_path}: {e}")
        return None


def _parse_mol2_name(mol2_path: str) -> str:
    """Extract molecule name from mol2 file (line after @<TRIPOS>MOLECULE)."""
    try:
        with open(mol2_path) as f:
            for line in f:
                if "@<TRIPOS>MOLECULE" in line:
                    name_line = next(f, "").strip()
                    return name_line if name_line else Path(mol2_path).stem
    except Exception:
        pass
    return Path(mol2_path).stem


def _try_parse_mol2_rdkit(mol2_path: str) -> Optional[Any]:
    """Try RDKit mol2 parsing, return None on failure (no crash)."""
    if not RDKIT_AVAILABLE:
        return None
    try:
        rd_logger = RDLogger.logger()
        rd_logger.setLevel(RDLogger.ERROR)
        mol = Chem.MolFromMol2File(mol2_path, removeHs=False, sanitize=False)
        if mol is not None:
            try:
                Chem.SanitizeMol(mol)
            except Exception:
                return None
        return mol
    except Exception:
        return None


def _parse_one_mol2(mol2_path: str, source_file: str, entry_idx: int = 0) -> Optional[Dict[str, Any]]:
    """
    Parse a single mol2 file with RDKit, falling back to OpenBabel.

    Handles DOCK6/AMBER atom types that crash RDKit.
    """
    name = _parse_mol2_name(mol2_path)

    # Try RDKit first
    mol = _try_parse_mol2_rdkit(mol2_path)
    if mol is not None:
        mol_noH = Chem.RemoveHs(mol)
        smiles = Chem.MolToSmiles(mol_noH)
        if smiles:
            return {
                "name": name,
                "mol": mol,
                "mol_noH": mol_noH,
                "smiles": smiles,
                "heavy_atoms": mol_noH.GetNumHeavyAtoms(),
                "rmsd": None,
                "all_ids": name,
                "source_file": source_file,
                "source_format": "mol2",
                "source_mol2_path": str(Path(mol2_path).resolve()),
                "entry_idx": entry_idx,
                "extra_props": {},
            }

    # Fallback to OpenBabel (handles DOCK6/AMBER atom types)
    logger.debug(f"    RDKit failed for {source_file}, trying obabel fallback")
    ob_result = _mol2_via_obabel(mol2_path)
    if ob_result:
        return {
            "name": name,
            "mol": ob_result["mol"],
            "mol_noH": ob_result["mol_noH"],
            "smiles": ob_result["smiles"],
            "heavy_atoms": ob_result["heavy_atoms"],
            "rmsd": None,
            "all_ids": name,
            "source_file": source_file,
            "source_format": "mol2 (obabel)",
            "source_mol2_path": str(Path(mol2_path).resolve()),
            "entry_idx": entry_idx,
            "extra_props": {},
        }

    return None


def parse_mol2_directory(mol2_dir: str) -> List[Dict[str, Any]]:
    """Parse a directory of mol2 files (one molecule per file)."""
    dir_path = Path(mol2_dir)
    mol2_files = sorted(dir_path.glob("*.mol2"))
    logger.info(f"  Parsing mol2 directory: {dir_path.name} ({len(mol2_files)} files)")

    molecules = []
    n_failed = 0

    for mol2_file in mol2_files:
        result = _parse_one_mol2(str(mol2_file), mol2_file.name)
        if result:
            molecules.append(result)
        else:
            n_failed += 1
            logger.debug(f"    Failed: {mol2_file.name}")

    logger.info(f"    Parsed: {len(molecules)} molecules, {n_failed} failed")
    return molecules


def parse_single_mol2(mol2_path: str) -> List[Dict[str, Any]]:
    """Parse a single mol2 file (may contain multiple molecules)."""
    path = Path(mol2_path)
    logger.info(f"  Parsing mol2: {path.name}")

    # Split by @<TRIPOS>MOLECULE for multi-mol mol2
    text = path.read_text()
    blocks = text.split("@<TRIPOS>MOLECULE")
    blocks = [b for b in blocks if b.strip()]

    if len(blocks) <= 1:
        # Single molecule — parse directly
        result = _parse_one_mol2(str(path), path.name, entry_idx=0)
        if result:
            logger.info(f"    Parsed: 1 molecule, 0 failed")
            return [result]
        else:
            logger.info(f"    Parsed: 0 molecules, 1 failed")
            return []

    # Multiple molecules — split and parse each
    import tempfile

    molecules = []
    n_failed = 0

    for idx, block in enumerate(blocks):
        full_block = "@<TRIPOS>MOLECULE" + block
        try:
            with tempfile.NamedTemporaryFile(suffix=".mol2", mode="w",
                                             delete=False) as tmp:
                tmp.write(full_block)
                tmp_path = tmp.name

            result = _parse_one_mol2(tmp_path, path.name, entry_idx=idx)
            Path(tmp_path).unlink(missing_ok=True)

            if result:
                molecules.append(result)
            else:
                n_failed += 1
        except Exception as e:
            n_failed += 1
            logger.debug(f"    Error in mol2 block {idx}: {e}")

    logger.info(f"    Parsed: {len(molecules)} molecules, {n_failed} failed")
    return molecules


# =============================================================================
# SMILES TEXT FILE PARSING
# =============================================================================

def parse_smiles_file(smiles_path: str) -> List[Dict[str, Any]]:
    """
    Parse a SMILES text file.
    Format: one molecule per line, either "SMILES" or "SMILES name".
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit required for SMILES parsing")

    path = Path(smiles_path)
    logger.info(f"  Parsing SMILES file: {path.name}")

    rd_logger = RDLogger.logger()
    rd_logger.setLevel(RDLogger.ERROR)

    molecules = []
    n_failed = 0

    for idx, line in enumerate(path.read_text().strip().split("\n")):
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        parts = line.split()
        smiles_raw = parts[0]
        name = parts[1] if len(parts) > 1 else f"SMI_{idx:04d}"

        mol = Chem.MolFromSmiles(smiles_raw)
        if mol is None:
            n_failed += 1
            continue

        smiles = Chem.MolToSmiles(mol)
        molecules.append({
            "name": name,
            "mol": Chem.AddHs(mol),
            "mol_noH": mol,
            "smiles": smiles,
            "heavy_atoms": mol.GetNumHeavyAtoms(),
            "rmsd": None,
            "all_ids": name,
            "source_file": path.name,
            "source_format": "smiles",
            "entry_idx": idx,
            "extra_props": {},
        })

    logger.info(f"    Parsed: {len(molecules)} molecules, {n_failed} failed")
    return molecules


# =============================================================================
# UNIFIED PARSER (dispatches by format)
# =============================================================================

def parse_input(
        input_path: Union[str, Path],
        name_source: str = "header_first_word",
        name_property: str = "_Name",
        name_column: Optional[str] = None,
        smiles_column: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """
    Parse any supported input into a unified molecule list.

    Handles single files, directories, and mixed formats.
    """
    path = Path(input_path)
    fmt = detect_input_format(path)
    logger.info(f"Input: {path} (detected: {fmt})")

    if fmt == "sdf":
        return parse_sdf(str(path), name_source, name_property)

    elif fmt in ("csv", "xlsx"):
        return parse_tabular(str(path), name_column, smiles_column)

    elif fmt == "smiles":
        return parse_smiles_file(str(path))

    elif fmt == "mol2":
        return parse_single_mol2(str(path))

    elif fmt == "mol2_dir":
        return parse_mol2_directory(str(path))

    elif fmt == "directory":
        # Scan directory and parse all supported files
        file_groups = discover_files_in_directory(path)
        all_molecules = []

        for file_fmt, files in file_groups.items():
            for f in files:
                try:
                    if file_fmt == "sdf":
                        mols = parse_sdf(str(f), name_source, name_property)
                    elif file_fmt in ("csv", "xlsx"):
                        mols = parse_tabular(str(f), name_column, smiles_column)
                    elif file_fmt == "smiles":
                        mols = parse_smiles_file(str(f))
                    elif file_fmt == "mol2":
                        mols = parse_single_mol2(str(f))
                    else:
                        continue
                    all_molecules.extend(mols)
                except Exception as e:
                    logger.warning(f"  Failed to parse {f.name}: {e}")

        logger.info(f"  Directory total: {len(all_molecules)} molecules "
                     f"from {sum(len(v) for v in file_groups.values())} files")
        return all_molecules

    else:
        raise ValueError(f"Unknown format: {fmt}")


# =============================================================================
# CONFORMER RESOLUTION
# =============================================================================

def resolve_conformers(
        molecules: List[Dict[str, Any]],
        strategy: str = "best_rmsd",
) -> List[Dict[str, Any]]:
    """
    Resolve multiple conformers per molecule to a single representative.

    Groups by name, then picks one per group.

    Args:
        molecules: List from parse_input
        strategy: "best_rmsd" | "first" | "all"

    Returns:
        List of unique molecule dicts (with n_conformers added)
    """
    if strategy == "all":
        for m in molecules:
            m["n_conformers"] = 1
        logger.info(f"  Strategy 'all': keeping all {len(molecules)} entries")
        return molecules

    # Group by name
    groups: Dict[str, List[Dict]] = defaultdict(list)
    for mol_data in molecules:
        groups[mol_data["name"]].append(mol_data)

    unique = []
    for name, conformers in groups.items():
        if strategy == "best_rmsd":
            with_rmsd = [c for c in conformers if c["rmsd"] is not None]
            if with_rmsd:
                best = min(with_rmsd, key=lambda x: x["rmsd"])
            else:
                best = conformers[0]
        elif strategy == "first":
            best = conformers[0]
        else:
            raise ValueError(f"Unknown conformer strategy: {strategy}")

        best["n_conformers"] = len(conformers)
        unique.append(best)

    logger.info(f"  Strategy '{strategy}': {len(molecules)} entries -> "
                f"{len(unique)} unique molecules")
    return unique


# =============================================================================
# SMILES DUPLICATE DETECTION
# =============================================================================

def detect_smiles_duplicates(molecules: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Flag molecules that have different names but identical canonical SMILES.

    Adds 'smiles_duplicate_of' field (None if unique, name of first occurrence if duplicate).
    Does NOT remove duplicates — just flags them for the user to decide.
    """
    smiles_seen: Dict[str, str] = {}  # smiles -> first name
    n_dupes = 0

    for mol_data in molecules:
        smi = mol_data["smiles"]
        if smi in smiles_seen:
            mol_data["smiles_duplicate_of"] = smiles_seen[smi]
            n_dupes += 1
        else:
            smiles_seen[smi] = mol_data["name"]
            mol_data["smiles_duplicate_of"] = None

    if n_dupes > 0:
        logger.warning(f"  {n_dupes} SMILES duplicates detected (different names, same structure)")

    return molecules


# =============================================================================
# MOLECULAR PROPERTIES
# =============================================================================

def compute_molecular_properties(mol) -> Dict[str, Any]:
    """
    Compute standard molecular properties for dock2profile compatibility.

    Returns dict with: MW, LogP, HBD, HBA, TPSA, Rotatable_Bonds, QED,
    Formal_Charge, Heavy_Atoms, Rings, Aromatic_Rings
    """
    if not RDKIT_AVAILABLE or mol is None:
        return {}

    try:
        return {
            "MW": round(Descriptors.MolWt(mol), 2),
            "LogP": round(Descriptors.MolLogP(mol), 3),
            "HBD": Descriptors.NumHDonors(mol),
            "HBA": Descriptors.NumHAcceptors(mol),
            "TPSA": round(Descriptors.TPSA(mol), 2),
            "Rotatable_Bonds": Descriptors.NumRotatableBonds(mol),
            "QED": round(Descriptors.qed(mol), 4),
            "Formal_Charge": Chem.GetFormalCharge(mol),
            "Heavy_Atoms": mol.GetNumHeavyAtoms(),
            "Rings": Descriptors.RingCount(mol),
            "Aromatic_Rings": Descriptors.NumAromaticRings(mol),
        }
    except Exception as e:
        logger.warning(f"  Property calculation failed: {e}")
        return {}


# =============================================================================
# OUTPUT GENERATION
# =============================================================================

def save_outputs(
        molecules: List[Dict[str, Any]],
        output_dir: Path,
        compute_props: bool = True,
) -> Dict[str, str]:
    """
    Save parsed molecules to CSV, SDF, and summary.

    Returns:
        Dict with output file paths
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Build DataFrame ---
    rows = []
    for mol_data in molecules:
        row = {
            "Name": mol_data["name"],
            "SMILES": mol_data["smiles"],
            "Heavy_Atoms": mol_data["heavy_atoms"],
            "Phermit_RMSD": mol_data.get("rmsd"),
            "N_Conformers": mol_data.get("n_conformers", 1),
            "Source_File": mol_data.get("source_file", ""),
            "Source_Format": mol_data.get("source_format", ""),
            "All_IDs": mol_data.get("all_ids", ""),
            "SMILES_Duplicate_Of": mol_data.get("smiles_duplicate_of"),
            "Source_mol2_path": mol_data.get("source_mol2_path", ""),
            "Source_SDF_path": mol_data.get("source_sdf_path", ""),
        }

        if compute_props:
            props = compute_molecular_properties(mol_data.get("mol_noH"))
            row.update(props)

        rows.append(row)

    df = pd.DataFrame(rows)

    # --- Save CSV ---
    csv_path = output_dir / "unique_molecules.csv"
    df.to_csv(csv_path, index=False, encoding="utf-8")
    logger.info(f"  Saved CSV: {csv_path} ({len(df)} molecules)")

    # --- Save SDF ---
    sdf_path = output_dir / "unique_molecules.sdf"
    if RDKIT_AVAILABLE:
        writer = Chem.SDWriter(str(sdf_path))
        for mol_data in molecules:
            mol = mol_data.get("mol")
            if mol is None:
                continue
            mol.SetProp("_Name", mol_data["name"])
            mol.SetProp("SMILES", mol_data["smiles"])
            if mol_data.get("rmsd") is not None:
                mol.SetProp("rmsd", f"{mol_data['rmsd']:.6f}")
            if mol_data.get("source_file"):
                mol.SetProp("source_file", mol_data["source_file"])
            writer.write(mol)
        writer.close()
        logger.info(f"  Saved SDF: {sdf_path}")

    # --- Summary text ---
    summary_path = output_dir / "parse_summary.txt"
    w = 70
    lines = [
        "=" * w,
        "00a MOLECULE PARSER - SUMMARY",
        "=" * w,
        "",
        f"Date:                {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"Unique molecules:    {len(molecules)}",
        "",
    ]

    # Source breakdown
    source_counts = defaultdict(int)
    for m in molecules:
        source_counts[m.get("source_file", "unknown")] += 1
    if source_counts:
        lines.append("Sources:")
        for src, count in sorted(source_counts.items()):
            lines.append(f"  {src}: {count} molecules")
        lines.append("")

    # Duplicate warning
    n_dupes = sum(1 for m in molecules if m.get("smiles_duplicate_of"))
    if n_dupes:
        lines.append(f"WARNING: {n_dupes} SMILES duplicates (different name, same structure)")
        lines.append("")

    # Molecule table
    lines.extend([
        "-" * w,
        f"{'Name':<30} {'HA':>4} {'MW':>8} {'LogP':>6} {'RMSD':>7} {'Source':<15}",
        "-" * w,
    ])

    for mol_data in molecules:
        name = mol_data["name"]
        if len(name) > 28:
            name = name[:25] + "..."
        rmsd_str = f"{mol_data['rmsd']:.3f}" if mol_data.get("rmsd") is not None else "—"
        # Get MW/LogP from computed props
        props = compute_molecular_properties(mol_data.get("mol_noH")) if compute_props else {}
        mw_str = f"{props.get('MW', 0):.1f}" if props.get("MW") else "—"
        logp_str = f"{props.get('LogP', 0):.2f}" if props.get("LogP") is not None else "—"
        src = mol_data.get("source_format", "")

        dupe_marker = " *" if mol_data.get("smiles_duplicate_of") else ""
        lines.append(
            f"{name:<30} {mol_data['heavy_atoms']:>4} {mw_str:>8} "
            f"{logp_str:>6} {rmsd_str:>7} {src:<15}{dupe_marker}"
        )

    if n_dupes:
        lines.extend(["", "* = SMILES duplicate of another molecule"])

    lines.extend(["", "=" * w])

    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    logger.info(f"  Saved summary: {summary_path}")

    return {
        "csv": str(csv_path),
        "sdf": str(sdf_path),
        "summary": str(summary_path),
    }


# =============================================================================
# MAIN PIPELINE FUNCTION
# =============================================================================

def run_molecule_parser(
        input_file: Union[str, Path],
        output_dir: Union[str, Path],
        name_source: str = "header_first_word",
        name_property: str = "_Name",
        name_column: Optional[str] = None,
        smiles_column: Optional[str] = None,
        conformer_strategy: str = "best_rmsd",
        compute_properties: bool = True,
        detect_smiles_duplicates: bool = True,
) -> Dict[str, Any]:
    """
    Run the complete molecule parsing pipeline.

    Accepts any supported input format (SDF, CSV, XLSX, mol2, SMILES txt,
    or a directory containing a mix of these).

    Pipeline:
        1. Detect format and parse all molecules
        2. Resolve conformers (by name)
        3. Detect SMILES duplicates (optional)
        4. Compute molecular properties (optional)
        5. Save CSV + SDF + summary

    Args:
        input_file: Path to input file or directory
        output_dir: Directory for output files
        name_source: How to extract name from SDF headers
        name_property: SDF property for name (if name_source="property")
        name_column: Column name in CSV/XLSX for molecule name
        smiles_column: Column name in CSV/XLSX for SMILES
        conformer_strategy: "best_rmsd" | "first" | "all"
        compute_properties: Calculate MW, LogP, HBD, HBA, TPSA, QED
        detect_smiles_duplicates: Flag duplicate SMILES with different names

    Returns:
        Dict with summary statistics and output file paths
    """
    logger.info("Pipeline: parse -> resolve -> deduplicate -> save")

    # --- Step 1: Parse ---
    molecules = parse_input(
        input_path=input_file,
        name_source=name_source,
        name_property=name_property,
        name_column=name_column,
        smiles_column=smiles_column,
    )

    if not molecules:
        logger.error("No molecules parsed. Check input file.")
        return {"n_molecules": 0, "error": "No molecules parsed"}

    logger.info(f"Total parsed: {len(molecules)} entries")

    # --- Step 2: Resolve conformers ---
    unique = resolve_conformers(molecules, strategy=conformer_strategy)

    # --- Step 3: Detect SMILES duplicates ---
    if detect_smiles_duplicates:
        unique = detect_smiles_duplicates_fn(unique)

    # --- Step 4: Save outputs ---
    paths = save_outputs(unique, Path(output_dir), compute_props=compute_properties)

    # --- Summary stats ---
    n_dupes = sum(1 for m in unique if m.get("smiles_duplicate_of"))
    result = {
        "n_total_parsed": len(molecules),
        "n_unique": len(unique),
        "n_smiles_duplicates": n_dupes,
        "conformer_strategy": conformer_strategy,
        "output_csv": paths["csv"],
        "output_sdf": paths["sdf"],
        "output_summary": paths["summary"],
    }

    logger.info(f"Done: {result['n_unique']} unique molecules")
    return result


# Alias to avoid name collision with the boolean parameter
detect_smiles_duplicates_fn = detect_smiles_duplicates
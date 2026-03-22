"""
Parse & Fragment (05a) — Parse docking poses + identify rigid fragments
=========================================================================
Combines pose coordinate extraction with molecular graph fragmentation
into a single step. Each molecule gets:
  - 3D coordinates for all poses (heavy atoms only)
  - Score values per pose
  - Rigid fragment definitions (from rotatable bond analysis)
  - Atom-to-fragment mapping

This is the foundation for all downstream analysis (05b-05f).

Supported engines:
  DOCK6: multi-pose scored mol2 (##########... + @<TRIPOS>MOLECULE blocks)
  GNINA: multi-pose SDF ($$$$-delimited, with score properties)

Fragment identification:
  1. Load molecular graph from GNINA output SDF (preserves pose atom order)
  2. Detect rotatable bonds (RDKit SMARTS, excludes rings/amides/terminal)
  3. Cut graph at rotatable bonds → rigid segments
  4. Merge small fragments (< min_fragment_size heavy atoms) into neighbors
  5. Label as ring or linker

Location: 01_src/molecular_docking/m05_gnina_analysis/parse_and_fragment.py
Project: molecular_docking
Module: 05a (core)
Version: 2.0
"""

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.logger().setLevel(RDLogger.ERROR)
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logger.warning("RDKit not available — fragmentation requires RDKit")


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass
class FragmentDef:
    """Definition of a rigid fragment within a molecule."""
    fragment_id: int
    label: str                              # "ring_0", "linker_1", etc.
    heavy_atom_indices: List[int]           # indices into heavy_atom arrays
    elements: List[str]                     # element symbols for this fragment
    n_atoms: int
    is_ring: bool = False
    is_linker: bool = False


@dataclass
class MoleculeData:
    """All parsed data for one molecule: poses + fragments."""
    name: str
    engine: str
    n_poses: int
    n_heavy_atoms: int
    heavy_atom_elements: List[str]          # ["C", "N", "O", "P", ...]
    coords: np.ndarray                      # shape (n_poses, n_heavy_atoms, 3)
    pose_indices: np.ndarray                # shape (n_poses,)
    score_keys: List[str]                   # ["vina_affinity", "cnn_score", ...]
    score_values: np.ndarray                # shape (n_poses, n_score_keys)
    fragments: List[FragmentDef]
    source_file: str = ""
    # Per-atom interaction terms from GNINA --atom_term_data
    # shape (n_poses, n_heavy_atoms, n_terms) or None if not available
    atom_terms: Optional[np.ndarray] = None
    atom_term_keys: Optional[List[str]] = None  # ["gauss1", "gauss2", "repulsion", ...]


# =============================================================================
# SYBYL ATOM TYPE MAPPING (for DOCK6 mol2)
# =============================================================================

_SYBYL_TO_ELEMENT = {
    "C.3": "C", "C.2": "C", "C.1": "C", "C.ar": "C", "C.cat": "C",
    "N.3": "N", "N.2": "N", "N.1": "N", "N.ar": "N", "N.am": "N",
    "N.pl3": "N", "N.4": "N",
    "O.3": "O", "O.2": "O", "O.co2": "O", "O.spc": "O", "O.t3p": "O",
    "S.3": "S", "S.2": "S", "S.o": "S", "S.o2": "S",
    "P.3": "P",
    "F": "F", "Cl": "Cl", "Br": "Br", "I": "I",
    "H": "H", "H.spc": "H", "H.t3p": "H",
    "Li": "Li", "Na": "Na", "K": "K",
    "Ca": "Ca", "Mg": "Mg", "Zn": "Zn", "Fe": "Fe",
    "Si": "Si", "Se": "Se",
}


def _sybyl_to_element(atom_type: str) -> str:
    if atom_type in _SYBYL_TO_ELEMENT:
        return _SYBYL_TO_ELEMENT[atom_type]
    base = atom_type.split(".")[0]
    if base in _SYBYL_TO_ELEMENT:
        return _SYBYL_TO_ELEMENT[base]
    return base.capitalize() if base else "X"


# =============================================================================
# DOCK6 MOL2 PARSING
# =============================================================================

def _extract_heavy_atoms_mol2(mol2_lines: List[str]):
    """Extract heavy atom coords + elements from mol2 ATOM block."""
    coords, elements = [], []
    in_atom = False
    for line in mol2_lines:
        if "@<TRIPOS>ATOM" in line:
            in_atom = True
            continue
        if line.startswith("@<TRIPOS>") and in_atom:
            break
        if not in_atom or not line.strip():
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        element = _sybyl_to_element(parts[5])
        if element == "H":
            continue
        try:
            coords.append([float(parts[2]), float(parts[3]), float(parts[4])])
            elements.append(element)
        except (ValueError, IndexError):
            continue
    if not coords:
        return np.empty((0, 3)), []
    return np.array(coords, dtype=np.float64), elements


def parse_dock6_mol2(mol2_path: str, name: str = ""):
    """Parse DOCK6 multi-pose scored mol2. Returns list of (coords, elements, scores)."""
    path = Path(mol2_path)
    if not path.exists() or path.stat().st_size == 0:
        return []
    if not name:
        name = path.stem.replace("_scored", "")

    text = path.read_text()
    raw_poses = []
    current_header, current_mol2 = [], []
    in_mol2 = False

    for line in text.split("\n"):
        if line.startswith("##########"):
            if in_mol2 and current_mol2:
                raw_poses.append((current_header, current_mol2))
                current_header, current_mol2 = [], []
                in_mol2 = False
            current_header.append(line)
        elif "@<TRIPOS>MOLECULE" in line:
            in_mol2 = True
            current_mol2.append(line)
        elif in_mol2:
            current_mol2.append(line)

    if current_mol2:
        raw_poses.append((current_header, current_mol2))

    results = []
    for idx, (header, mol2_lines) in enumerate(raw_poses):
        scores = {}
        for line in header:
            if ":" in line:
                parts = line.split(":", 1)
                key = parts[0].replace("#", "").strip()
                try:
                    scores[key] = float(parts[1].strip())
                except ValueError:
                    pass
        coords, elements = _extract_heavy_atoms_mol2(mol2_lines)
        if len(coords) > 0:
            results.append({"idx": idx, "coords": coords, "elements": elements,
                            "scores": scores, "atom_terms": None})

    return results


# =============================================================================
# GNINA SDF PARSING
# =============================================================================

_GNINA_PROP_MAP = {
    "minimizedAffinity": "vina_affinity",
    "minimizedRMSD": "vina_rmsd",
    "CNNscore": "cnn_score",
    "CNNaffinity": "cnn_affinity",
    "CNN_VS": "cnn_vs",
    "CNNaffinity_variance": "cnn_affinity_variance",
}

_COMMON_ELEMENTS = {"C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "H",
                     "Si", "Se", "B", "Li", "Na", "K", "Ca", "Mg", "Zn", "Fe"}


def _extract_heavy_atoms_sdf(block: str):
    """Extract heavy atom coords + elements from SDF mol block."""
    lines = block.split("\n")
    counts_idx = None
    for i in range(min(10, len(lines))):
        parts = lines[i].split()
        if len(parts) >= 2:
            try:
                na, nb = int(parts[0]), int(parts[1])
                if 1 <= na <= 999 and 0 <= nb <= 9999:
                    counts_idx = i
                    break
            except ValueError:
                continue
    if counts_idx is None:
        return np.empty((0, 3)), []

    n_atoms = int(lines[counts_idx].split()[0])
    coords, elements = [], []
    for j in range(counts_idx + 1, min(counts_idx + 1 + n_atoms, len(lines))):
        parts = lines[j].split()
        if len(parts) < 4:
            continue
        try:
            x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
        except ValueError:
            continue
        elem = parts[3].strip()
        if elem == "H":
            continue
        if elem not in _COMMON_ELEMENTS:
            elem = elem.capitalize()
        coords.append([x, y, z])
        elements.append(elem)

    if not coords:
        return np.empty((0, 3)), []
    return np.array(coords, dtype=np.float64), elements


def _parse_sdf_properties(block: str) -> Dict[str, float]:
    """Extract > <property> values from SDF block."""
    props = {}
    lines = block.split("\n")
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("> <"):
            prop_name = line.split("> <")[1].rstrip(">").strip()
            if i + 1 < len(lines):
                try:
                    props[prop_name] = float(lines[i + 1].strip())
                except ValueError:
                    pass
            i += 2
        else:
            i += 1
    return props


def parse_gnina_sdf(sdf_path: str, name: str = ""):
    """Parse GNINA multi-pose SDF. Returns list of (coords, elements, scores)."""
    path = Path(sdf_path)
    if not path.exists() or path.stat().st_size == 0:
        return []
    if not name:
        name = path.stem.replace("_gnina", "")

    text = path.read_text()
    blocks = text.split("$$$$")
    results = []

    for idx, block in enumerate(blocks):
        block = block.strip()
        if not block:
            continue
        coords, elements = _extract_heavy_atoms_sdf(block)
        if len(coords) == 0:
            continue
        raw_props = _parse_sdf_properties(block)
        scores = {}
        for raw_key, canon_key in _GNINA_PROP_MAP.items():
            if raw_key in raw_props:
                scores[canon_key] = raw_props[raw_key]
        for k, v in raw_props.items():
            if k not in _GNINA_PROP_MAP:
                scores[k] = v
        # Parse atomic interaction terms if present
        atom_terms = _parse_atomic_interaction_terms(block, len(coords))

        results.append({"idx": idx, "coords": coords, "elements": elements,
                        "scores": scores, "atom_terms": atom_terms})

    return results


# =============================================================================
# ATOM INTERACTION TERMS PARSING (GNINA --atom_term_data)
# =============================================================================

# Vina scoring function weights (fixed, from GNINA/smina source)
VINA_WEIGHTS = {
    "gauss1":      -0.035579,
    "gauss2":      -0.005156,
    "repulsion":    0.840245,
    "hydrophobic": -0.035069,
    "hbond":       -0.587439,
}

ATOM_TERM_KEYS = ["gauss1", "gauss2", "repulsion", "hydrophobic", "hbond"]


def _parse_atomic_interaction_terms(block: str, n_heavy_atoms: int):
    """
    Parse atomic_interaction_terms from an SDF property block.

    Format (from GNINA --atom_term_data):
        > <atomic_interaction_terms>
        atomid el pos gauss(...) gauss(...) repulsion(...) hydrophobic(...) non_dir_h_bond(...)
         <x,y,z> val1 val2 val3 val4 val5
         <x,y,z> val1 val2 val3 val4 val5
         ...

    Returns:
        np.ndarray of shape (n_heavy_atoms, 5) with term values for heavy atoms only,
        or None if not found/parseable.

    Note: GNINA outputs terms for ALL atoms (including H). We parse all,
    then return only the first n_heavy_atoms entries (heavy atoms come first
    in the SDF atom block since we filter H during coord extraction).
    """
    lines = block.split("\n")

    # Find the property
    start_idx = None
    for i, line in enumerate(lines):
        if "> <atomic_interaction_terms>" in line:
            start_idx = i + 1
            break

    if start_idx is None:
        return None

    # Skip header line (atomid el pos gauss...)
    if start_idx < len(lines) and "gauss" in lines[start_idx]:
        start_idx += 1

    # Parse atom lines
    all_terms = []
    for i in range(start_idx, len(lines)):
        line = lines[i].strip()
        if not line or line.startswith(">") or line.startswith("$$$$"):
            break

        # Format: <x,y,z> val1 val2 val3 val4 val5
        # or:     atomid el <x,y,z> val1 val2 val3 val4 val5
        parts = line.split()

        # Find the position token (starts with <)
        pos_idx = None
        for j, p in enumerate(parts):
            if p.startswith("<"):
                pos_idx = j
                break

        if pos_idx is None:
            continue

        # Values follow after position
        val_parts = parts[pos_idx + 1:]
        if len(val_parts) < 5:
            continue

        try:
            terms = [float(v) for v in val_parts[:5]]
            all_terms.append(terms)
        except ValueError:
            continue

    if not all_terms:
        return None

    all_arr = np.array(all_terms, dtype=np.float64)

    # Return only heavy atom rows (first n_heavy_atoms)
    # This works because in the SDF, heavy atoms come before H atoms
    # in the atom block, and GNINA preserves this order
    if len(all_arr) >= n_heavy_atoms:
        return all_arr[:n_heavy_atoms]
    else:
        # Fewer terms than heavy atoms — might be heavy-only output
        return all_arr


# =============================================================================
# FRAGMENT IDENTIFICATION (RDKit)
# =============================================================================

def _get_heavy_atom_map(mol) -> List[int]:
    """RDKit atom indices for heavy atoms only."""
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]


def _get_rotatable_bonds(mol) -> List[Tuple[int, int]]:
    """Rotatable bond atom pairs (excludes rings, terminal, amides)."""
    rot_smarts = Chem.MolFromSmarts(
        "[!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]"
    )
    return list(mol.GetSubstructMatches(rot_smarts))


def _build_rigid_segments(mol, rot_bonds: List[Tuple[int, int]]) -> List[Set[int]]:
    """Connected components after removing rotatable bonds."""
    rot_set = {(min(a, b), max(a, b)) for a, b in rot_bonds}
    adj = {i: set() for i in range(mol.GetNumAtoms())}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if (min(a, b), max(a, b)) not in rot_set:
            adj[a].add(b)
            adj[b].add(a)

    visited, segments = set(), []
    for start in range(mol.GetNumAtoms()):
        if start in visited:
            continue
        comp = set()
        queue = [start]
        while queue:
            node = queue.pop(0)
            if node in visited:
                continue
            visited.add(node)
            comp.add(node)
            for nb in adj[node]:
                if nb not in visited:
                    queue.append(nb)
        segments.append(comp)
    return segments


def _merge_small(segments, mol, rot_bonds, min_heavy=4):
    """Merge segments with < min_heavy heavy atoms into largest neighbor."""
    def hcount(seg):
        return sum(1 for i in seg if mol.GetAtomWithIdx(i).GetAtomicNum() > 1)

    seg_for_atom = {}
    for i, seg in enumerate(segments):
        for a in seg:
            seg_for_atom[a] = i

    merged = [set(s) for s in segments]
    changed = True
    while changed:
        changed = False
        skip = set()
        seg_for_atom_local = {}
        for j, s in enumerate(merged):
            if j not in skip:
                for a in s:
                    seg_for_atom_local[a] = j

        for i, seg in enumerate(merged):
            if i in skip:
                continue
            if hcount(seg) < min_heavy:
                best_nb, best_sz = None, -1
                for a, b in rot_bonds:
                    sa, sb = seg_for_atom_local.get(a), seg_for_atom_local.get(b)
                    if sa == i and sb is not None and sb != i and sb not in skip:
                        if hcount(merged[sb]) > best_sz:
                            best_sz, best_nb = hcount(merged[sb]), sb
                    elif sb == i and sa is not None and sa != i and sa not in skip:
                        if hcount(merged[sa]) > best_sz:
                            best_sz, best_nb = hcount(merged[sa]), sa
                if best_nb is not None:
                    merged[best_nb] = merged[best_nb] | seg
                    skip.add(i)
                    changed = True
        if changed:
            merged = [s for j, s in enumerate(merged) if j not in skip]
    return merged


def identify_fragments(mol, min_heavy_atoms: int = 4) -> List[FragmentDef]:
    """
    Identify rigid fragments by cutting at rotatable bonds.
    Returns list sorted by size descending (largest = anchor).
    """
    # Find rotatable bonds (with robust ring info handling)
    try:
        try:
            Chem.FastFindRings(mol)
        except Exception:
            pass
        rot_bonds = _get_rotatable_bonds(mol)
    except Exception as e:
        logger.debug(f"  Rotatable bond detection failed: {e}")
        heavy_map = _get_heavy_atom_map(mol)
        if not heavy_map:
            return []
        elems = [mol.GetAtomWithIdx(i).GetSymbol() for i in heavy_map]
        return [FragmentDef(0, "rigid_0", list(range(len(heavy_map))),
                            elems, len(heavy_map), False, False)]

    if not rot_bonds:
        heavy_map = _get_heavy_atom_map(mol)
        elems = [mol.GetAtomWithIdx(i).GetSymbol() for i in heavy_map]
        return [FragmentDef(0, "rigid_0", list(range(len(heavy_map))),
                            elems, len(heavy_map), True, False)]

    segments = _build_rigid_segments(mol, rot_bonds)
    segments = _merge_small(segments, mol, rot_bonds, min_heavy_atoms)

    heavy_map = _get_heavy_atom_map(mol)
    rdkit_to_heavy = {ri: hi for hi, ri in enumerate(heavy_map)}

    # Ring atoms for labeling
    ring_atoms = set()
    try:
        for ring in mol.GetRingInfo().AtomRings():
            ring_atoms.update(ring)
    except Exception:
        pass

    fragments = []
    for si, seg in enumerate(segments):
        h_idx = sorted(rdkit_to_heavy[a] for a in seg if a in rdkit_to_heavy)
        if not h_idx:
            continue
        elems = [mol.GetAtomWithIdx(heavy_map[i]).GetSymbol() for i in h_idx]
        n_ring = sum(1 for a in seg if a in ring_atoms and a in rdkit_to_heavy)
        is_ring = n_ring > len(h_idx) * 0.5
        fragments.append(FragmentDef(
            si, f"{'ring' if is_ring else 'linker'}_{si}",
            h_idx, elems, len(h_idx), is_ring, not is_ring,
        ))

    fragments.sort(key=lambda f: f.n_atoms, reverse=True)
    for i, f in enumerate(fragments):
        f.fragment_id = i
        f.label = f"{'ring' if f.is_ring else 'linker'}_{i}"

    return fragments


# =============================================================================
# MOLECULE GRAPH LOADING (RDKit)
# =============================================================================

def load_mol_graph(sdf_path: Optional[str] = None,
                    mol2_path: Optional[str] = None):
    """
    Load molecular graph for fragmentation.
    Priority: SDF (matches GNINA atom order) > mol2.
    Handles sanitization failures gracefully.
    """
    if not RDKIT_AVAILABLE:
        return None

    for path, loader in [(sdf_path, "sdf"), (mol2_path, "mol2")]:
        if not path or not Path(path).exists():
            continue
        try:
            if loader == "sdf":
                suppl = Chem.SDMolSupplier(path, removeHs=False, sanitize=False)
                mol = next(iter(suppl), None)
            else:
                mol = Chem.MolFromMol2File(path, removeHs=False, sanitize=False)

            if mol is not None:
                try:
                    Chem.SanitizeMol(mol)
                except Exception:
                    try:
                        Chem.FastFindRings(mol)
                    except Exception:
                        pass
                return mol
        except Exception:
            continue
    return None


def verify_elements(mol, pose_elements: List[str]) -> bool:
    """Check RDKit heavy atom elements match pose elements (handles np.str_)."""
    heavy_map = _get_heavy_atom_map(mol)
    rdkit_elems = [mol.GetAtomWithIdx(i).GetSymbol() for i in heavy_map]
    pose_elems = [str(e) for e in pose_elements]
    return len(rdkit_elems) == len(pose_elems) and rdkit_elems == pose_elems


# =============================================================================
# MOLECULE SOURCE RESOLUTION
# =============================================================================

def resolve_sources(
        molecules_csv: Optional[str],
        antechamber_dir: Optional[str] = None,
        ionization_sdf_dir: Optional[str] = None,
        gnina_poses_dir: Optional[str] = None,
) -> Dict[str, Dict[str, str]]:
    """
    Build {name: {"sdf": path, "mol2": path}} for each molecule.
    GNINA output SDF overrides all others (guarantees atom order match).
    """
    sources = {}

    # From CSV metadata
    if molecules_csv and Path(molecules_csv).exists():
        df = pd.read_csv(molecules_csv)
        ncol = "Name" if "Name" in df.columns else "name"
        for _, row in df.iterrows():
            name = row.get(ncol, "")
            if not name:
                continue
            src = {}
            if "Source_mol2_path" in df.columns and pd.notna(row.get("Source_mol2_path")):
                src["mol2"] = str(row["Source_mol2_path"])
            if "Source_SDF_path" in df.columns and pd.notna(row.get("Source_SDF_path")):
                src["sdf"] = str(row["Source_SDF_path"])
            sources[name] = src

    # Antechamber mol2
    if antechamber_dir and Path(antechamber_dir).exists():
        for f in Path(antechamber_dir).glob("*.mol2"):
            sources.setdefault(f.stem, {}).setdefault("mol2", str(f))

    # Ionization SDF
    if ionization_sdf_dir and Path(ionization_sdf_dir).exists():
        for f in Path(ionization_sdf_dir).glob("*.sdf"):
            sources.setdefault(f.stem, {}).setdefault("sdf", str(f))

    # GNINA output SDF — FORCE OVERRIDE (matches pose atom order)
    if gnina_poses_dir and Path(gnina_poses_dir).exists():
        for f in Path(gnina_poses_dir).glob("*_gnina.sdf"):
            name = f.stem.replace("_gnina", "")
            sources.setdefault(name, {})["sdf"] = str(f)

    return sources


# =============================================================================
# SERIALIZATION
# =============================================================================

def save_molecule_data(mol_data: MoleculeData, output_dir: Path):
    """Save one molecule: .npz (coords/scores) + .json (fragments)."""
    # NPZ: coords, elements, scores, atom_terms
    save_dict = {
        "coords": mol_data.coords,
        "elements": np.array(mol_data.heavy_atom_elements, dtype="U4"),
        "pose_indices": mol_data.pose_indices,
        "score_keys": np.array(mol_data.score_keys, dtype="U40"),
        "score_values": mol_data.score_values,
    }
    if mol_data.atom_terms is not None:
        save_dict["atom_terms"] = mol_data.atom_terms
        save_dict["atom_term_keys"] = np.array(mol_data.atom_term_keys, dtype="U20")
    np.savez_compressed(output_dir / f"{mol_data.name}.npz", **save_dict)

    # JSON: fragment definitions
    frag_json = {
        "name": mol_data.name,
        "engine": mol_data.engine,
        "n_poses": mol_data.n_poses,
        "n_heavy_atoms": mol_data.n_heavy_atoms,
        "n_fragments": len(mol_data.fragments),
        "fragments": [
            {
                "fragment_id": f.fragment_id,
                "label": f.label,
                "heavy_atom_indices": f.heavy_atom_indices,
                "elements": f.elements,
                "n_atoms": f.n_atoms,
                "is_ring": f.is_ring,
                "is_linker": f.is_linker,
            }
            for f in mol_data.fragments
        ],
    }
    with open(output_dir / f"{mol_data.name}.json", "w") as fh:
        json.dump(frag_json, fh, indent=2)


def load_molecule_data(name: str, data_dir: str) -> Optional[MoleculeData]:
    """Load one molecule from .npz + .json saved by save_molecule_data."""
    data_path = Path(data_dir)
    npz_file = data_path / f"{name}.npz"
    json_file = data_path / f"{name}.json"

    if not npz_file.exists() or not json_file.exists():
        return None

    npz = np.load(npz_file, allow_pickle=False)
    with open(json_file) as fh:
        frag_json = json.load(fh)

    fragments = [
        FragmentDef(
            fragment_id=fd["fragment_id"],
            label=fd["label"],
            heavy_atom_indices=fd["heavy_atom_indices"],
            elements=fd["elements"],
            n_atoms=fd["n_atoms"],
            is_ring=fd["is_ring"],
            is_linker=fd["is_linker"],
        )
        for fd in frag_json["fragments"]
    ]

    # Atom terms (optional — only if GNINA ran with --atom_term_data)
    atom_terms = npz["atom_terms"] if "atom_terms" in npz else None
    atom_term_keys = list(npz["atom_term_keys"]) if "atom_term_keys" in npz else None

    return MoleculeData(
        name=frag_json["name"],
        engine=frag_json["engine"],
        n_poses=frag_json["n_poses"],
        n_heavy_atoms=frag_json["n_heavy_atoms"],
        heavy_atom_elements=list(npz["elements"]),
        coords=npz["coords"],
        pose_indices=npz["pose_indices"],
        score_keys=list(npz["score_keys"]),
        score_values=npz["score_values"],
        fragments=fragments,
        atom_terms=atom_terms,
        atom_term_keys=atom_term_keys,
    )


def load_all_molecules(data_dir: str,
                        molecule_names: Optional[List[str]] = None,
                        ) -> Dict[str, MoleculeData]:
    """Load all molecules from 05a output directory."""
    data_path = Path(data_dir)
    if not data_path.exists():
        return {}

    results = {}
    for json_file in sorted(data_path.glob("*.json")):
        name = json_file.stem
        if molecule_names and name not in molecule_names:
            continue
        mol_data = load_molecule_data(name, data_dir)
        if mol_data:
            results[name] = mol_data

    logger.info(f"  Loaded {len(results)} molecules from {data_dir}")
    return results


# =============================================================================
# PIPELINE ENTRY POINT
# =============================================================================

def run_parse_and_fragment(
        output_dir: str,
        engine: str = "gnina",
        dock6_dir: Optional[str] = None,
        gnina_poses_dir: Optional[str] = None,
        molecules_csv: Optional[str] = None,
        antechamber_dir: Optional[str] = None,
        ionization_sdf_dir: Optional[str] = None,
        min_fragment_size: int = 4,
        molecule_names: Optional[List[str]] = None,
) -> Dict[str, any]:
    """
    Parse poses + identify fragments for all molecules.

    For each molecule:
      1. Parse multi-pose file (mol2 or SDF) → coords + scores
      2. Load molecular graph from SDF → fragment identification
      3. Verify atom order match between graph and poses
      4. Save .npz (coords/scores) + .json (fragments)
    """
    logger.info("=" * 60)
    logger.info("05a PARSE & FRAGMENT v2.0")
    logger.info("=" * 60)

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # Resolve molecule sources for fragmentation
    sources = resolve_sources(
        molecules_csv, antechamber_dir, ionization_sdf_dir, gnina_poses_dir
    )
    logger.info(f"  Molecule sources: {len(sources)}")

    # Find pose files
    pose_files = {}
    if engine == "dock6" and dock6_dir:
        dock6_path = Path(dock6_dir)
        if dock6_path.exists():
            for d in sorted(dock6_path.iterdir()):
                if d.is_dir():
                    scored = d / f"{d.name}_scored.mol2"
                    if scored.exists():
                        pose_files[d.name] = str(scored)
    elif engine == "gnina" and gnina_poses_dir:
        gn_path = Path(gnina_poses_dir)
        if gn_path.exists():
            for f in sorted(gn_path.glob("*_gnina.sdf")):
                name = f.stem.replace("_gnina", "")
                pose_files[name] = str(f)

    if molecule_names:
        pose_files = {k: v for k, v in pose_files.items() if k in molecule_names}

    logger.info(f"  Pose files: {len(pose_files)} ({engine})")

    # Process each molecule
    n_ok, n_skip, n_no_frag = 0, 0, 0
    summary_rows = []

    for name in sorted(pose_files.keys()):
        pose_path = pose_files[name]

        # 1. Parse poses
        if engine == "dock6":
            raw_poses = parse_dock6_mol2(pose_path, name)
        else:
            raw_poses = parse_gnina_sdf(pose_path, name)

        if not raw_poses:
            n_skip += 1
            continue

        n_poses = len(raw_poses)
        elements = raw_poses[0]["elements"]
        n_atoms = len(elements)

        # Build arrays
        coords = np.zeros((n_poses, n_atoms, 3), dtype=np.float64)
        pose_indices = np.zeros(n_poses, dtype=np.int32)
        score_keys = sorted(raw_poses[0]["scores"].keys())
        score_values = np.full((n_poses, len(score_keys)), np.nan, dtype=np.float64)

        # Atom interaction terms (GNINA --atom_term_data)
        has_atom_terms = raw_poses[0].get("atom_terms") is not None
        if has_atom_terms:
            n_terms = 5  # gauss1, gauss2, repulsion, hydrophobic, hbond
            atom_terms_arr = np.zeros((n_poses, n_atoms, n_terms), dtype=np.float64)
        else:
            atom_terms_arr = None

        for i, rp in enumerate(raw_poses):
            if len(rp["coords"]) == n_atoms:
                coords[i] = rp["coords"]
            pose_indices[i] = rp["idx"]
            for j, key in enumerate(score_keys):
                if key in rp["scores"]:
                    score_values[i, j] = rp["scores"][key]
            # Atom terms
            if has_atom_terms and rp.get("atom_terms") is not None:
                at = rp["atom_terms"]
                if len(at) == n_atoms:
                    atom_terms_arr[i] = at
                elif len(at) > n_atoms:
                    atom_terms_arr[i] = at[:n_atoms]

        # 2. Load molecular graph for fragmentation
        src = sources.get(name, {})
        mol = load_mol_graph(sdf_path=src.get("sdf"), mol2_path=src.get("mol2"))

        fragments = []
        if mol is not None:
            # Verify atom order
            if verify_elements(mol, elements):
                try:
                    fragments = identify_fragments(mol, min_fragment_size)
                except Exception as e:
                    logger.debug(f"  {name}: fragmentation failed: {e}")
            else:
                logger.debug(f"  {name}: atom order mismatch, skipping fragmentation")

        if not fragments:
            # Fallback: whole molecule as one fragment
            fragments = [FragmentDef(
                0, "whole_0", list(range(n_atoms)), elements, n_atoms, False, False
            )]
            n_no_frag += 1

        # 3. Save
        mol_data = MoleculeData(
            name=name, engine=engine,
            n_poses=n_poses, n_heavy_atoms=n_atoms,
            heavy_atom_elements=elements,
            coords=coords, pose_indices=pose_indices,
            score_keys=score_keys, score_values=score_values,
            fragments=fragments, source_file=pose_path,
            atom_terms=atom_terms_arr,
            atom_term_keys=ATOM_TERM_KEYS if has_atom_terms else None,
        )
        save_molecule_data(mol_data, out_path)
        n_ok += 1

        summary_rows.append({
            "Name": name,
            "engine": engine,
            "n_poses": n_poses,
            "n_heavy_atoms": n_atoms,
            "n_fragments": len(fragments),
            "fragment_sizes": [f.n_atoms for f in fragments],
            "fragment_labels": [f.label for f in fragments],
            "has_atom_terms": has_atom_terms,
        })

    # Summary CSV
    df = pd.DataFrame(summary_rows)
    if not df.empty:
        csv_path = out_path / "parse_summary.csv"
        df.to_csv(csv_path, index=False)
        logger.info(f"  Summary: {csv_path}")

    logger.info(f"\n{'=' * 60}")
    logger.info("PARSE & FRAGMENT COMPLETE")
    logger.info(f"{'=' * 60}")
    logger.info(f"  Parsed: {n_ok} molecules")
    logger.info(f"  Skipped: {n_skip} (no poses)")
    logger.info(f"  No fragmentation (whole mol): {n_no_frag}")
    if summary_rows:
        n_frags = [r["n_fragments"] for r in summary_rows]
        logger.info(f"  Fragments/molecule: mean={np.mean(n_frags):.1f}, "
                     f"range={min(n_frags)}-{max(n_frags)}")

    return {
        "success": True,
        "n_molecules": n_ok,
        "n_skipped": n_skip,
        "n_no_fragmentation": n_no_frag,
        "output_dir": str(out_path),
    }

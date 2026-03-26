"""
Contact Mapping (05f) — Receptor residue contacts per fragment cluster
=========================================================================
For each fragment in each cluster position, identifies which receptor
residues are within contact distance. This connects the geometric
clustering (05b) with binding site biology.

Algorithm:
    1. Load receptor PDB → extract residue atoms with coordinates
    2. For each molecule's fragment cluster medoid pose:
       - Extract fragment atom coordinates
       - Find receptor residues with any atom < cutoff distance
       - Record: residue name, chain, number, min distance, contact type
    3. Aggregate across molecules per hotspot (from 05d):
       - Which residues are contacted by most molecules in a hotspot?
       - Consensus contacts = pharmacophoric interactions

Contact types (simple distance-based):
    < 2.5 Å: strong H-bond / salt bridge
    < 3.5 Å: H-bond / polar contact
    < 4.5 Å: van der Waals / hydrophobic
    < 6.0 Å: weak / solvent-mediated

Reads:  05a + 05b + 05d (hotspots) + receptor PDB
Saves:  per-molecule contact JSON, contact_summary.csv, hotspot_contacts.csv

Location: 01_src/molecular_docking/m05_gnina_analysis/contact_mapping.py
Project: molecular_docking
Module: 05f (core)
Version: 3.0 — atomic output + PLIP for top N (2026-03-26)
"""

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .parse_and_fragment import MoleculeData, load_all_molecules
from .fragment_clustering import MoleculeClusterResult, load_cluster_result

try:
    from plip.structure.preparation import PDBComplex
    HAS_PLIP = True
except ImportError:
    HAS_PLIP = False

logger = logging.getLogger(__name__)


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass
class ReceptorAtom:
    """One heavy atom from the receptor PDB."""
    index: int
    atom_name: str          # " CA ", " OD1", etc.
    residue_name: str       # "ASP", "ARG", etc.
    chain: str
    residue_number: int
    coords: np.ndarray      # (3,)
    element: str


@dataclass
class ResidueContact:
    """Contact between a fragment and a receptor residue."""
    residue_name: str
    chain: str
    residue_number: int
    residue_id: str         # "ASP73.A"
    min_distance: float     # closest atom-atom distance
    n_contacts: int         # number of atom pairs within cutoff
    contact_type: str       # "strong_hbond", "hbond", "vdw", "weak"
    closest_receptor_atom: str
    closest_ligand_atom: str


@dataclass
class FragmentContactResult:
    """All contacts for one fragment in one cluster position."""
    molecule_name: str
    fragment_id: int
    fragment_label: str
    cluster_id: int
    medoid_pose_index: int
    contacts: List[ResidueContact]


# =============================================================================
# RECEPTOR PDB PARSING
# =============================================================================

_ELEMENT_FROM_NAME = {
    "C": "C", "N": "N", "O": "O", "S": "S", "P": "P",
    "FE": "Fe", "ZN": "Zn", "MG": "Mg", "CA": "Ca", "MN": "Mn",
    "CL": "Cl", "BR": "Br", "SE": "Se",
}


def parse_receptor_pdb(pdb_path: str) -> List[ReceptorAtom]:
    """
    Parse receptor PDB, extract heavy atoms with coordinates.
    Skips HOH, hydrogens, and HETATM (ligands).
    """
    atoms = []
    path = Path(pdb_path)
    if not path.exists():
        logger.error(f"Receptor PDB not found: {pdb_path}")
        return []

    idx = 0
    for line in path.read_text().split("\n"):
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue

        # Skip water
        res_name = line[17:20].strip()
        if res_name in ("HOH", "WAT", "TIP", "SOL"):
            continue

        # Skip HETATM unless it's a cofactor/ion we want
        if line.startswith("HETATM"):
            # Keep common cofactors and ions
            if res_name not in ("ZN", "MG", "CA", "MN", "FE", "CU",
                                "NAG", "NAD", "FAD", "ATP", "ADP",
                                "GDP", "GTP", "UDP"):
                continue

        atom_name = line[12:16].strip()

        # Element from columns 76-78 or from atom name
        element = line[76:78].strip() if len(line) > 77 else ""
        if not element:
            element = atom_name[0] if atom_name else "X"
        element = element.strip().capitalize()

        # Skip hydrogen
        if element == "H" or atom_name.startswith("H") or atom_name.startswith("1H"):
            continue

        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except (ValueError, IndexError):
            continue

        chain = line[21:22].strip() or "A"
        try:
            res_num = int(line[22:26].strip())
        except ValueError:
            continue

        atoms.append(ReceptorAtom(
            index=idx, atom_name=atom_name,
            residue_name=res_name, chain=chain,
            residue_number=res_num,
            coords=np.array([x, y, z], dtype=np.float64),
            element=element,
        ))
        idx += 1

    logger.info(f"  Receptor: {len(atoms)} heavy atoms from {pdb_path}")
    return atoms


def _build_receptor_arrays(receptor_atoms: List[ReceptorAtom]):
    """Pre-compute receptor coordinate array for fast distance calc."""
    coords = np.array([a.coords for a in receptor_atoms], dtype=np.float64)
    return coords


# =============================================================================
# CONTACT DETECTION
# =============================================================================

def _classify_contact(distance: float) -> str:
    """Classify contact by distance."""
    if distance < 2.5:
        return "strong_hbond"
    elif distance < 3.5:
        return "hbond"
    elif distance < 4.5:
        return "vdw"
    else:
        return "weak"


def find_fragment_contacts(
        frag_coords: np.ndarray,        # (n_frag_atoms, 3)
        frag_elements: List[str],
        receptor_atoms: List[ReceptorAtom],
        receptor_coords: np.ndarray,     # (n_rec_atoms, 3)
        cutoff: float = 4.5,
) -> List[ResidueContact]:
    """
    Find receptor residues within cutoff of fragment atoms.

    Uses vectorized distance computation for speed.
    """
    n_frag = frag_coords.shape[0]
    n_rec = receptor_coords.shape[0]

    if n_frag == 0 or n_rec == 0:
        return []

    # Compute all pairwise distances: (n_frag, n_rec)
    # Using broadcasting: diff[i,j] = frag[i] - rec[j]
    diff = frag_coords[:, np.newaxis, :] - receptor_coords[np.newaxis, :, :]
    dist_matrix = np.sqrt(np.sum(diff ** 2, axis=2))  # (n_frag, n_rec)

    # Find residue-level contacts
    residue_contacts = {}  # residue_id -> best contact info

    # Get all pairs within cutoff
    frag_idx, rec_idx = np.where(dist_matrix < cutoff)

    for fi, ri in zip(frag_idx, rec_idx):
        d = dist_matrix[fi, ri]
        ra = receptor_atoms[ri]
        res_id = f"{ra.residue_name}{ra.residue_number}.{ra.chain}"

        if res_id not in residue_contacts:
            residue_contacts[res_id] = {
                "residue_name": ra.residue_name,
                "chain": ra.chain,
                "residue_number": ra.residue_number,
                "min_distance": d,
                "n_contacts": 1,
                "closest_receptor_atom": ra.atom_name,
                "closest_ligand_atom": frag_elements[fi] if fi < len(frag_elements) else "X",
            }
        else:
            rc = residue_contacts[res_id]
            rc["n_contacts"] += 1
            if d < rc["min_distance"]:
                rc["min_distance"] = d
                rc["closest_receptor_atom"] = ra.atom_name
                rc["closest_ligand_atom"] = frag_elements[fi] if fi < len(frag_elements) else "X"

    # Convert to ResidueContact objects
    result = []
    for res_id, info in sorted(residue_contacts.items(), key=lambda x: x[1]["min_distance"]):
        result.append(ResidueContact(
            residue_name=info["residue_name"],
            chain=info["chain"],
            residue_number=info["residue_number"],
            residue_id=res_id,
            min_distance=round(info["min_distance"], 2),
            n_contacts=info["n_contacts"],
            contact_type=_classify_contact(info["min_distance"]),
            closest_receptor_atom=info["closest_receptor_atom"],
            closest_ligand_atom=info["closest_ligand_atom"],
        ))

    return result


# =============================================================================
# MEDOID FINDER (reused from 05e)
# =============================================================================

def _find_medoid(rmsd_matrix: np.ndarray, members: np.ndarray) -> int:
    if len(members) <= 1:
        return int(members[0])
    sub = rmsd_matrix[np.ix_(members, members)]
    return int(members[int(np.argmin(np.mean(sub, axis=1)))])


# =============================================================================
# PER-MOLECULE CONTACT ANALYSIS
# =============================================================================

def analyze_molecule_contacts(
        mol_data: MoleculeData,
        cluster_result: MoleculeClusterResult,
        receptor_atoms: List[ReceptorAtom],
        receptor_coords: np.ndarray,
        cutoff: float = 4.5,
        max_clusters_per_fragment: int = 3,
) -> List[FragmentContactResult]:
    """Compute contacts for all fragment clusters of one molecule."""
    results = []

    for fr in cluster_result.fragment_results:
        frag = mol_data.fragments[fr.fragment_id]

        unique_cls = sorted(set(fr.labels) - {-1},
                           key=lambda c: -int(np.sum(fr.labels == c)))

        for cl in unique_cls[:max_clusters_per_fragment]:
            members = np.where(fr.labels == cl)[0]
            if len(members) == 0:
                continue

            medoid_idx = _find_medoid(fr.rmsd_matrix, members)
            frag_coords = mol_data.coords[medoid_idx, frag.heavy_atom_indices, :]
            frag_elements = [mol_data.heavy_atom_elements[i] for i in frag.heavy_atom_indices]

            contacts = find_fragment_contacts(
                frag_coords, frag_elements,
                receptor_atoms, receptor_coords,
                cutoff=cutoff,
            )

            results.append(FragmentContactResult(
                molecule_name=mol_data.name,
                fragment_id=fr.fragment_id,
                fragment_label=fr.label,
                cluster_id=cl,
                medoid_pose_index=medoid_idx,
                contacts=contacts,
            ))

    return results


# =============================================================================
# HOTSPOT CONTACT AGGREGATION
# =============================================================================

def aggregate_hotspot_contacts(
        all_contacts: Dict[str, List[FragmentContactResult]],
        hotspot_map: Optional[Dict],
) -> pd.DataFrame:
    """
    For each hotspot, aggregate which residues are contacted by which molecules.
    Returns a DataFrame of consensus contacts.
    """
    if not hotspot_map:
        return pd.DataFrame()

    rows = []
    for hs in hotspot_map.get("hotspots", []):
        hs_id = hs["hotspot_id"]

        # Collect all contacts from hotspot members
        residue_votes = {}  # residue_id -> {molecules: set, distances: [], types: []}

        for member in hs.get("members", []):
            mol_name = member["molecule"]
            frag_id = member["fragment_id"]
            cl_id = member["cluster_id"]

            mol_contacts = all_contacts.get(mol_name, [])
            for fc in mol_contacts:
                if fc.fragment_id == frag_id and fc.cluster_id == cl_id:
                    for c in fc.contacts:
                        if c.residue_id not in residue_votes:
                            residue_votes[c.residue_id] = {
                                "residue_name": c.residue_name,
                                "chain": c.chain,
                                "residue_number": c.residue_number,
                                "molecules": set(),
                                "distances": [],
                                "types": [],
                            }
                        rv = residue_votes[c.residue_id]
                        rv["molecules"].add(mol_name)
                        rv["distances"].append(c.min_distance)
                        rv["types"].append(c.contact_type)

        for res_id, rv in sorted(residue_votes.items(),
                                  key=lambda x: -len(x[1]["molecules"])):
            n_mols = len(rv["molecules"])
            rows.append({
                "hotspot_id": hs_id,
                "hotspot_n_molecules": hs["n_molecules"],
                "residue_id": res_id,
                "residue_name": rv["residue_name"],
                "chain": rv["chain"],
                "residue_number": rv["residue_number"],
                "n_molecules_contacting": n_mols,
                "contact_fraction": round(n_mols / hs["n_molecules"], 3),
                "mean_distance": round(float(np.mean(rv["distances"])), 2),
                "min_distance": round(float(np.min(rv["distances"])), 2),
                "dominant_contact_type": max(set(rv["types"]), key=rv["types"].count),
            })

    return pd.DataFrame(rows)


# =============================================================================
# SERIALIZATION
# =============================================================================

def save_contact_results(
        all_contacts: Dict[str, List[FragmentContactResult]],
        hotspot_contacts: pd.DataFrame,
        output_dir: Path,
):
    """Save contact analysis: per-molecule JSON + summary CSVs."""

    # Per-molecule JSONs
    for name, contact_list in all_contacts.items():
        mol_json = {
            "name": name,
            "fragments": [],
        }
        for fc in contact_list:
            mol_json["fragments"].append({
                "fragment_id": fc.fragment_id,
                "fragment_label": fc.fragment_label,
                "cluster_id": fc.cluster_id,
                "medoid_pose_index": fc.medoid_pose_index,
                "n_contacts": len(fc.contacts),
                "contacts": [
                    {
                        "residue_id": c.residue_id,
                        "residue_name": c.residue_name,
                        "chain": c.chain,
                        "residue_number": c.residue_number,
                        "min_distance": c.min_distance,
                        "n_contacts": c.n_contacts,
                        "contact_type": c.contact_type,
                        "closest_receptor_atom": c.closest_receptor_atom,
                        "closest_ligand_atom": c.closest_ligand_atom,
                    }
                    for c in fc.contacts
                ],
            })
        with open(output_dir / f"{name}_contacts.json", "w") as fh:
            def _json_default(o):
                if isinstance(o, (np.integer,)):
                    return int(o)
                if isinstance(o, (np.floating,)):
                    return float(o)
                if isinstance(o, np.ndarray):
                    return o.tolist()
                return str(o)

            json.dump(mol_json, fh, indent=2, default=_json_default)

    # Flat summary CSV: one row per contact
    rows = []
    for name, contact_list in sorted(all_contacts.items()):
        for fc in contact_list:
            for c in fc.contacts:
                rows.append({
                    "Name": name,
                    "fragment_id": fc.fragment_id,
                    "fragment_label": fc.fragment_label,
                    "cluster_id": fc.cluster_id,
                    "residue_id": c.residue_id,
                    "residue_name": c.residue_name,
                    "chain": c.chain,
                    "residue_number": c.residue_number,
                    "min_distance": c.min_distance,
                    "n_atom_contacts": c.n_contacts,
                    "contact_type": c.contact_type,
                    "closest_receptor_atom": c.closest_receptor_atom,
                    "closest_ligand_atom": c.closest_ligand_atom,
                })

    df = pd.DataFrame(rows)
    if not df.empty:
        df.to_csv(output_dir / "contact_summary.csv", index=False)

    # Hotspot contacts
    if not hotspot_contacts.empty:
        hotspot_contacts.to_csv(output_dir / "hotspot_contacts.csv", index=False)


# =============================================================================
# PLIP INTERACTION ANALYSIS
# =============================================================================

def run_plip_on_medoid(
        receptor_pdb: str,
        ligand_coords: np.ndarray,
        ligand_elements: List[str],
        mol_name: str,
        output_dir: Path,
) -> List[Dict]:
    """
    Run PLIP on receptor + ligand medoid to get geometric interactions.

    Creates temporary PDB: receptor + ligand as HETATM.
    Returns list of interaction dicts.
    """
    if not HAS_PLIP:
        return []

    import tempfile

    # Build ligand HETATM lines
    lig_lines = []
    for i, (coord, elem) in enumerate(zip(ligand_coords, ligand_elements)):
        atom_name = f"{elem}{i+1:>3d}"[:4]
        lig_lines.append(
            f"HETATM{i+1:5d} {atom_name:<4s} LIG A 999    "
            f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
            f"  1.00  0.00          {elem:>2s}"
        )

    # Combine receptor + ligand
    rec_text = Path(receptor_pdb).read_text()
    # Remove END line from receptor
    rec_lines = [l for l in rec_text.split("\n") if not l.startswith("END")]
    combined = "\n".join(rec_lines + lig_lines + ["END"])

    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
            tmp.write(combined)
            tmp_path = tmp.name

        mol = PDBComplex()
        mol.load_pdb(tmp_path)
        mol.analyze()

        interactions = []
        for bsite_id in mol.interaction_sets:
            iset = mol.interaction_sets[bsite_id]

            for hb in iset.hbonds_pdon + iset.hbonds_ldon:
                interactions.append({
                    "type": "hbond",
                    "ligand_atom": hb.atype,
                    "receptor_residue": f"{hb.restype}{hb.resnr}.{hb.reschain}",
                    "receptor_atom": hb.d.type if hasattr(hb, 'd') else "",
                    "distance": round(hb.distance_ad, 2),
                    "angle": round(hb.angle, 1) if hasattr(hb, 'angle') else None,
                })

            for hp in iset.hydrophobic_contacts:
                interactions.append({
                    "type": "hydrophobic",
                    "ligand_atom": hp.ligatom_orig_idx,
                    "receptor_residue": f"{hp.restype}{hp.resnr}.{hp.reschain}",
                    "receptor_atom": hp.bsatom_orig_idx,
                    "distance": round(hp.distance, 2),
                    "angle": None,
                })

            for sb in iset.saltbridge_lneg + iset.saltbridge_pneg:
                interactions.append({
                    "type": "salt_bridge",
                    "ligand_atom": "",
                    "receptor_residue": f"{sb.restype}{sb.resnr}.{sb.reschain}",
                    "receptor_atom": "",
                    "distance": round(sb.distance, 2),
                    "angle": None,
                })

            for ps in iset.pistacking:
                interactions.append({
                    "type": "pi_stacking",
                    "ligand_atom": "",
                    "receptor_residue": f"{ps.restype}{ps.resnr}.{ps.reschain}",
                    "receptor_atom": "",
                    "distance": round(ps.distance, 2),
                    "angle": round(ps.angle, 1) if hasattr(ps, 'angle') else None,
                })

        # Save per-molecule PLIP JSON
        if interactions:
            plip_json = output_dir / f"{mol_name}_plip.json"
            with open(plip_json, "w") as f:
                json.dump({"name": mol_name, "interactions": interactions}, f, indent=2)

        return interactions

    except Exception as e:
        logger.debug(f"  PLIP failed for {mol_name}: {e}")
        return []
    finally:
        try:
            Path(tmp_path).unlink(missing_ok=True)
        except Exception:
            pass


# =============================================================================
# PIPELINE
# =============================================================================

def run_contact_mapping(
        parsed_dir: str,
        cluster_dir: str,
        output_dir: str,
        receptor_path: str,
        hotspot_dir: Optional[str] = None,
        cutoff: float = 4.5,
        max_clusters_per_fragment: int = 3,
        molecule_names: Optional[List[str]] = None,
        plip_enabled: bool = False,
        plip_top_n: int = 20,
        gnina_scores_csv: Optional[str] = None,
) -> Dict[str, any]:
    """
    Contact mapping for all molecules.

    Reads: 05a + 05b + receptor PDB + 05d (optional)
    Saves: per-molecule JSONs + summary CSVs
    """
    logger.info("=" * 60)
    logger.info("05f CONTACT MAPPING v3.0")
    logger.info("=" * 60)

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # Parse receptor
    receptor_atoms = parse_receptor_pdb(receptor_path)
    if not receptor_atoms:
        return {"success": False, "error": "No receptor atoms parsed"}

    receptor_coords = _build_receptor_arrays(receptor_atoms)

    # Load hotspot map (optional)
    hotspot_map = None
    if hotspot_dir:
        hm_path = Path(hotspot_dir) / "hotspot_map.json"
        if hm_path.exists():
            with open(hm_path) as fh:
                hotspot_map = json.load(fh)

    # Load molecules
    all_mols = load_all_molecules(parsed_dir, molecule_names)
    if not all_mols:
        logger.error("No parsed molecules found")
        return {"success": False, "error": "No data"}

    all_contacts = {}
    n_ok = 0

    for name in sorted(all_mols.keys()):
        mol_data = all_mols[name]
        cluster_result = load_cluster_result(name, cluster_dir)
        if cluster_result is None:
            continue

        try:
            contacts = analyze_molecule_contacts(
                mol_data, cluster_result,
                receptor_atoms, receptor_coords,
                cutoff=cutoff,
                max_clusters_per_fragment=max_clusters_per_fragment,
            )
        except Exception as e:
            logger.warning(f"  {name}: contact analysis failed — {e}")
            continue

        all_contacts[name] = contacts
        n_ok += 1

    # Aggregate hotspot contacts
    hotspot_contacts = aggregate_hotspot_contacts(all_contacts, hotspot_map)

    # Save
    save_contact_results(all_contacts, hotspot_contacts, out_path)

    # PLIP analysis for top N molecules
    plip_rows = []
    if plip_enabled and HAS_PLIP:
        # Determine top N by vina_affinity (from gnina_scores.csv) or by contact count
        top_molecules = []
        if gnina_scores_csv and Path(gnina_scores_csv).exists():
            df_gnina = pd.read_csv(gnina_scores_csv)
            if 'name' in df_gnina.columns and 'vina_affinity' in df_gnina.columns:
                df_ok = df_gnina[df_gnina['success'] == True].nsmallest(plip_top_n, 'vina_affinity')
                top_molecules = df_ok['name'].tolist()

        if not top_molecules:
            # Fallback: top N by number of contacts
            if all_contacts:
                mol_contact_counts = {name: sum(len(fc.contacts) for fc in clist)
                                      for name, clist in all_contacts.items()}
                top_molecules = sorted(mol_contact_counts, key=mol_contact_counts.get, reverse=True)[:plip_top_n]

        logger.info(f"  Running PLIP on top {len(top_molecules)} molecules...")
        for mol_name in top_molecules:
            if mol_name not in all_mols:
                continue
            mol_data = all_mols[mol_name]
            cluster_result = load_cluster_result(mol_name, cluster_dir)
            if cluster_result is None or not cluster_result.fragment_results:
                continue

            # Use medoid of dominant cluster of fragment 0 (typically largest ring)
            fr0 = cluster_result.fragment_results[0]
            dom_members = np.where(fr0.labels == fr0.dominant_cluster)[0]
            if len(dom_members) == 0:
                continue
            medoid_idx = _find_medoid(fr0.rmsd_matrix, dom_members)

            # Get ALL atom coords for this pose (not just fragment 0)
            all_coords = mol_data.coords[medoid_idx]  # (n_atoms, 3)
            all_elements = mol_data.heavy_atom_elements

            interactions = run_plip_on_medoid(
                receptor_path, all_coords, all_elements,
                mol_name, out_path,
            )

            for ix in interactions:
                plip_rows.append({"Name": mol_name, **ix})

        if plip_rows:
            df_plip = pd.DataFrame(plip_rows)
            df_plip.to_csv(out_path / "plip_interactions_top20.csv", index=False)
            logger.info(f"  PLIP interactions: {len(plip_rows)} for {df_plip['Name'].nunique()} molecules")
    elif plip_enabled and not HAS_PLIP:
        logger.warning("  PLIP requested but not installed (pip install plip)")

    # Stats
    total_contacts = sum(len(c.contacts) for cl in all_contacts.values() for c in cl)
    logger.info(f"\n{'=' * 60}")
    logger.info("CONTACT MAPPING COMPLETE")
    logger.info(f"{'=' * 60}")
    logger.info(f"  Molecules: {n_ok}")
    logger.info(f"  Total residue contacts: {total_contacts}")
    if not hotspot_contacts.empty:
        logger.info(f"  Hotspot contacts: {len(hotspot_contacts)} rows")
        # Top consensus contacts
        top = hotspot_contacts.nlargest(5, "contact_fraction")
        for _, row in top.iterrows():
            logger.info(f"    H{row['hotspot_id']}: {row['residue_id']} "
                        f"({row['n_molecules_contacting']}/{row['hotspot_n_molecules']} mol, "
                        f"{row['mean_distance']:.1f}A, {row['dominant_contact_type']})")

    return {
        "success": True,
        "n_molecules": n_ok,
        "total_contacts": total_contacts,
        "output_dir": str(out_path),
    }

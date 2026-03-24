"""
Pharmit Pharmacophore Generator (06a) v4.0
============================================
Generate ranked sub-pocket pharmacophore queries from PLIP interactions.

PLIP IS the pharmacophore. Each unique interaction point from PLIP
becomes a pharmacophore feature. No RDKit intermediate layer.

Pipeline:
    1. Load PLIP interactions JSON (from 03a)
    2. Group interactions by ligand coordinates → pharmacophore points
    3. Load DOCK6 footprint → assign energy per point
    4. Filter repulsive (energy > cutoff)
    5. Classify sub-pocket from contacted residues
    6. Generate Pharmit JSONs by strategy + ranking table

Input:
    Required: PLIP interactions JSON (from 03a)
    Optional: DOCK6 residue_consensus.csv (from 04b)
    Optional: receptor PDB (for GNINA rescore)
    Optional: ligand SDF/mol2 (for GNINA rescore)
Output:
    pharmacophore_{strategy}.json  — Pharmit queries (5 strategies)
    pharmacophore_ranking.csv      — Decision table
    pharmacophore_ranking.txt      — Human-readable summary

Location: 01_src/molecular_docking/m06_pharmit/pharmit_pharmacophore.py
Project: molecular_docking
Module: 06a (core)
Version: 4.0
"""

import csv
import json
import logging
import subprocess
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

logger = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

# PLIP interaction type → Pharmit feature type
PLIP_TO_PHARMIT = {
    ("hbond", False):    "HydrogenAcceptor",   # ligand accepts H
    ("hbond", True):     "HydrogenDonor",       # ligand donates H
    ("salt_bridge", "negative"): "NegativeIon",  # ligand is negative
    ("salt_bridge", "positive"): "PositiveIon",  # ligand is positive
    ("pi_stack", None):  "Aromatic",
    ("hydrophobic", None): "Hydrophobic",
    ("pi_cation", None): "Aromatic",
    ("water_bridge", None): "HydrogenAcceptor",
    ("halogen_bond", None): "HydrogenAcceptor",
}

# Sub-pocket classification by receptor residue
XYLOSE_RESIDUES = {"TRP495", "ASP494", "TRP392", "SER575", "TYR565", "CYS574"}
URACIL_RESIDUES = {"ASP361", "THR390", "ARG363"}
PHOSPHATE_RESIDUES = {"ARG598", "LYS599"}
# HIS335 is dual: salt bridge (phosphate) + hbond/pi-stack (ribose/uracil)

# Strategies for Pharmit query generation
STRATEGIES = {
    "xylose":    {"include": ["xylose"], "phosphate_as": None},
    "uracil":    {"include": ["uracil"], "phosphate_as": None},
    "combined":  {"include": ["xylose", "uracil", "ribose"], "phosphate_as": None},
    "analogues": {"include": ["xylose", "uracil", "ribose", "phosphate"], "phosphate_as": "NegativeIon"},
    "druglike":  {"include": ["xylose", "uracil", "ribose", "phosphate"], "phosphate_as": "HydrogenAcceptor"},
}


# ═══════════════════════════════════════════════════════════════════════
# DATA STRUCTURES
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class PharmacophorePoint:
    """One pharmacophore point derived from PLIP interactions."""
    index: int
    pharmit_type: str           # HydrogenAcceptor, HydrogenDonor, etc.
    x: float
    y: float
    z: float

    # PLIP evidence
    residues: List[str] = field(default_factory=list)
    interaction_types: List[str] = field(default_factory=list)
    n_interactions: int = 0
    plip_details: List[Dict] = field(default_factory=list)

    # Direction vector (from ligand → receptor for acceptors, reverse for donors)
    svector: Dict[str, float] = field(default_factory=lambda: {"x": 0.0, "y": 0.0, "z": 1.0})

    # DOCK6 footprint energy
    dock6_energy: float = 0.0
    dock6_vdw: float = 0.0
    dock6_es: float = 0.0
    dock6_freq: float = 0.0

    # GNINA rescore (supplementary)
    gnina_energy: Optional[float] = None

    # Sub-pocket
    sub_pocket: str = "unknown"

    # Ranking
    rank: int = 0
    priority: str = "DISABLED"


# ═══════════════════════════════════════════════════════════════════════
# STEP 1: LOAD PLIP → GROUP → PHARMACOPHORE POINTS
# ═══════════════════════════════════════════════════════════════════════

def _classify_interaction(inter: Dict) -> str:
    """Determine Pharmit feature type from a PLIP interaction."""
    itype = inter.get("interaction_type", "")

    if itype == "hbond":
        is_donor = inter.get("ligand_is_donor", False)
        return PLIP_TO_PHARMIT.get(("hbond", is_donor), "HydrogenAcceptor")

    if itype == "salt_bridge":
        charge = inter.get("ligand_charge", "negative")
        return PLIP_TO_PHARMIT.get(("salt_bridge", charge), "NegativeIon")

    if itype == "pi_stack":
        return "Aromatic"

    if itype == "hydrophobic":
        return "Hydrophobic"

    if itype == "pi_cation":
        return "Aromatic"

    if itype == "water_bridge":
        is_donor = inter.get("ligand_is_donor", False)
        return "HydrogenDonor" if is_donor else "HydrogenAcceptor"

    return "HydrogenAcceptor"


def _compute_direction_vector(lig_coords: List[float],
                              rec_coords: List[float],
                              pharmit_type: str) -> Dict[str, float]:
    """Compute direction vector from PLIP coordinates."""
    if not lig_coords or not rec_coords or len(lig_coords) < 3 or len(rec_coords) < 3:
        return {"x": 0.0, "y": 0.0, "z": 1.0}

    lc = np.array(lig_coords[:3])
    rc = np.array(rec_coords[:3])

    if pharmit_type in ("HydrogenAcceptor", "NegativeIon"):
        # Vector FROM ligand atom TOWARD receptor (lone pair direction)
        vec = rc - lc
    elif pharmit_type == "HydrogenDonor":
        # Vector FROM ligand atom AWAY from receptor (H points toward receptor)
        vec = lc - rc
    elif pharmit_type == "Aromatic":
        # Normal to the plane — use a default for now
        return {"x": 0.0, "y": 0.0, "z": 1.0}
    else:
        vec = rc - lc

    norm = np.linalg.norm(vec)
    if norm < 1e-6:
        return {"x": 0.0, "y": 0.0, "z": 1.0}
    vec /= norm
    return {"x": round(float(vec[0]), 6),
            "y": round(float(vec[1]), 6),
            "z": round(float(vec[2]), 6)}


def load_plip_as_pharmacophore(json_path: str,
                               group_tolerance: float = 1.0,
                               include_hydrophobic: bool = False,
                               ) -> List[PharmacophorePoint]:
    """
    Load PLIP interactions and convert directly to pharmacophore points.

    Groups interactions by ligand coordinates (within tolerance).
    Each group = one pharmacophore point.

    Args:
        json_path: Path to PLIP interactions JSON (from 03a)
        group_tolerance: Distance threshold for grouping coordinates (Angstrom)
        include_hydrophobic: Include hydrophobic contacts as features
    """
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    interactions = data.get("interactions", [])
    logger.info(f"  Loaded {len(interactions)} PLIP interactions")

    # Group by ligand coordinates
    # Each group collects interactions at the same ligand position
    groups = []  # list of (centroid, [interactions])

    for inter in interactions:
        coords = inter.get("ligand_coords", inter.get("coords", []))
        if not coords or len(coords) < 3:
            continue

        itype = inter.get("interaction_type", "")
        if itype == "hydrophobic" and not include_hydrophobic:
            continue

        pos = np.array(coords[:3], dtype=np.float64)

        # Find existing group within tolerance
        merged = False
        for g_centroid, g_interactions in groups:
            if np.linalg.norm(pos - g_centroid) <= group_tolerance:
                g_interactions.append(inter)
                # Update centroid as mean
                all_coords = [np.array(i.get("ligand_coords", i.get("coords", [0,0,0]))[:3])
                              for i in g_interactions]
                g_centroid[:] = np.mean(all_coords, axis=0)
                merged = True
                break

        if not merged:
            groups.append((pos.copy(), [inter]))

    logger.info(f"  Grouped into {len(groups)} pharmacophore points")

    # Convert groups to PharmacophorePoint objects
    points = []
    for idx, (centroid, group_inters) in enumerate(groups):
        # Determine Pharmit type from the strongest interaction in the group
        # Priority: salt_bridge > hbond > pi_stack > hydrophobic
        type_priority = {"salt_bridge": 3, "hbond": 2, "pi_stack": 1,
                         "pi_cation": 1, "hydrophobic": 0, "water_bridge": 1,
                         "halogen_bond": 2}

        best_inter = max(group_inters,
                         key=lambda i: type_priority.get(i.get("interaction_type", ""), 0))
        pharmit_type = _classify_interaction(best_inter)

        # Collect all residues
        residues = []
        interaction_types = []
        for inter in group_inters:
            res = inter.get("residue", "")
            if res and res not in residues:
                residues.append(res)
            it = inter.get("interaction_type", "")
            if it and it not in interaction_types:
                interaction_types.append(it)

        # Direction vector from first interaction with receptor coords
        svector = {"x": 0.0, "y": 0.0, "z": 1.0}
        for inter in group_inters:
            rc = inter.get("receptor_coords", [])
            if rc and len(rc) >= 3:
                svector = _compute_direction_vector(
                    centroid.tolist(), rc, pharmit_type)
                break

        # PLIP details for CSV/TXT
        details = []
        for inter in group_inters:
            details.append({
                "type": inter.get("interaction_type", ""),
                "residue": inter.get("residue", ""),
                "receptor_atom": inter.get("receptor_atom", ""),
                "distance": inter.get("distance", 0.0),
            })

        points.append(PharmacophorePoint(
            index=idx,
            pharmit_type=pharmit_type,
            x=round(float(centroid[0]), 4),
            y=round(float(centroid[1]), 4),
            z=round(float(centroid[2]), 4),
            residues=residues,
            interaction_types=interaction_types,
            n_interactions=len(group_inters),
            plip_details=details,
            svector=svector,
        ))

    # Summary
    type_counts = Counter(p.pharmit_type for p in points)
    for ftype, n in sorted(type_counts.items()):
        logger.info(f"    {ftype}: {n}")

    return points


# ═══════════════════════════════════════════════════════════════════════
# STEP 2: DOCK6 FOOTPRINT ENERGY
# ═══════════════════════════════════════════════════════════════════════

def load_dock6_footprint(csv_path: str) -> Dict[str, Dict[str, float]]:
    """Load DOCK6 residue_consensus.csv → ref_vdw + ref_es per residue."""
    import pandas as pd
    df = pd.read_csv(csv_path)

    energy = {}
    for _, row in df.iterrows():
        rid = row["residue_id"]
        ref_vdw = float(row.get("ref_vdw", 0) or 0)
        ref_es = float(row.get("ref_es", 0) or 0)
        frac = float(row.get("frac_contributing", 0) or 0)

        energy[rid] = {
            "vdw": ref_vdw, "es": ref_es,
            "total": ref_vdw + ref_es, "freq": frac,
        }

    logger.info(f"  Loaded DOCK6 footprint: {len(energy)} residues")
    return energy


def assign_dock6_energy(points: List[PharmacophorePoint],
                        dock6_energy: Dict[str, Dict[str, float]]) -> None:
    """Assign DOCK6 energy to each point from its contacted residues."""
    for pt in points:
        total_vdw, total_es = 0.0, 0.0
        freqs = []

        for res in pt.residues:
            # Try RES###.A format
            for suffix in [".A", ".B", ""]:
                key = f"{res}{suffix}"
                if key in dock6_energy:
                    e = dock6_energy[key]
                    total_vdw += e["vdw"]
                    total_es += e["es"]
                    freqs.append(e["freq"])
                    break

        pt.dock6_vdw = round(total_vdw, 4)
        pt.dock6_es = round(total_es, 4)
        pt.dock6_energy = round(total_vdw + total_es, 4)
        pt.dock6_freq = round(float(np.mean(freqs)), 4) if freqs else 0.0


# ═══════════════════════════════════════════════════════════════════════
# STEP 3: SUB-POCKET CLASSIFICATION
# ═══════════════════════════════════════════════════════════════════════

def classify_sub_pockets(points: List[PharmacophorePoint]) -> None:
    """Classify each point's sub-pocket from its PLIP residues."""
    for pt in points:
        res_set = set(pt.residues)

        # Check overlap with known sub-pocket residues
        n_xylose = len(res_set & XYLOSE_RESIDUES)
        n_uracil = len(res_set & URACIL_RESIDUES)
        n_phosphate = len(res_set & PHOSPHATE_RESIDUES)

        # HIS335 is special: if interaction is salt_bridge → phosphate
        # If hbond or pi_stack → ribose/uracil
        his335_as_phosphate = ("HIS335" in res_set and
                               "salt_bridge" in pt.interaction_types)
        his335_as_base = ("HIS335" in res_set and
                          ("hbond" in pt.interaction_types or
                           "pi_stack" in pt.interaction_types))

        if his335_as_phosphate:
            n_phosphate += 1

        # NegativeIon → always phosphate
        if pt.pharmit_type == "NegativeIon":
            pt.sub_pocket = "phosphate"
        elif n_phosphate > 0 and n_xylose == 0 and n_uracil == 0 and not his335_as_base:
            pt.sub_pocket = "phosphate"
        elif n_xylose > n_uracil:
            pt.sub_pocket = "xylose"
        elif n_uracil > 0:
            pt.sub_pocket = "uracil"
        elif his335_as_base:
            pt.sub_pocket = "ribose"
        else:
            pt.sub_pocket = "ribose"

    counts = Counter(p.sub_pocket for p in points)
    for sp, n in sorted(counts.items()):
        logger.info(f"    {sp}: {n}")


# ═══════════════════════════════════════════════════════════════════════
# STEP 4: RANKING + FILTERING
# ═══════════════════════════════════════════════════════════════════════

def rank_and_filter(points: List[PharmacophorePoint],
                    energy_cutoff: float = 0.0,
                    n_required: int = 3,
                    n_optional: int = 2) -> List[PharmacophorePoint]:
    """Filter repulsive, rank by DOCK6 energy."""
    has_dock6 = any(p.dock6_energy != 0 for p in points)

    # Filter repulsive
    before = len(points)
    if has_dock6:
        filtered = [p for p in points if p.dock6_energy <= energy_cutoff]
    else:
        filtered = list(points)

    n_excluded = before - len(filtered)
    if n_excluded:
        logger.info(f"  Excluded {n_excluded} repulsive features (energy > {energy_cutoff})")

    # Rank
    if has_dock6:
        logger.info("  Ranking by: DOCK6 footprint energy")
        filtered.sort(key=lambda p: p.dock6_energy)
    else:
        logger.info("  Ranking by: number of PLIP interactions")
        filtered.sort(key=lambda p: p.n_interactions, reverse=True)

    for i, pt in enumerate(filtered):
        pt.rank = i + 1
        if i < n_required:
            pt.priority = "REQUIRED"
        elif i < n_required + n_optional:
            pt.priority = "OPTIONAL"
        else:
            pt.priority = "DISABLED"

    return filtered


# ═══════════════════════════════════════════════════════════════════════
# STEP 5: GNINA RESCORE (optional supplementary)
# ═══════════════════════════════════════════════════════════════════════

def _find_gnina() -> Optional[str]:
    for p in ["gnina", "/usr/local/bin/gnina",
              str(Path.home() / "bin" / "gnina"), "/home/bprieto/bin/gnina"]:
        try:
            r = subprocess.run([p, "--version"], capture_output=True,
                               text=True, timeout=10)
            if r.returncode == 0:
                return p
        except (subprocess.SubprocessError, FileNotFoundError):
            continue
    return None


def _ensure_sdf(input_path: str, work_dir: Path) -> Optional[str]:
    path = Path(input_path)
    if path.suffix.lower() in (".sdf", ".mol"):
        return str(path)
    elif path.suffix.lower() == ".mol2":
        sdf_out = str(work_dir / (path.stem + ".sdf"))
        try:
            r = subprocess.run(["obabel", str(path), "-O", sdf_out],
                               capture_output=True, text=True, timeout=30)
            if r.returncode == 0 and Path(sdf_out).exists():
                return sdf_out
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass
        return str(path)  # gnina also accepts mol2
    return str(path)


def run_gnina_rescore(ligand_path: str, receptor_path: str,
                      output_sdf: str, timeout: int = 60) -> bool:
    gnina = _find_gnina()
    if gnina is None:
        logger.warning("  GNINA not found — skipping rescore")
        return False

    cmd = [gnina, "--receptor", receptor_path, "--ligand", ligand_path,
           "--score_only", "--atom_term_data", "--out", output_sdf]
    logger.info(f"  GNINA: scoring {Path(ligand_path).name}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        if result.returncode != 0:
            logger.warning(f"  GNINA failed (rc={result.returncode})")
            return False
        if not Path(output_sdf).exists():
            return False
        for line in result.stdout.split("\n"):
            if "affinity" in line.lower() or "CNNscore" in line:
                logger.info(f"    {line.strip()}")
        return True
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return False


# ═══════════════════════════════════════════════════════════════════════
# STEP 6: PHARMIT JSON GENERATION
# ═══════════════════════════════════════════════════════════════════════

def generate_pharmit_json(points: List[PharmacophorePoint],
                          strategy_name: str,
                          strategy: Dict,
                          radius: float = 1.0,
                          n_required: int = 3,
                          n_optional: int = 2) -> Dict[str, Any]:
    """Generate one Pharmit JSON for a given strategy."""
    include_pockets = strategy["include"]
    phosphate_as = strategy.get("phosphate_as")

    selected = [p for p in points if p.sub_pocket in include_pockets]

    # Re-rank within selection
    has_dock6 = any(p.dock6_energy != 0 for p in selected)
    if has_dock6:
        selected.sort(key=lambda p: p.dock6_energy)
    else:
        selected.sort(key=lambda p: p.n_interactions, reverse=True)

    json_points = []
    n_enabled = 0
    for i, pt in enumerate(selected):
        name = pt.pharmit_type

        # Convert phosphate type if strategy says so
        if pt.sub_pocket == "phosphate" and phosphate_as is not None:
            name = phosphate_as

        enabled = i < (n_required + n_optional)
        if enabled:
            n_enabled += 1

        has_vec = name in ("HydrogenDonor", "HydrogenAcceptor", "Aromatic")

        json_points.append({
            "name": name,
            "hasvec": has_vec,
            "x": pt.x, "y": pt.y, "z": pt.z,
            "radius": radius,
            "enabled": enabled,
            "vector_on": 0,
            "svector": pt.svector,
            "minsize": "", "maxsize": "",
            "selected": False,
        })

    return {
        "points": json_points,
        "_meta": {
            "strategy": strategy_name,
            "n_features": len(json_points),
            "n_enabled": n_enabled,
            "sub_pockets": include_pockets,
            "phosphate_as": phosphate_as,
        },
    }


# ═══════════════════════════════════════════════════════════════════════
# OUTPUT WRITERS
# ═══════════════════════════════════════════════════════════════════════

def write_ranking_csv(points: List[PharmacophorePoint], path: Path):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Rank", "Priority", "Feature", "Sub_Pocket", "x", "y", "z",
            "Residues", "N_Interactions", "Interaction_Types",
            "DOCK6_Energy", "DOCK6_vdW", "DOCK6_ES", "DOCK6_Freq",
            "GNINA_Energy", "Details",
        ])
        for pt in points:
            details = "; ".join(
                f"{d['type']}→{d['residue']}:{d['receptor_atom']}@{d['distance']:.1f}A"
                for d in pt.plip_details
            )
            writer.writerow([
                pt.rank, pt.priority, pt.pharmit_type, pt.sub_pocket,
                pt.x, pt.y, pt.z,
                ";".join(pt.residues),
                pt.n_interactions,
                ";".join(pt.interaction_types),
                pt.dock6_energy, pt.dock6_vdw, pt.dock6_es, pt.dock6_freq,
                pt.gnina_energy if pt.gnina_energy is not None else "",
                details,
            ])
    logger.info(f"  Saved: {path}")


def write_ranking_txt(points: List[PharmacophorePoint],
                      ligand_name: str, strategies: List[str],
                      path: Path):
    w = 90
    has_dock6 = any(p.dock6_energy != 0 for p in points)

    lines = [
        "=" * w,
        "06a PHARMIT PHARMACOPHORE — FEATURE RANKING (v4.0)",
        "=" * w,
        "",
        f"Date:        {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"Ligand:      {ligand_name}",
        f"Source:       PLIP interactions (from 03a)",
        f"Ranking by:  {'DOCK6 footprint energy' if has_dock6 else 'PLIP interaction count'}",
        f"Features:    {len(points)}",
        f"Strategies:  {', '.join(strategies)}",
        "",
    ]

    by_pri = Counter(p.priority for p in points)
    lines.append(f"  REQUIRED: {by_pri.get('REQUIRED', 0)}")
    lines.append(f"  OPTIONAL: {by_pri.get('OPTIONAL', 0)}")
    lines.append(f"  DISABLED: {by_pri.get('DISABLED', 0)}")
    lines.append("")

    lines.append("-" * w)
    lines.append(
        f"{'Rk':>2s}  {'Pri':<8s}  {'Feature':<18s}  {'Pocket':<10s}  "
        f"{'Residues':<25s}  {'N':>2s}  {'DOCK6':>8s}  {'Freq':>5s}"
    )
    lines.append("-" * w)

    for pt in points:
        star = {"REQUIRED": "***", "OPTIONAL": "** ", "DISABLED": "   "}[pt.priority]
        res_str = ", ".join(pt.residues)
        if len(res_str) > 23:
            res_str = res_str[:23] + ".."
        dock6_str = f"{pt.dock6_energy:+8.2f}" if pt.dock6_energy != 0 else "       -"
        freq_str = f"{pt.dock6_freq:5.2f}" if pt.dock6_freq > 0 else "    -"

        lines.append(
            f"{pt.rank:2d}  {star} {pt.priority:<8s}  {pt.pharmit_type:<18s}  "
            f"{pt.sub_pocket:<10s}  {res_str:<25s}  {pt.n_interactions:2d}  "
            f"{dock6_str}  {freq_str}"
        )

        for d in pt.plip_details:
            lines.append(
                f"{'':>38s}{d['type']} → {d['residue']}:{d['receptor_atom']} "
                f"@ {d['distance']:.1f}Å"
            )

    lines.extend(["", "=" * w])
    path.write_text("\n".join(lines), encoding="utf-8")
    logger.info(f"  Saved: {path}")


# ═══════════════════════════════════════════════════════════════════════
# MAIN PIPELINE
# ═══════════════════════════════════════════════════════════════════════

def generate_pharmit_pharmacophore(
        plip_json_path: str,
        output_dir: str,
        output_name: str = "pharmacophore",
        dock6_csv_path: Optional[str] = None,
        ligand_path: Optional[str] = None,
        receptor_path: Optional[str] = None,
        include_hydrophobic: bool = False,
        energy_cutoff: float = 0.0,
        n_required: int = 3,
        n_optional: int = 2,
        radius: float = 1.0,
        group_tolerance: float = 1.0,
        strategies: Optional[List[str]] = None,
        gnina_timeout: int = 60,
) -> Dict[str, Any]:
    """
    Full pipeline: PLIP JSON → ranked sub-pocket Pharmit JSONs.

    PLIP IS the pharmacophore. Each interaction point becomes a feature.
    DOCK6 footprint provides energy ranking. GNINA optional.
    """
    logger.info("=" * 60)
    logger.info("PHARMIT PHARMACOPHORE GENERATOR (06a) v4.0")
    logger.info("=" * 60)

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"  PLIP:     {Path(plip_json_path).name}")
    if dock6_csv_path:
        logger.info(f"  DOCK6:    {Path(dock6_csv_path).name}")

    if strategies is None:
        strategies = list(STRATEGIES.keys())

    # Get ligand name from PLIP JSON
    with open(plip_json_path) as f:
        plip_meta = json.load(f)
    ligand_name = plip_meta.get("ligand_name", "unknown")

    # ── Step 1: PLIP → pharmacophore points ──
    logger.info("")
    logger.info("Step 1: PLIP → pharmacophore points")
    points = load_plip_as_pharmacophore(
        plip_json_path,
        group_tolerance=group_tolerance,
        include_hydrophobic=include_hydrophobic,
    )

    if not points:
        return {"success": False, "error": "No pharmacophore points from PLIP"}

    # ── Step 2: DOCK6 footprint energy ──
    logger.info("")
    if dock6_csv_path and Path(dock6_csv_path).exists():
        logger.info("Step 2: DOCK6 footprint energy")
        dock6_energy = load_dock6_footprint(dock6_csv_path)
        assign_dock6_energy(points, dock6_energy)

        n_with_energy = sum(1 for p in points if p.dock6_energy != 0)
        logger.info(f"  Energy assigned to {n_with_energy}/{len(points)} points")
    else:
        logger.info("Step 2: DOCK6 — not provided")

    # ── Step 3: Sub-pocket classification ──
    logger.info("")
    logger.info("Step 3: Sub-pocket classification")
    classify_sub_pockets(points)

    # ── Step 4: GNINA rescore (optional) ──
    if ligand_path and receptor_path:
        lig_file = Path(ligand_path)
        rec_file = Path(receptor_path)
        if lig_file.exists() and rec_file.exists():
            logger.info("")
            logger.info("Step 4: GNINA rescore (supplementary)")
            ligand_sdf = _ensure_sdf(str(lig_file), out_dir)
            rescore_sdf = str(out_dir / f"{output_name}_rescore.sdf")
            run_gnina_rescore(ligand_sdf, str(rec_file), rescore_sdf, gnina_timeout)

    # ── Step 5: Rank and filter ──
    logger.info("")
    logger.info("Step 5: Ranking + filtering")
    points = rank_and_filter(points, energy_cutoff, n_required, n_optional)

    # ── Step 6: Generate outputs ──
    logger.info("")
    logger.info("Step 6: Generating outputs")

    csv_path = out_dir / f"{output_name}_ranking.csv"
    write_ranking_csv(points, csv_path)

    txt_path = out_dir / f"{output_name}_ranking.txt"
    write_ranking_txt(points, ligand_name, strategies, txt_path)

    generated = {}
    for strat_name in strategies:
        if strat_name not in STRATEGIES:
            continue
        strat = STRATEGIES[strat_name]
        pharmit_data = generate_pharmit_json(
            points, strat_name, strat,
            radius=radius, n_required=n_required, n_optional=n_optional,
        )
        json_path = out_dir / f"{output_name}_{strat_name}.json"
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(pharmit_data, f, indent=4)
        meta = pharmit_data["_meta"]
        logger.info(f"  {strat_name}: {meta['n_features']} features "
                     f"({meta['n_enabled']} enabled) → {json_path.name}")
        generated[strat_name] = str(json_path)

    # Summary
    logger.info("")
    logger.info("=" * 60)
    logger.info(f"  DONE — {len(points)} points, {len(generated)} queries")
    logger.info("=" * 60)

    for pt in points[:n_required + n_optional]:
        res = ", ".join(pt.residues[:3])
        dock6 = f"DOCK6={pt.dock6_energy:+.1f}" if pt.dock6_energy != 0 else ""
        logger.info(f"  #{pt.rank} [{pt.priority}] {pt.pharmit_type} "
                     f"({pt.sub_pocket}) → {res} {dock6}")

    return {
        "success": True,
        "ranking_csv": str(csv_path),
        "ranking_txt": str(txt_path),
        "queries": generated,
        "n_features": len(points),
        "strategies": strategies,
        "features": [
            {
                "rank": p.rank, "priority": p.priority,
                "name": p.pharmit_type, "sub_pocket": p.sub_pocket,
                "x": p.x, "y": p.y, "z": p.z,
                "residues": p.residues,
                "n_interactions": p.n_interactions,
                "dock6_energy": p.dock6_energy,
            }
            for p in points
        ],
    }
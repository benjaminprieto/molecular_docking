"""
07b Hit Export - Core
=======================
Read a selection Excel, match molecules against gnina_scores.csv,
and export protonated SDF files for hit_validation.

Project: molecular_docking
Module: 07b
Version: 1.0
"""

import logging
import shutil
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd

logger = logging.getLogger(__name__)


def _find_header_row(xlsx_path: str, max_rows: int = 10) -> int:
    """Scan first rows of Excel to find the one containing 'Name' as a value."""
    df_raw = pd.read_excel(xlsx_path, header=None, nrows=max_rows)
    for idx, row in df_raw.iterrows():
        if any(str(v).strip() == "Name" for v in row.values):
            return int(idx)
    raise ValueError(
        f"Could not find a header row containing 'Name' in the first "
        f"{max_rows} rows of {xlsx_path}"
    )


def run_hit_export(
    selection_xlsx: str,
    gnina_scores_csv: str,
    output_dir: str,
) -> Dict[str, Any]:
    """Export selected molecules as individual SDF files.

    Parameters
    ----------
    selection_xlsx : str
        Absolute path to the selection Excel file.
    gnina_scores_csv : str
        Path to gnina_scores.csv (output of module 02c).
    output_dir : str
        Directory for outputs (sdf/ subfolder + summary CSV).

    Returns
    -------
    dict with keys: success, n_selected, n_found, n_exported, n_missing,
                    missing, summary_csv, sdf_dir
    """
    selection_path = Path(selection_xlsx)
    scores_path = Path(gnina_scores_csv)
    out_path = Path(output_dir)

    # --- Validate inputs ---
    if not scores_path.exists():
        raise FileNotFoundError(
            f"gnina_scores.csv not found: {scores_path}\n"
            "Run module 02c (gnina_scores) first."
        )
    if not selection_path.exists():
        raise FileNotFoundError(f"Selection Excel not found: {selection_path}")

    # --- Read selection Excel ---
    header_row = _find_header_row(str(selection_path))
    df_sel = pd.read_excel(str(selection_path), header=header_row)
    # Normalize column names
    df_sel.columns = [str(c).strip() for c in df_sel.columns]
    if "Name" not in df_sel.columns:
        raise ValueError(
            f"Selection Excel has no 'Name' column. "
            f"Found columns: {list(df_sel.columns)}"
        )
    selected_names: List[str] = (
        df_sel["Name"].dropna().astype(str).str.strip().tolist()
    )
    selected_names = [n for n in selected_names if n not in ("UDX", "UDX_control", "UDX_reference")]
    n_selected = len(selected_names)
    logger.info(f"Selection Excel: {n_selected} molecules from {selection_path.name}")

    # --- Read gnina_scores.csv ---
    df_scores = pd.read_csv(scores_path)
    logger.info(f"gnina_scores.csv: {len(df_scores)} rows")

    # Build lookup: name -> row (prefer non-duplicate entries)
    # Separate originals from duplicates
    is_dup = df_scores["SMILES_Duplicate_Of"].notna() if "SMILES_Duplicate_Of" in df_scores.columns else pd.Series(False, index=df_scores.index)
    df_originals = df_scores[~is_dup].set_index("name")
    df_duplicates = df_scores[is_dup]

    # Build duplicate -> original mapping
    dup_to_original = {}
    if "SMILES_Duplicate_Of" in df_scores.columns:
        for _, row in df_duplicates.iterrows():
            dup_to_original[str(row["name"])] = str(row["SMILES_Duplicate_Of"])

    # --- Prepare output dirs ---
    sdf_dir = out_path / "sdf"
    sdf_dir.mkdir(parents=True, exist_ok=True)

    # --- Process each selected molecule ---
    summary_rows = []
    missing = []
    n_found = 0
    n_exported = 0

    for mol_name in selected_names:
        # Look up in gnina_scores
        score_row = None

        if mol_name in df_originals.index:
            score_row = df_originals.loc[mol_name]
        elif mol_name in dup_to_original:
            # This molecule is a SMILES duplicate — use the original's data
            original_name = dup_to_original[mol_name]
            if original_name in df_originals.index:
                score_row = df_originals.loc[original_name]
                logger.info(f"  {mol_name} is SMILES duplicate of {original_name}, using original data")
        else:
            # Try reverse: maybe mol_name is the original of a duplicate
            # Search duplicates whose original is mol_name
            pass  # Already covered by df_originals lookup

        if score_row is None:
            logger.warning(f"  {mol_name}: NOT FOUND in gnina_scores.csv")
            # Try fallback: mol2 from 01a_antechamber
            campaign_results = Path(gnina_scores_csv).parent.parent
            mol2_fallback = campaign_results / "01a_antechamber" / "mol2" / f"{mol_name}.mol2"
            if mol2_fallback.exists():
                mol2_dest = sdf_dir / f"{mol_name}.mol2"
                shutil.copy2(mol2_fallback, mol2_dest)
                n_exported += 1
                logger.info(f"  {mol_name}: exported mol2 (fallback from 01a, no gnina scores)")
                summary_rows.append({
                    "name": mol_name,
                    "SMILES": "",
                    "MW": "",
                    "TPSA": "",
                    "QED": "",
                    "vina_affinity": "",
                    "cnn_score": "",
                    "source_file": "",
                    "sdf_exported": True,
                    "sdf_path": str(mol2_dest),
                })
            else:
                missing.append(mol_name)
                summary_rows.append({
                    "name": mol_name,
                    "SMILES": "",
                    "MW": "",
                    "TPSA": "",
                    "QED": "",
                    "vina_affinity": "",
                    "cnn_score": "",
                    "source_file": "",
                    "sdf_exported": False,
                    "sdf_path": "",
                })
            continue

        n_found += 1

        # Handle potential Series vs scalar for duplicate index
        if isinstance(score_row, pd.DataFrame):
            score_row = score_row.iloc[0]

        ligand_file = str(score_row.get("ligand_file", ""))
        smiles = str(score_row.get("SMILES", ""))
        mw = score_row.get("MW", "")
        tpsa = score_row.get("TPSA", "")
        qed = score_row.get("QED", "")
        vina = score_row.get("vina_affinity", "")
        cnn = score_row.get("cnn_score", "")
        source = str(score_row.get("Source_File", ""))

        # Copy SDF
        sdf_exported = False
        sdf_dest = sdf_dir / f"{mol_name}.sdf"

        if ligand_file and Path(ligand_file).exists():
            shutil.copy2(ligand_file, sdf_dest)
            sdf_exported = True
            n_exported += 1
            logger.info(f"  {mol_name}: exported SDF")
        else:
            # Fallback: mol2 from 01a_antechamber
            campaign_results = Path(gnina_scores_csv).parent.parent  # 05_results/{campaign_id}/
            mol2_fallback = campaign_results / "01a_antechamber" / "mol2" / f"{mol_name}.mol2"
            if mol2_fallback.exists():
                mol2_dest = sdf_dir / f"{mol_name}.mol2"
                shutil.copy2(mol2_fallback, mol2_dest)
                sdf_exported = True
                n_exported += 1
                logger.info(f"  {mol_name}: exported mol2 (fallback from 01a)")
            else:
                logger.warning(f"  {mol_name}: no SDF (00c) nor mol2 (01a) found")

        summary_rows.append({
            "name": mol_name,
            "SMILES": smiles,
            "MW": mw,
            "TPSA": tpsa,
            "QED": qed,
            "vina_affinity": vina,
            "cnn_score": cnn,
            "source_file": source,
            "sdf_exported": sdf_exported,
            "sdf_path": str(sdf_dest) if sdf_exported else "",
        })

    # --- Write summary CSV ---
    summary_csv_path = out_path / "hit_export_summary.csv"
    df_summary = pd.DataFrame(summary_rows)
    df_summary.to_csv(summary_csv_path, index=False)
    logger.info(f"Summary written: {summary_csv_path}")

    n_missing = len(missing)
    success = n_exported > 0

    logger.info(f"Results: {n_selected} selected, {n_found} found, "
                f"{n_exported} exported, {n_missing} missing")
    if missing:
        logger.warning(f"Missing molecules: {missing}")

    return {
        "success": success,
        "n_selected": n_selected,
        "n_found": n_found,
        "n_exported": n_exported,
        "n_missing": n_missing,
        "missing": missing,
        "summary_csv": str(summary_csv_path),
        "sdf_dir": str(sdf_dir),
    }

"""
Score Collector - Core Module (02a)
=====================================
Parsea scored mol2 de DOCK6, construye Excel dock2profile-compatible,
y separa mol2 de best pose por molecula.

Location: 01_src/molecular_docking/m02_collection/score_collector.py
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Union

logger = logging.getLogger(__name__)


def run_score_collection(
        docking_dir: Union[str, Path],
        molecules_csv: Union[str, Path],
        output_dir: Union[str, Path],
        score_key: str = "Grid_Score",
        max_molecules: int = 500,
        extract_best_pose_mol2: bool = True,
        keep_all_poses: bool = False,
        compute_properties: bool = True,
        scores_filename: str = "01_top_500_molecules.xlsx",
        mol2_dirname: str = "docked_molecules",
        source_label: Optional[str] = None,
) -> Dict[str, Any]:
    """Run the complete score collection pipeline."""
    raise NotImplementedError("Module 02a core — to be implemented")


def parse_scored_mol2(mol2_path: str) -> List[Dict[str, Any]]:
    """Parse DOCK6 scored mol2, extract scores for each pose."""
    raise NotImplementedError


def get_best_pose(poses: List[Dict], score_key: str) -> Optional[Dict]:
    """Get pose with best (lowest) score."""
    raise NotImplementedError


def extract_single_pose_mol2(scored_mol2: str, pose_index: int,
                              output_mol2: str) -> bool:
    """Extract single pose from multi-pose scored mol2."""
    raise NotImplementedError

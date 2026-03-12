#!/usr/bin/env python3
"""
setup_project.py - Bootstrap molecular_docking project
========================================================
Creates the complete project structure: directories, __init__.py files,
config YAMLs, CLI scripts, core module stubs, tests, and packaging files.

Run ONCE from the project root directory:

    cd ~/projects/molecular_docking
    python setup_project.py

After running:
    conda activate molecular_docking_env
    conda env update -f environment.yaml --prune
    pip install -e ".[dev]"
    bash check_dependencies.sh

Project: molecular_docking
Version: 1.0.0
"""

import os
from pathlib import Path

# =============================================================================
# PROJECT ROOT
# =============================================================================

ROOT = Path(__file__).parent
VERSION = "1.0.0"


def write_file(relative_path: str, content: str):
    """Create a file with content, creating parent dirs as needed."""
    path = ROOT / relative_path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    print(f"  Created: {relative_path}")


# =============================================================================
# 1. DIRECTORIES (with .gitkeep for empty ones)
# =============================================================================

def create_directories():
    """Create all project directories."""
    print("\n=== Creating directories ===")
    dirs = [
        "01_src/molecular_docking/m00_preparation",
        "01_src/molecular_docking/m01_docking",
        "01_src/molecular_docking/m02_collection",
        "02_scripts",
        "03_configs",
        "04_data/campaigns/example_campaign/receptor",
        "04_data/campaigns/example_campaign/molecules",
        "04_data/campaigns/example_campaign/grids",
        "05_results",
        "docs",
        "tests",
    ]
    for d in dirs:
        (ROOT / d).mkdir(parents=True, exist_ok=True)
        print(f"  Dir: {d}/")

    # .gitkeep for empty dirs
    for gk in [
        "04_data/campaigns/example_campaign/receptor/.gitkeep",
        "04_data/campaigns/example_campaign/molecules/.gitkeep",
        "04_data/campaigns/example_campaign/grids/.gitkeep",
        "05_results/.gitkeep",
    ]:
        (ROOT / gk).touch()


# =============================================================================
# 2. PACKAGING FILES
# =============================================================================

def create_packaging_files():
    """Create pyproject.toml, setup.py, environment.yaml, etc."""
    print("\n=== Creating packaging files ===")

    # --- pyproject.toml ---
    write_file("pyproject.toml", f'''\
[build-system]
requires = ["setuptools>=68.0", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "molecular_docking"
version = "{VERSION}"
description = "Generic DOCK6 molecular docking pipeline. Produces outputs compatible with dock2profile."
readme = "README.md"
license = "MIT"
requires-python = ">=3.9"
authors = [
    {{name = "Environ Bio"}},
]
keywords = ["molecular-docking", "dock6", "drug-discovery", "virtual-screening"]
classifiers = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: Chemistry",
    "Programming Language :: Python :: 3",
    "Operating System :: POSIX :: Linux",
]

dependencies = [
    "pandas>=1.5.0",
    "numpy>=1.23.0",
    "openpyxl>=3.0.0",
    "pyyaml>=6.0",
]

[project.optional-dependencies]
protonation = ["dimorphite_dl>=1.3.2", "pdb2pqr>=3.6.0"]
dev = ["pytest>=7.0", "pytest-cov>=4.0"]

[tool.setuptools.packages.find]
where = ["01_src"]

[tool.setuptools.package-dir]
"" = "01_src"

[tool.pytest.ini_options]
testpaths = ["tests"]
addopts = "-v --tb=short"
''')

    # --- setup.py ---
    write_file("setup.py", '''\
#!/usr/bin/env python3
"""
Backwards-compatible setup.py. Metadata lives in pyproject.toml.

    pip install -e ".[dev]"
"""
from setuptools import setup, find_packages
from pathlib import Path

version = "''' + VERSION + '''"
init_file = Path(__file__).parent / "01_src" / "molecular_docking" / "__init__.py"
if init_file.exists():
    for line in init_file.read_text().splitlines():
        if line.startswith("__version__"):
            version = line.split('"')[1]
            break

readme = Path(__file__).parent / "README.md"
long_description = readme.read_text(encoding="utf-8") if readme.exists() else ""

setup(
    name="molecular_docking",
    version=version,
    description="Generic DOCK6 molecular docking pipeline compatible with dock2profile",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Environ Bio",
    license="MIT",
    python_requires=">=3.9",
    package_dir={"": "01_src"},
    packages=find_packages(where="01_src"),
)
''')

    # --- environment.yaml ---
    write_file("environment.yaml", '''\
# =============================================================================
# molecular_docking - Conda Environment
# =============================================================================
# conda env create -f environment.yaml
# conda activate molecular_docking_env
# pip install -e ".[dev]"
# =============================================================================

name: molecular_docking_env

channels:
  - conda-forge
  - defaults

dependencies:
  - python>=3.9,<3.13
  - pandas>=1.5
  - numpy>=1.23
  - scipy>=1.9
  - openpyxl>=3.0
  - pyyaml>=6.0
  - rdkit>=2023.03.1
  - openbabel>=3.1.1
  - ambertools>=22.0
  - matplotlib>=3.6
  - seaborn>=0.12
  - pytest>=7.0
  - pytest-cov>=4.0
  - tqdm>=4.64
  - lxml>=4.9
  - pip
  - pip:
      - dimorphite_dl>=1.3.2
      - pdb2pqr>=3.6.0
''')

    # --- requirements.txt ---
    write_file("requirements.txt", '''\
# pip-only dependencies. For full install use: conda env create -f environment.yaml
pandas>=1.5.0
numpy>=1.23.0
scipy>=1.9.0
openpyxl>=3.0.0
pyyaml>=6.0
tqdm>=4.64.0
lxml>=4.9.0
dimorphite_dl>=1.3.2
pdb2pqr>=3.6.0
matplotlib>=3.6.0
seaborn>=0.12.0
pytest>=7.0
pytest-cov>=4.0
''')

    # --- MANIFEST.in ---
    write_file("MANIFEST.in", '''\
include README.md
include pyproject.toml
include environment.yaml
include requirements.txt
include check_dependencies.sh
recursive-include 02_scripts *.py
recursive-include 03_configs *.yaml
recursive-include 04_data *.yaml .gitkeep
recursive-include docs *.md
recursive-include tests *.py
''')

    # --- .gitignore ---
    write_file(".gitignore", '''\
05_results/*/
__pycache__/
*.pyc
*.egg-info/
dist/
build/
.idea/
.vscode/
.DS_Store
''')


# =============================================================================
# 3. __init__.py FILES
# =============================================================================

def create_init_files():
    """Create package __init__.py files."""
    print("\n=== Creating __init__.py files ===")

    write_file("01_src/molecular_docking/__init__.py", f'''\
"""
molecular_docking - DOCK6 Docking Pipeline
============================================
Generic pipeline for molecular docking with DOCK6.
Produces outputs compatible with dock2profile.

Modules:
    m00_preparation  - Parse molecules, prepare receptor & ligands
    m01_docking      - Grid generation & DOCK6 execution
    m02_collection   - Score collection & Excel generation
"""
__version__ = "{VERSION}"
''')

    write_file("01_src/molecular_docking/m00_preparation/__init__.py",
               '"""Phase 0: Preparation (molecule parsing, receptor/ligand prep)."""\n')

    write_file("01_src/molecular_docking/m01_docking/__init__.py",
               '"""Phase 1: Docking (grid generation, DOCK6 execution)."""\n')

    write_file("01_src/molecular_docking/m02_collection/__init__.py",
               '"""Phase 2: Collection (score parsing, Excel generation)."""\n')

    write_file("tests/__init__.py", "")


# =============================================================================
# 4. CORE MODULE STUBS
# =============================================================================

def create_core_modules():
    """Create core module stubs with function signatures."""
    print("\n=== Creating core modules ===")

    # --- 00a molecule_parser ---
    write_file("01_src/molecular_docking/m00_preparation/molecule_parser.py", '''\
"""
Molecule Parser - Core Module (00a)
=====================================
Parsea moleculas de cualquier formato → tabla normalizada + SDF limpio.
Formatos: SDF, CSV/XLSX, mol2 directory, SMILES txt.

Location: 01_src/molecular_docking/m00_preparation/molecule_parser.py
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Union

logger = logging.getLogger(__name__)


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
    """Run the complete molecule parsing pipeline."""
    raise NotImplementedError("Module 00a core — to be implemented")


def detect_input_format(input_file: Union[str, Path]) -> str:
    """Detect format: 'sdf' | 'csv' | 'xlsx' | 'smiles' | 'mol2_dir'."""
    raise NotImplementedError


def parse_sdf(sdf_path: str, name_source: str, name_property: str) -> List[Dict]:
    """Parse SDF file into list of molecule dicts."""
    raise NotImplementedError


def resolve_conformers(molecules: List[Dict], strategy: str) -> List[Dict]:
    """Resolve multiple conformers per molecule to single representative."""
    raise NotImplementedError


def compute_molecular_properties(mol) -> Dict[str, float]:
    """Compute MW, LogP, HBD, HBA, TPSA, QED, etc."""
    raise NotImplementedError
''')

    # --- 00b receptor_preparation ---
    write_file("01_src/molecular_docking/m00_preparation/receptor_preparation.py", '''\
"""
Receptor Preparation - Core Module (00b)
==========================================
Prepara receptor para DOCK6: protonacion + cargas + mol2.
OPCIONAL: skip si receptor.protonation.enabled=false.

Location: 01_src/molecular_docking/m00_preparation/receptor_preparation.py
"""

import logging
from pathlib import Path
from typing import Dict, Optional, Any, Union

logger = logging.getLogger(__name__)


def run_receptor_preparation(
        receptor_pdb: Union[str, Path],
        output_dir: Union[str, Path],
        docking_ph: float = 7.2,
        protonation_tool: str = "chimera",
        force_field: str = "AMBER",
        chain: Optional[str] = None,
        remove_water: bool = True,
        remove_hetatm: bool = True,
        remove_alt_conformations: bool = True,
) -> Dict[str, Any]:
    """Run receptor preparation pipeline."""
    raise NotImplementedError("Module 00b core — to be implemented")


def validate_prepared_mol2(mol2_path: Union[str, Path]) -> Dict[str, Any]:
    """Validate existing receptor mol2 for DOCK6 compatibility."""
    raise NotImplementedError
''')

    # --- 00c ligand_preparation ---
    write_file("01_src/molecular_docking/m00_preparation/ligand_preparation.py", '''\
"""
Ligand Preparation - Core Module (00c)
========================================
Protona ligandos + genera mol2 con cargas AM1-BCC via antechamber.
Protonacion OPCIONAL segun campaign_config.

Location: 01_src/molecular_docking/m00_preparation/ligand_preparation.py
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Union

logger = logging.getLogger(__name__)


def run_ligand_preparation(
        molecules_csv: Union[str, Path],
        output_dir: Union[str, Path],
        docking_ph: float = 7.2,
        protonate_enabled: bool = True,
        protonate_tool: str = "dimorphite_dl",
        obabel_fallback: bool = True,
        charge_method: str = "bcc",
        atom_type: str = "gaff2",
        antechamber_timeout: int = 300,
        n_3d_attempts: int = 5,
) -> Dict[str, Any]:
    """Run ligand preparation for all molecules."""
    raise NotImplementedError("Module 00c core — to be implemented")


def protonate_molecule(smiles: str, ph: float, tool: str, fallback: bool) -> str:
    """Protonate SMILES at given pH."""
    raise NotImplementedError


def generate_3d_conformer(smiles: str, n_attempts: int = 5):
    """Generate 3D conformer from SMILES using RDKit ETKDG."""
    raise NotImplementedError


def run_antechamber(input_sdf: Path, output_mol2: Path, charge_method: str,
                    atom_type: str, net_charge: int, timeout: int) -> bool:
    """Run antechamber to generate mol2 with AM1-BCC charges."""
    raise NotImplementedError


def validate_mol2(mol2_path: Path) -> Dict[str, Any]:
    """Validate mol2 file: atoms, bonds, charges."""
    raise NotImplementedError
''')

    # --- 01a grid_generation ---
    write_file("01_src/molecular_docking/m01_docking/grid_generation.py", '''\
"""
Grid Generation - Core Module (01a)
=====================================
Genera grids DOCK6. OPCIONAL: skip si grids pre-existentes son validos.
Pipeline: DMS → sphgen → sphere_selector → showbox → grid

Location: 01_src/molecular_docking/m01_docking/grid_generation.py
"""

import logging
from pathlib import Path
from typing import Dict, Optional, Any, Union, List

logger = logging.getLogger(__name__)


def run_grid_generation(
        receptor_noH_pdb: Union[str, Path],
        receptor_charged_mol2: Union[str, Path],
        output_dir: Union[str, Path],
        box_center: Optional[List[float]] = None,
        box_radius: float = 15.0,
        grid_spacing: float = 0.3,
        probe_radius: float = 1.4,
        max_spheres: int = 50,
) -> Dict[str, Any]:
    """Run complete DOCK6 grid generation pipeline."""
    raise NotImplementedError("Module 01a core — to be implemented")


def validate_existing_grids(grid_dir: Union[str, Path],
                            spheres_file: str, energy_grid: str,
                            bump_grid: str) -> bool:
    """Check if pre-existing grids are valid and complete."""
    grid_dir = Path(grid_dir)
    required = [
        grid_dir / spheres_file,
        grid_dir / energy_grid,
        grid_dir / bump_grid,
    ]
    return all(f.exists() and f.stat().st_size > 0 for f in required)
''')

    # --- 01b dock6_runner ---
    write_file("01_src/molecular_docking/m01_docking/dock6_runner.py", '''\
"""
DOCK6 Runner - Core Module (01b)
==================================
Genera input files de DOCK6 y ejecuta docking para cada ligando.
Soporta flex (anchor-and-grow) y rigid docking.

Location: 01_src/molecular_docking/m01_docking/dock6_runner.py
"""

import logging
import os
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Any, Union

logger = logging.getLogger(__name__)


# =============================================================================
# DOCK6 INPUT TEMPLATES
# =============================================================================

DOCK6_FLEX_TEMPLATE = """\\
conformer_search_type                                        flex
user_specified_anchor                                        no
limit_max_anchors                                            no
min_anchor_size                                              {min_anchor_size}
pruning_use_clustering                                       yes
pruning_max_orients                                          {pruning_max_orients}
pruning_clustering_cutoff                                    {pruning_clustering_cutoff}
pruning_conformer_score_cutoff                               {pruning_conformer_score_cutoff}
pruning_conformer_score_scaling_factor                       1.0
use_clash_overlap                                            no
write_growth_tree                                            no
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             {ligand_mol2}
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                yes
automated_matching                                           yes
receptor_site_file                                           {spheres_file}
max_orientations                                             {max_orientations}
critical_points                                              no
chemical_matching                                            no
use_ligand_spheres                                           no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           yes
grid_score_secondary                                         no
grid_score_rep_rad_scale                                     1
grid_score_vdw_scale                                         1
grid_score_es_scale                                          1
grid_score_grid_prefix                                       {grid_prefix}
multigrid_score_secondary                                    no
dock3.5_score_secondary                                      no
continuous_score_secondary                                   no
footprint_similarity_score_secondary                         no
pharmacophore_score_secondary                                no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              {minimize}
simplex_max_iterations                                       {simplex_max_iterations}
simplex_tors_premin_iterations                               0
simplex_max_cycles                                           {simplex_max_cycles}
simplex_score_converge                                       {simplex_score_converge}
simplex_cycle_converge                                       {simplex_cycle_converge}
simplex_trans_step                                           {simplex_trans_step}
simplex_rot_step                                             {simplex_rot_step}
simplex_tors_step                                            {simplex_tors_step}
simplex_random_seed                                          0
simplex_restraint_min                                        no
atom_model                                                   all
vdw_defn_file                                                {vdw_defn_file}
flex_defn_file                                               {flex_defn_file}
flex_drive_file                                              {flex_drive_file}
ligand_outfile_prefix                                        {output_prefix}
write_orientations                                           {write_orientations}
num_scored_conformers                                        {num_scored_conformers}
rank_ligands                                                 no
"""

DOCK6_RIGID_TEMPLATE = """\\
conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             {ligand_mol2}
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                yes
automated_matching                                           yes
receptor_site_file                                           {spheres_file}
max_orientations                                             {max_orientations}
critical_points                                              no
chemical_matching                                            no
use_ligand_spheres                                           no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           yes
grid_score_secondary                                         no
grid_score_rep_rad_scale                                     1
grid_score_vdw_scale                                         1
grid_score_es_scale                                          1
grid_score_grid_prefix                                       {grid_prefix}
multigrid_score_secondary                                    no
dock3.5_score_secondary                                      no
continuous_score_secondary                                   no
footprint_similarity_score_secondary                         no
pharmacophore_score_secondary                                no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              {minimize}
simplex_max_iterations                                       {simplex_max_iterations}
simplex_max_cycles                                           {simplex_max_cycles}
simplex_score_converge                                       {simplex_score_converge}
simplex_cycle_converge                                       {simplex_cycle_converge}
simplex_trans_step                                           {simplex_trans_step}
simplex_rot_step                                             {simplex_rot_step}
simplex_tors_step                                            {simplex_tors_step}
simplex_random_seed                                          0
simplex_restraint_min                                        no
atom_model                                                   all
vdw_defn_file                                                {vdw_defn_file}
flex_defn_file                                               {flex_defn_file}
flex_drive_file                                              {flex_drive_file}
ligand_outfile_prefix                                        {output_prefix}
write_orientations                                           {write_orientations}
num_scored_conformers                                        {num_scored_conformers}
rank_ligands                                                 no
"""


# =============================================================================
# PUBLIC API
# =============================================================================

def run_dock6_batch(
        ligand_mol2_dir: Union[str, Path],
        spheres_file: Union[str, Path],
        grid_prefix: str,
        output_dir: Union[str, Path],
        search_method: str = "flex",
        max_orientations: int = 1000,
        num_scored_conformers: int = 20,
        minimize: bool = True,
        simplex_max_iterations: int = 500,
        timeout_per_molecule: int = 600,
        molecule_filter: Optional[List[str]] = None,
        dry_run: bool = False,
        **kwargs,
) -> Dict[str, Any]:
    """Run DOCK6 for all prepared ligands."""
    raise NotImplementedError("Module 01b core — to be implemented")


def find_dock6_params() -> Dict[str, str]:
    """Find DOCK6 parameter files (vdw_AMBER_parm99.defn, etc.)."""
    param_files = {
        "vdw_defn_file": "vdw_AMBER_parm99.defn",
        "flex_defn_file": "flex.defn",
        "flex_drive_file": "flex_drive.tbl",
    }
    search_paths = []
    for var in ["DOCK_HOME", "DOCK6_HOME", "DOCK_BASE"]:
        val = os.environ.get(var)
        if val:
            search_paths.append(Path(val) / "parameters")
    try:
        dock6_path = subprocess.run(
            ["which", "dock6"], capture_output=True, text=True
        ).stdout.strip()
        if dock6_path:
            search_paths.append(Path(dock6_path).parent.parent / "parameters")
    except Exception:
        pass

    result = {}
    for key, filename in param_files.items():
        result[key] = filename
        for search_dir in search_paths:
            candidate = search_dir / filename
            if candidate.exists():
                result[key] = str(candidate)
                break
    return result


def resolve_grid_prefix(grid_dir: str, energy_grid: str) -> str:
    """Resolve grid prefix for DOCK6 (path without .nrg extension)."""
    prefix = energy_grid.replace(".nrg", "").replace(".NRG", "")
    return str(Path(grid_dir) / prefix)
''')

    # --- 02a score_collector ---
    write_file("01_src/molecular_docking/m02_collection/score_collector.py", '''\
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
''')


# =============================================================================
# 5. CLI SCRIPTS
# =============================================================================

def create_cli_scripts():
    """Create CLI scripts for each module."""
    print("\n=== Creating CLI scripts ===")

    # --- 00a ---
    write_file("02_scripts/00a_molecule_parser.py", '''\
#!/usr/bin/env python3
"""
00a Molecule Parser - CLI
Usage:
    python 02_scripts/00a_molecule_parser.py \\
        --config 03_configs/00a_molecule_parser.yaml \\
        --campaign 04_data/campaigns/phermit_groove/campaign_config.yaml
"""
import argparse, logging, sys, yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))
from molecular_docking.m00_preparation.molecule_parser import run_molecule_parser

logger = logging.getLogger(__name__)

def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f: return yaml.safe_load(f)

def main():
    parser = argparse.ArgumentParser(description="Parse molecules")
    parser.add_argument("--config", type=str, default=None)
    parser.add_argument("--campaign", type=str, default=None)
    parser.add_argument("--input", type=str, default=None)
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args()

    # Defaults
    params = {"conformer_strategy": "best_rmsd", "compute_properties": True,
              "detect_smiles_duplicates": True}
    name_source, name_property, name_column, smiles_column = "header_first_word", "_Name", None, None
    input_file, output_dir, campaign_id = None, None, "direct"

    if args.campaign:
        cc = load_yaml(args.campaign)
        campaign_dir = Path(args.campaign).parent
        campaign_id = cc.get("campaign_id", campaign_dir.name)
        mc = cc.get("molecules", {})
        input_file = str(campaign_dir / mc.get("input_file", ""))
        name_column = mc.get("name_column")
        smiles_column = mc.get("smiles_column")
        name_source = mc.get("name_source", name_source)
        name_property = mc.get("name_property", name_property)
        output_dir = str(Path("05_results") / campaign_id / "00a_molecule_parser")

    if args.config:
        mc = load_yaml(args.config).get("parameters", {})
        params.update(mc)

    if args.input: input_file = args.input
    if args.output: output_dir = args.output
    if not input_file: parser.error("Provide --campaign or --input")
    if not output_dir: output_dir = "05_results/00a_molecule_parser"

    logger.info("=" * 60)
    logger.info("  MOLECULAR_DOCKING - Module 00a: Molecule Parser")
    logger.info("=" * 60)

    result = run_molecule_parser(
        input_file=input_file, output_dir=output_dir,
        name_source=name_source, name_property=name_property,
        name_column=name_column, smiles_column=smiles_column,
        **params,
    )
    return 0

if __name__ == "__main__":
    sys.exit(main())
''')

    # --- 00b ---
    write_file("02_scripts/00b_receptor_preparation.py", '''\
#!/usr/bin/env python3
"""
00b Receptor Preparation - CLI
OPTIONAL: skips if receptor.protonation.enabled=false.
"""
import argparse, logging, sys, yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))
from molecular_docking.m00_preparation.receptor_preparation import run_receptor_preparation, validate_prepared_mol2

logger = logging.getLogger(__name__)

def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f: return yaml.safe_load(f)

def main():
    parser = argparse.ArgumentParser(description="Prepare receptor for DOCK6")
    parser.add_argument("--config", type=str, default=None)
    parser.add_argument("--campaign", type=str, default=None)
    parser.add_argument("--ph", type=float, default=None)
    args = parser.parse_args()

    logger.info("=" * 60)
    logger.info("  MOLECULAR_DOCKING - Module 00b: Receptor Preparation")
    logger.info("=" * 60)

    if args.campaign:
        cc = load_yaml(args.campaign)
        rc = cc.get("receptor", {})
        pc = rc.get("protonation", {})
        if not pc.get("enabled", False):
            prepared = rc.get("prepared_mol2")
            if prepared:
                logger.info(f"Protonation DISABLED. Using: {prepared}")
            else:
                logger.info("Protonation DISABLED, no prepared_mol2. SKIPPED.")
            return 0
    logger.info("Receptor protonation — to be implemented")
    return 0

if __name__ == "__main__":
    sys.exit(main())
''')

    # --- 00c ---
    write_file("02_scripts/00c_ligand_preparation.py", '''\
#!/usr/bin/env python3
"""00c Ligand Preparation - CLI"""
import argparse, logging, sys, yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))
from molecular_docking.m00_preparation.ligand_preparation import run_ligand_preparation

logger = logging.getLogger(__name__)

def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f: return yaml.safe_load(f)

def main():
    parser = argparse.ArgumentParser(description="Prepare ligands for DOCK6")
    parser.add_argument("--config", type=str, default=None)
    parser.add_argument("--campaign", type=str, default=None)
    parser.add_argument("--input-csv", type=str, default=None)
    parser.add_argument("--output", type=str, default=None)
    parser.add_argument("--ph", type=float, default=None)
    args = parser.parse_args()

    docking_ph, protonate_enabled, protonate_tool = 7.2, True, "dimorphite_dl"
    charge_method, atom_type = "bcc", "gaff2"
    molecules_csv, output_dir, campaign_id = None, None, "direct"

    if args.campaign:
        cc = load_yaml(args.campaign)
        campaign_dir = Path(args.campaign).parent
        campaign_id = cc.get("campaign_id", campaign_dir.name)
        docking_ph = cc.get("docking_ph", docking_ph)
        mc = cc.get("molecules", {}).get("protonation", {})
        protonate_enabled = mc.get("enabled", True)
        protonate_tool = mc.get("tool", protonate_tool)
        molecules_csv = str(Path("05_results") / campaign_id / "00a_molecule_parser" / "unique_molecules.csv")
        output_dir = str(Path("05_results") / campaign_id / "00c_ligand_preparation")

    if args.config:
        p = load_yaml(args.config).get("parameters", {})
        charge_method = p.get("charge_method", charge_method)
        atom_type = p.get("atom_type", atom_type)

    if args.input_csv: molecules_csv = args.input_csv
    if args.output: output_dir = args.output
    if args.ph: docking_ph = args.ph
    if not molecules_csv: parser.error("Provide --campaign or --input-csv")

    logger.info("=" * 60)
    logger.info("  MOLECULAR_DOCKING - Module 00c: Ligand Preparation")
    logger.info("=" * 60)

    result = run_ligand_preparation(
        molecules_csv=molecules_csv, output_dir=output_dir,
        docking_ph=docking_ph, protonate_enabled=protonate_enabled,
        protonate_tool=protonate_tool, charge_method=charge_method,
        atom_type=atom_type,
    )
    return 0

if __name__ == "__main__":
    sys.exit(main())
''')

    # --- 01a ---
    write_file("02_scripts/01a_grid_generation.py", '''\
#!/usr/bin/env python3
"""01a Grid Generation - CLI. OPTIONAL: skips if grids exist."""
import argparse, logging, sys, yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))
from molecular_docking.m01_docking.grid_generation import validate_existing_grids

logger = logging.getLogger(__name__)

def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f: return yaml.safe_load(f)

def main():
    parser = argparse.ArgumentParser(description="Generate DOCK6 grids")
    parser.add_argument("--config", type=str, default=None)
    parser.add_argument("--campaign", type=str, default=None)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    logger.info("=" * 60)
    logger.info("  MOLECULAR_DOCKING - Module 01a: Grid Generation")
    logger.info("=" * 60)

    if args.campaign:
        cc = load_yaml(args.campaign)
        campaign_dir = Path(args.campaign).parent
        gc = cc.get("grids", {})
        grid_dir = gc.get("grid_dir", "grids/")
        grid_path = Path(grid_dir) if Path(grid_dir).is_absolute() else campaign_dir / grid_dir

        if not args.force and validate_existing_grids(
            str(grid_path), gc.get("spheres_file", "selected_spheres.sph"),
            gc.get("energy_grid", "grid.nrg"), gc.get("bump_grid", "grid.bmp"),
        ):
            logger.info(f"Valid grids at: {grid_path}")
            logger.info("SKIPPED (use --force to regenerate)")
            return 0

    logger.info("Grid generation — to be implemented. Use pre-existing grids.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
''')

    # --- 01b ---
    write_file("02_scripts/01b_dock6_run.py", '''\
#!/usr/bin/env python3
"""01b DOCK6 Run - CLI"""
import argparse, logging, sys, yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))
from molecular_docking.m01_docking.dock6_runner import run_dock6_batch, resolve_grid_prefix

logger = logging.getLogger(__name__)

def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f: return yaml.safe_load(f)

def main():
    parser = argparse.ArgumentParser(description="Run DOCK6 docking")
    parser.add_argument("--config", type=str, required=True)
    parser.add_argument("--campaign", type=str, required=True)
    parser.add_argument("--name", type=str, default=None, help="Dock only this molecule")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--timeout", type=int, default=None)
    args = parser.parse_args()

    cc = load_yaml(args.campaign)
    campaign_dir = Path(args.campaign).parent
    campaign_id = cc.get("campaign_id", campaign_dir.name)

    gc = cc.get("grids", {})
    grid_dir = gc.get("grid_dir", "grids/")
    grid_path = Path(grid_dir) if Path(grid_dir).is_absolute() else campaign_dir / grid_dir
    spheres = str(grid_path / gc.get("spheres_file", "selected_spheres.sph"))
    grid_prefix = resolve_grid_prefix(str(grid_path), gc.get("energy_grid", "grid.nrg"))

    ligand_dir = str(Path("05_results") / campaign_id / "00c_ligand_preparation" / "mol2")
    output_dir = str(Path("05_results") / campaign_id / "01b_dock6_run")

    params = load_yaml(args.config).get("parameters", {})
    molecule_filter = [args.name] if args.name else None

    logger.info("=" * 60)
    logger.info("  MOLECULAR_DOCKING - Module 01b: DOCK6 Run")
    logger.info("=" * 60)

    result = run_dock6_batch(
        ligand_mol2_dir=ligand_dir, spheres_file=spheres,
        grid_prefix=grid_prefix, output_dir=output_dir,
        search_method=params.get("search_method", "flex"),
        max_orientations=params.get("max_orientations", 1000),
        num_scored_conformers=params.get("num_scored_conformers", 20),
        minimize=params.get("minimize", True),
        simplex_max_iterations=params.get("simplex_max_iterations", 500),
        timeout_per_molecule=args.timeout or params.get("timeout_per_molecule", 600),
        molecule_filter=molecule_filter, dry_run=args.dry_run,
    )
    return 0

if __name__ == "__main__":
    sys.exit(main())
''')

    # --- 02a ---
    write_file("02_scripts/02a_score_collection.py", '''\
#!/usr/bin/env python3
"""02a Score Collection - CLI"""
import argparse, logging, sys, yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))
from molecular_docking.m02_collection.score_collector import run_score_collection

logger = logging.getLogger(__name__)

def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f: return yaml.safe_load(f)

def main():
    parser = argparse.ArgumentParser(description="Collect scores → dock2profile Excel")
    parser.add_argument("--config", type=str, required=True)
    parser.add_argument("--campaign", type=str, required=True)
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args()

    cc = load_yaml(args.campaign)
    campaign_id = cc.get("campaign_id", Path(args.campaign).parent.name)
    oc = cc.get("output", {})

    docking_dir = str(Path("05_results") / campaign_id / "01b_dock6_run")
    molecules_csv = str(Path("05_results") / campaign_id / "00a_molecule_parser" / "unique_molecules.csv")
    output_dir = args.output or str(Path("05_results") / campaign_id / "02a_score_collection")
    params = load_yaml(args.config).get("parameters", {})

    logger.info("=" * 60)
    logger.info("  MOLECULAR_DOCKING - Module 02a: Score Collection")
    logger.info("=" * 60)

    result = run_score_collection(
        docking_dir=docking_dir, molecules_csv=molecules_csv,
        output_dir=output_dir,
        score_key=params.get("score_key", "Grid_Score"),
        max_molecules=params.get("max_molecules", 500),
        scores_filename=oc.get("scores_filename", "01_top_500_molecules.xlsx"),
        mol2_dirname=oc.get("mol2_dirname", "docked_molecules"),
        source_label=cc.get("metadata", {}).get("source"),
    )
    return 0

if __name__ == "__main__":
    sys.exit(main())
''')


# =============================================================================
# 6. CONFIG YAMLS
# =============================================================================

def create_config_yamls():
    """Create module configuration YAML files."""
    print("\n=== Creating config YAMLs ===")

    write_file("03_configs/00a_molecule_parser.yaml", '''\
# 00a Molecule Parser
parameters:
  conformer_strategy: "best_rmsd"   # best_rmsd | first | all
  compute_properties: true
  detect_smiles_duplicates: true
  log_level: "INFO"
outputs:
  subdir: "00a_molecule_parser"
''')

    write_file("03_configs/00b_receptor_preparation.yaml", '''\
# 00b Receptor Preparation (OPTIONAL)
parameters:
  protonation_tool: null              # null = use campaign_config
  force_field: "AMBER"
  remove_alt_conformations: true
  log_level: "INFO"
outputs:
  subdir: "00b_receptor_preparation"
''')

    write_file("03_configs/00c_ligand_preparation.yaml", '''\
# 00c Ligand Preparation
parameters:
  charge_method: "bcc"                # bcc (AM1-BCC) | gas (Gasteiger)
  atom_type: "gaff2"
  antechamber_timeout: 300
  obabel_fallback: true
  n_3d_attempts: 5
  log_level: "INFO"
outputs:
  subdir: "00c_ligand_preparation"
''')

    write_file("03_configs/01a_grid_generation.yaml", '''\
# 01a Grid Generation (OPTIONAL - skip if grids exist)
parameters:
  probe_radius: 1.4
  max_spheres: 50
  box_margin: 10.0
  grid_spacing: 0.3
  energy_cutoff_distance: 9999.0
  bump_overlap: 0.75
  log_level: "INFO"
outputs:
  subdir: "01a_grid_generation"
''')

    write_file("03_configs/01b_dock6_run.yaml", '''\
# 01b DOCK6 Run
parameters:
  search_method: "flex"               # flex | rigid
  min_anchor_size: 5
  pruning_max_orients: 1000
  pruning_clustering_cutoff: 100
  pruning_conformer_score_cutoff: 100.0
  max_orientations: 1000
  minimize: true
  simplex_max_iterations: 500
  simplex_max_cycles: 1
  simplex_score_converge: 0.1
  simplex_cycle_converge: 1.0
  simplex_trans_step: 1.0
  simplex_rot_step: 0.1
  simplex_tors_step: 10.0
  num_scored_conformers: 20
  timeout_per_molecule: 600
  write_orientations: false
  log_level: "INFO"
outputs:
  subdir: "01b_dock6_run"
''')

    write_file("03_configs/02a_score_collection.yaml", '''\
# 02a Score Collection
parameters:
  score_key: "Grid_Score"
  max_molecules: 500                  # 0 = all
  extract_best_pose_mol2: true
  keep_all_poses: false
  compute_properties: true
  log_level: "INFO"
outputs:
  subdir: "02a_score_collection"
''')


# =============================================================================
# 7. CAMPAIGN CONFIG TEMPLATE
# =============================================================================

def create_campaign_template():
    """Create example campaign_config.yaml."""
    print("\n=== Creating campaign template ===")

    write_file("04_data/campaigns/example_campaign/campaign_config.yaml", '''\
# =============================================================================
# Campaign Configuration Template
# =============================================================================
# One CAMPAIGN = one receptor + one set of molecules + docking conditions.
# Copy to 04_data/campaigns/{your_campaign}/ and edit.
# =============================================================================

campaign_id: "example_campaign"
description: "Example docking campaign"

receptor:
  pdb: "receptor/receptor.pdb"
  chain: null
  protonation:
    enabled: false                    # true = protonate with tool below
    tool: "chimera"                   # chimera | pdb2pqr | reduce | obabel
    force_field: "AMBER"
  prepared_mol2: null                 # Pre-protonated mol2 (skip 00b)
  remove_water: true
  remove_hetatm: true

docking_ph: 7.2                       # Affects protonation of receptor & ligands

molecules:
  input_file: "molecules/input_molecules.sdf"
  name_column: "Name"
  smiles_column: "Smile"
  name_source: "header_first_word"
  name_property: "_Name"
  protonation:
    enabled: true
    tool: "dimorphite_dl"
    obabel_fallback: true

grids:
  grid_dir: "grids/"
  spheres_file: "selected_spheres.sph"
  energy_grid: "grid.nrg"
  bump_grid: "grid.bmp"
  generate: false

output:
  scores_filename: "01_top_500_molecules.xlsx"
  mol2_dirname: "docked_molecules"

metadata:
  source: null                        # "Phermit" | "HTS1710" | etc.
  date: null
  notes: null
''')


# =============================================================================
# 8. TESTS & DOCS
# =============================================================================

def create_tests_and_docs():
    """Create test stubs and documentation."""
    print("\n=== Creating tests & docs ===")

    write_file("tests/test_pipeline.py", '''\
"""Tests for molecular_docking pipeline."""
import pytest

class TestMoleculeParser:
    def test_import(self):
        from molecular_docking.m00_preparation.molecule_parser import run_molecule_parser
        assert callable(run_molecule_parser)

class TestLigandPreparation:
    def test_import(self):
        from molecular_docking.m00_preparation.ligand_preparation import run_ligand_preparation
        assert callable(run_ligand_preparation)

class TestDock6Runner:
    def test_import(self):
        from molecular_docking.m01_docking.dock6_runner import run_dock6_batch
        assert callable(run_dock6_batch)

    def test_resolve_grid_prefix(self):
        from molecular_docking.m01_docking.dock6_runner import resolve_grid_prefix
        result = resolve_grid_prefix("/path/to/grids", "grid.nrg")
        assert result == "/path/to/grids/grid"

class TestGridValidation:
    def test_missing_grids(self, tmp_path):
        from molecular_docking.m01_docking.grid_generation import validate_existing_grids
        assert validate_existing_grids(str(tmp_path), "s.sph", "g.nrg", "g.bmp") is False

class TestScoreCollector:
    def test_import(self):
        from molecular_docking.m02_collection.score_collector import run_score_collection
        assert callable(run_score_collection)
''')

    write_file("README.md", '''\
# MOLECULAR_DOCKING

Generic DOCK6 molecular docking pipeline. Produces outputs compatible with **dock2profile**.

## Setup

```bash
conda env create -f environment.yaml
conda activate molecular_docking_env
pip install -e ".[dev]"
bash check_dependencies.sh
```

## Usage

```bash
# Create campaign
cp -r 04_data/campaigns/example_campaign 04_data/campaigns/my_campaign
# Edit campaign_config.yaml

# Run pipeline
python 02_scripts/00a_molecule_parser.py      --config 03_configs/00a_molecule_parser.yaml      --campaign 04_data/campaigns/my_campaign/campaign_config.yaml
python 02_scripts/00b_receptor_preparation.py --config 03_configs/00b_receptor_preparation.yaml --campaign 04_data/campaigns/my_campaign/campaign_config.yaml
python 02_scripts/00c_ligand_preparation.py   --config 03_configs/00c_ligand_preparation.yaml   --campaign 04_data/campaigns/my_campaign/campaign_config.yaml
python 02_scripts/01a_grid_generation.py      --config 03_configs/01a_grid_generation.yaml      --campaign 04_data/campaigns/my_campaign/campaign_config.yaml
python 02_scripts/01b_dock6_run.py            --config 03_configs/01b_dock6_run.yaml            --campaign 04_data/campaigns/my_campaign/campaign_config.yaml
python 02_scripts/02a_score_collection.py     --config 03_configs/02a_score_collection.yaml     --campaign 04_data/campaigns/my_campaign/campaign_config.yaml
```

## Structure

```
01_src/    Core modules (logic, no CLI)
02_scripts/ CLI scripts (argparse + YAML → calls core)
03_configs/ YAML per module (algorithmic params)
04_data/    Campaigns (receptor + molecules + grids)
05_results/ Outputs per campaign/module
```
''')


# =============================================================================
# 9. DEPENDENCY CHECKER
# =============================================================================

def create_check_script():
    """Create check_dependencies.sh."""
    print("\n=== Creating dependency checker ===")

    write_file("check_dependencies.sh", '''\
#!/bin/bash
echo "============================================================"
echo "  MOLECULAR_DOCKING - Dependency Check"
echo "============================================================"
echo ""
PASS=0; FAIL=0; WARN=0
check_cmd() {
    local name="$1" cmd="$2" req="$3"
    if eval "$cmd" > /dev/null 2>&1; then
        echo "  [OK]   $name"; ((PASS++))
    elif [ "$req" = "required" ]; then
        echo "  [FAIL] $name: NOT FOUND"; ((FAIL++))
    else
        echo "  [WARN] $name: not found (optional)"; ((WARN++))
    fi
}
echo "--- Python ---"
check_cmd "python3" "python3 --version" "required"
echo ""
echo "--- Python Libraries ---"
for lib in rdkit pandas numpy openpyxl yaml; do
    if python3 -c "import $lib" 2>/dev/null; then
        echo "  [OK]   $lib"; ((PASS++))
    else
        echo "  [FAIL] $lib"; ((FAIL++))
    fi
done
for lib in dimorphite_dl pdb2pqr openbabel; do
    if python3 -c "import $lib" 2>/dev/null; then
        echo "  [OK]   $lib (optional)"; ((PASS++))
    else
        echo "  [WARN] $lib (optional)"; ((WARN++))
    fi
done
echo ""
echo "--- DOCK6 ---"
check_cmd "dock6" "which dock6" "required"
check_cmd "grid" "which grid" "optional"
check_cmd "sphgen" "which sphgen" "optional"
check_cmd "sphere_selector" "which sphere_selector" "optional"
check_cmd "showbox" "which showbox" "optional"
echo ""
echo "--- AmberTools ---"
check_cmd "antechamber" "which antechamber" "required"
check_cmd "parmchk2" "which parmchk2" "optional"
check_cmd "tleap" "which tleap" "optional"
check_cmd "reduce" "which reduce" "optional"
echo ""
echo "--- OpenBabel ---"
check_cmd "obabel" "obabel -V" "required"
echo ""
echo "--- Optional ---"
check_cmd "chimera" "which chimera" "optional"
echo ""
echo "============================================================"
echo "  SUMMARY: $PASS passed, $FAIL failed, $WARN warnings"
echo "============================================================"
[ $FAIL -gt 0 ] && echo "  Fix FAIL items before running pipeline" && exit 1
echo "  Ready!" && exit 0
''')

    # Make executable
    os.chmod(ROOT / "check_dependencies.sh", 0o755)


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 60)
    print("  MOLECULAR_DOCKING - Project Setup v" + VERSION)
    print("=" * 60)
    print(f"  Root: {ROOT}")

    create_directories()
    create_packaging_files()
    create_init_files()
    create_core_modules()
    create_cli_scripts()
    create_config_yamls()
    create_campaign_template()
    create_tests_and_docs()
    create_check_script()

    # Count files created
    all_files = list(ROOT.rglob("*"))
    n_files = sum(1 for f in all_files if f.is_file() and f.name != "setup_project.py")

    print(f"\n{'=' * 60}")
    print(f"  DONE: {n_files} files created")
    print(f"{'=' * 60}")
    print(f"""
Next steps:
  1. conda activate molecular_docking_env
  2. conda env update -f environment.yaml --prune
  3. pip install -e ".[dev]"
  4. bash check_dependencies.sh
  5. pytest tests/ -v
""")


if __name__ == "__main__":
    main()

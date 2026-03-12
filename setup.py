#!/usr/bin/env python3
"""
molecular_docking - Setup Script
===================================
Backwards-compatible setup.py. Metadata also lives in pyproject.toml;
this file exists for `pip install -e .` compatibility and for
environments that don't support PEP 621 yet.

Installation:
    # Development (editable):
    pip install -e ".[dev]"

    # Production:
    pip install .

    # With protonation tools:
    pip install -e ".[protonation]"

NOTE:
    RDKit, OpenBabel, and AmberTools should be installed via conda
    BEFORE running pip install. See environment.yaml.

    DOCK6 requires a separate academic license from UCSF.
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read version from __init__.py
version = "1.0.0"
init_file = Path(__file__).parent / "01_src" / "molecular_docking" / "__init__.py"
if init_file.exists():
    for line in init_file.read_text().splitlines():
        if line.startswith("__version__"):
            version = line.split('"')[1]
            break

# Read README for long description
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

    # --- Package discovery ---
    package_dir={"": "01_src"},
    packages=find_packages(where="01_src"),

    # --- Core dependencies (pip-installable only) ---
    # RDKit, OpenBabel, AmberTools → install via conda (see environment.yaml)
    install_requires=[
        "pandas>=1.5.0",
        "numpy>=1.23.0",
        "openpyxl>=3.0.0",
        "pyyaml>=6.0",
    ],

    # --- Optional dependency groups ---
    extras_require={
        # Protonation tools (pip-installable)
        "protonation": [
            "dimorphite_dl>=1.3.2",
            "pdb2pqr>=3.6.0",
        ],
        # Testing
        "dev": [
            "pytest>=7.0",
            "pytest-cov>=4.0",
        ],
        # Full install (everything pip can handle)
        "all": [
            "dimorphite_dl>=1.3.2",
            "pdb2pqr>=3.6.0",
            "pytest>=7.0",
            "pytest-cov>=4.0",
        ],
    },

    # --- Metadata ---
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Microsoft :: Windows",
    ],
    keywords=[
        "molecular-docking",
        "dock6",
        "drug-discovery",
        "virtual-screening",
        "computational-chemistry",
    ],
)

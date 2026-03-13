#!/usr/bin/env python3
"""
Backwards-compatible setup.py. Metadata lives in pyproject.toml.

    pip install -e ".[dev]"
"""
from setuptools import setup, find_packages
from pathlib import Path

version = "2.0.0"
init_file = Path(__file__).parent / "01_src" / "molecular_docking" / "__init__.py"
if init_file.exists():
    for line in init_file.read_text().splitlines():
        if line.startswith("__version__"):
            version = line.split('"')[1]
            break

setup(
    name="molecular_docking",
    version=version,
    description="DOCK6 molecular docking pipeline compatible with dock2profile",
    author="Benjamin Prieto",
    license="MIT",
    python_requires=">=3.9,<3.13",
    package_dir={"": "01_src"},
    packages=find_packages(where="01_src"),
    install_requires=[
        "pandas>=1.5.0",
        "numpy>=1.23.0",
        "openpyxl>=3.0.0",
        "pyyaml>=6.0",
        "tqdm>=4.64.0",
    ],
    extras_require={
        "protonation": ["pdb2pqr>=3.6.0"],
        "dev": ["pytest>=7.0", "pytest-cov>=4.0"],
    },
)
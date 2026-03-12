#!/usr/bin/env python3
"""
Backwards-compatible setup.py. Metadata lives in pyproject.toml.

    pip install -e ".[dev]"
"""
from setuptools import setup, find_packages
from pathlib import Path

version = "1.0.0"
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

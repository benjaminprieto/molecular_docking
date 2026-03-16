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
        from molecular_docking.m01_docking.score_collector import run_score_collection
        assert callable(run_score_collection)

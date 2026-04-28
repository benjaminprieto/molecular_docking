"""Tests for molecular_docking pipeline."""
import pytest


class TestMoleculeParser:
    def test_import(self):
        from molecular_docking.m00_preparation.molecule_parser import run_molecule_parser
        assert callable(run_molecule_parser)


class TestBindingSiteDefinition:
    def test_import(self):
        from molecular_docking.m00_preparation.binding_site_definition import run_binding_site_definition
        assert callable(run_binding_site_definition)


class TestAntechamberPreparation:
    def test_import(self):
        from molecular_docking.m01_docking.antechamber_preparation import run_antechamber_preparation
        assert callable(run_antechamber_preparation)


class TestGridGeneration:
    def test_import(self):
        from molecular_docking.m01_docking.grid_generation import validate_existing_grids
        assert callable(validate_existing_grids)

    def test_missing_grids(self, tmp_path):
        from molecular_docking.m01_docking.grid_generation import validate_existing_grids
        assert validate_existing_grids(str(tmp_path), "s.sph", "g.nrg", "g.bmp") is False


class TestDock6Runner:
    def test_import(self):
        from molecular_docking.m01_docking.dock6_runner import run_dock6_batch
        assert callable(run_dock6_batch)

    def test_resolve_grid_prefix(self):
        from molecular_docking.m01_docking.dock6_runner import resolve_grid_prefix
        result = resolve_grid_prefix("/path/to/grids", "grid.nrg")
        assert result == "/path/to/grids/grid"


class TestScoreCollector:
    def test_import(self):
        from molecular_docking.m01_docking.score_collector import run_score_collection
        assert callable(run_score_collection)


class TestGninaPreparation:
    def test_import(self):
        from molecular_docking.m02_gnina.gnina_preparation import run_gnina_preparation
        assert callable(run_gnina_preparation)

    def test_binding_box(self):
        from molecular_docking.m02_gnina.gnina_preparation import GninaBindingBox
        box = GninaBindingBox(center_x=10.0, center_y=20.0, center_z=30.0,
                              size_x=25.0, size_y=25.0, size_z=25.0)
        assert box.volume == 25.0 * 25.0 * 25.0
        args = box.to_gnina_args()
        assert "--center_x" in args
        assert "10.000" in args


class TestGninaRunner:
    def test_import(self):
        from molecular_docking.m02_gnina.gnina_runner import run_gnina_docking
        assert callable(run_gnina_docking)

    def test_parse_gnina_output(self):
        from molecular_docking.m02_gnina.gnina_runner import parse_gnina_output
        sample = """
mode |   affinity | CNN score | CNN affinity
-----+------------+-----------+--------------
   1       -8.320      0.847       -7.920
   2       -7.952      0.721       -7.450
   3       -7.128      0.612       -6.880
"""
        poses = parse_gnina_output(sample)
        assert len(poses) == 3
        assert poses[0]["vina_affinity"] == -8.320
        assert poses[0]["cnn_score"] == 0.847
        assert poses[2]["mode"] == 3


class TestGninaScoreCollector:
    def test_import(self):
        from molecular_docking.m02_gnina.gnina_score_collector import run_gnina_score_collection
        assert callable(run_gnina_score_collection)
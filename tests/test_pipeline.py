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


class TestVinaPreparation:
    def test_import(self):
        from molecular_docking.m02_gnina.vina_preparation import run_vina_preparation
        assert callable(run_vina_preparation)

    def test_binding_box(self):
        from molecular_docking.m02_gnina.vina_preparation import VinaBindingBox
        box = VinaBindingBox(center_x=10.0, center_y=20.0, center_z=30.0,
                             size_x=25.0, size_y=25.0, size_z=25.0)
        assert box.volume == 25.0 * 25.0 * 25.0
        args = box.to_vina_args()
        assert "--center_x" in args
        assert "10.000" in args


class TestVinaRunner:
    def test_import(self):
        from molecular_docking.m02_gnina.vina_runner import run_vina_docking
        assert callable(run_vina_docking)

    def test_parse_vina_output(self):
        from molecular_docking.m02_gnina.vina_runner import parse_vina_output
        sample = """
Scoring function : vina
Detected 4 CPUs
mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1       -8.320      0.000      0.000
   2       -7.952      1.834      2.456
   3       -7.128      2.103      4.221
"""
        poses = parse_vina_output(sample)
        assert len(poses) == 3
        assert poses[0]["affinity"] == -8.320
        assert poses[1]["rmsd_lb"] == 1.834
        assert poses[2]["mode"] == 3


class TestVinaScoreCollector:
    def test_import(self):
        from molecular_docking.m02_gnina.vina_score_collector import run_vina_score_collection
        assert callable(run_vina_score_collection)
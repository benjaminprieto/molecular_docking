#!/bin/bash
# =============================================================================
# Test footprint re-scoring on UDX
# =============================================================================
# This re-scores existing poses (no re-docking) with per-residue decomposition
# =============================================================================

cd /tmp
mkdir -p footprint_test
cd footprint_test

# Symlinks (80-char path fix)
ln -sf /home/bprieto/projects/molecular_docking/05_results/UDX_pharmit_pH63/01c_dock6_run/04_data_campaigns_example_campaign_molecules_UDX.pdb/04_data_campaigns_example_campaign_molecules_UDX.pdb_scored.mol2 poses.mol2
ln -sf /home/bprieto/projects/molecular_docking/04_data/campaigns/UDX_pharmit_pH63/reference/UDX.mol2 reference.mol2
ln -sf /home/bprieto/projects/molecular_docking/05_results/UDX_pharmit_pH63/00b_receptor_preparation/rec_charged.mol2 receptor.mol2
ln -sf /opt/dock6/parameters/vdw_AMBER_parm99.defn vdw.defn
ln -sf /opt/dock6/parameters/flex.defn flex.defn
ln -sf /opt/dock6/parameters/flex_drive.tbl flex_drive.tbl

# Write dock.in for footprint re-scoring
cat > dock6_fps.in << 'DOCK_IN'
conformer_search_type                            rigid
use_internal_energy                              no
ligand_atom_file                                 poses.mol2
limit_max_ligands                                no
skip_molecule                                    no
read_mol_solvation                               no
calculate_rmsd                                   no
use_database_filter                              no
orient_ligand                                    no
bump_filter                                      no
score_molecules                                  yes
contact_score_primary                            no
contact_score_secondary                          no
grid_score_primary                               no
grid_score_secondary                             no
multigrid_score_primary                          no
multigrid_score_secondary                        no
dock3.5_score_primary                            no
dock3.5_score_secondary                          no
continuous_score_primary                         no
continuous_score_secondary                       no
footprint_similarity_score_primary               yes
footprint_similarity_score_secondary             no
fps_score_use_footprint_reference_mol2           yes
fps_score_footprint_reference_mol2_filename      reference.mol2
fps_score_foot_compare_type                      Euclidean
fps_score_normalize_foot                         no
fps_score_foot_comp_all_residue                  yes
fps_score_receptor_filename                      receptor.mol2
fps_score_vdw_att_exp                            6
fps_score_vdw_rep_exp                            9
fps_score_vdw_rep_rad_scale                      1
fps_score_use_distance_dependent_dielectric      yes
fps_score_dielectric                             4.0
fps_score_vdw_fp_scale                           1
fps_score_es_fp_scale                            1
fps_score_hb_fp_scale                            0
pharmacophore_score_secondary                    no
descriptor_score_secondary                       no
gbsa_zou_score_secondary                         no
gbsa_hawkins_score_secondary                     no
SASA_score_secondary                             no
amber_score_secondary                            no
minimize_ligand                                  no
atom_model                                       all
vdw_defn_file                                    vdw.defn
flex_defn_file                                   flex.defn
flex_drive_file                                  flex_drive.tbl
ligand_outfile_prefix                            UDX_footprint
write_orientations                               no
num_scored_conformers                            100
rank_ligands                                     no
DOCK_IN

echo "Running footprint re-scoring..."
dock6 -i dock6_fps.in -o dock6_fps.out

echo ""
echo "=== Checking output ==="
ls -la UDX_footprint_scored.mol2 2>/dev/null && echo "Output exists" || echo "NO OUTPUT"

echo ""
echo "=== Footprint data in output ==="
grep -i "fps\|footprint\|FPS\|vdw_fp\|es_fp" UDX_footprint_scored.mol2 2>/dev/null | head -20

echo ""
echo "=== Header of first pose ==="
head -30 UDX_footprint_scored.mol2 2>/dev/null

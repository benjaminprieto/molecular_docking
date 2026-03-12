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

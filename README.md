# MOLECULAR_DOCKING

Pipeline de docking molecular con DOCK6. Produce outputs compatibles con **dock2profile**.

## Requisitos

| Dependencia | Instalación | Notas |
|---|---|---|
| Python 3.9-3.12 | conda | |
| RDKit, OpenBabel, AmberTools | conda | Via `environment.yaml` |
| PDB2PQR | pip | Protonación pH-aware del receptor |
| **DOCK6** | Manual | Licencia académica gratuita de [UCSF](https://dock.compbio.ucsf.edu/DOCK_6/index.htm) |
| **ChimeraX** | Manual | [Descargar](https://www.cgl.ucsf.edu/chimerax/download.html) — necesario para preparar el receptor |

## Instalación

```bash
# 1. Clonar
git clone https://github.com/benjaminprieto/molecular_docking.git
cd molecular_docking

# 2. Crear entorno conda
conda env create -f environment.yaml
conda activate molecular_docking_env

# 3. Instalar paquete en modo editable
pip install -e ".[dev]"

# 4. Verificar dependencias
bash check_dependencies.sh
```

### DOCK6

DOCK6 requiere licencia académica gratuita. Después de obtenerla:

```bash
# Compilar e instalar en /opt/dock6/
tar -xzf dock.6.X.tar.gz
cd dock6
./configure gnu
make

# Agregar al PATH (en ~/.bashrc)
export PATH=/opt/dock6/bin:$PATH
```

### ChimeraX

Descargar desde https://www.cgl.ucsf.edu/chimerax/download.html e instalar. El pipeline busca el binario en `/usr/bin/chimerax-daily`, luego `chimerax` en PATH.

## Estructura del proyecto

```
molecular_docking/
├── 01_src/molecular_docking/       Core modules (lógica, sin CLI)
│   ├── m00_preparation/
│   │   ├── molecule_parser.py          00a — parseo de moléculas
│   │   ├── receptor_preparation.py     00b — receptor → mol2 DOCK6-ready
│   │   ├── ionization_profiling.py     00c — protonación de ligandos al pH
│   │   ├── antechamber_preparation.py  00d — mol2 con cargas AM1-BCC
│   │   └── binding_site_definition.py  00e — recorte del receptor (opcional)
│   ├── m01_docking/
│   │   ├── grid_generation.py          01a — DMS → spheres → grids
│   │   └── dock6_runner.py             01b — dock6 por molécula
│   └── m02_collection/
│       └── score_collector.py          02a — parseo de scores → Excel
├── 02_scripts/                     CLI scripts (argparse + YAML → core)
├── 03_configs/                     YAML por módulo (parámetros algorítmicos)
├── 04_data/campaigns/              Campañas (receptor + moléculas + grids)
│   └── example_campaign/
│       ├── campaign_config.yaml        Fuente de verdad de la campaña
│       ├── receptor/                   PDB del receptor
│       └── molecules/                  Moléculas de entrada
├── 05_results/                     Outputs por campaña/módulo
├── environment.yaml
├── pyproject.toml
└── check_dependencies.sh
```

## Flujo del pipeline

```
00a molecule_parser       → unique_molecules.csv + .sdf
00b receptor_preparation  → rec_charged.mol2 + rec_noH.pdb
00c ionization_profiling  → SDF protonados por pH
00d antechamber           → mol2 con AM1-BCC charges
00e binding_site_def      → rec_noH_site.pdb (opcional)
01a grid_generation       → DMS, spheres, box, grid.nrg/bmp
01b dock6_run             → scored mol2 por molécula
02a score_collection      → Excel compatible dock2profile
```

## Uso

### 1. Crear campaña

```bash
cp -r 04_data/campaigns/example_campaign 04_data/campaigns/mi_campana
```

Editar `campaign_config.yaml` con el receptor, moléculas, y pH de docking.

### 2. Correr pipeline

Desde la raíz del proyecto (o como Run Configurations en PyCharm):

```bash
# 00a — Parsear moléculas
python 02_scripts/00a_molecule_parser.py --config 03_configs/00a_molecule_parser.yaml --campaign 04_data/campaigns/mi_campana/campaign_config.yaml

# 00b — Preparar receptor (ChimeraX + PDB2PQR)
python 02_scripts/00b_receptor_preparation.py --config 03_configs/00b_receptor_preparation.yaml --campaign 04_data/campaigns/mi_campana/campaign_config.yaml

# 00c — Protonar ligandos
python 02_scripts/00c_ionization_profiling.py --config 03_configs/00c_ionization_profiling.yaml --campaign 04_data/campaigns/mi_campana/campaign_config.yaml

# 00d — Antechamber (AM1-BCC charges)
python 02_scripts/00d_antechamber_preparation.py --config 03_configs/00d_antechamber_preparation.yaml --campaign 04_data/campaigns/mi_campana/campaign_config.yaml

# 01a — Generar grids
python 02_scripts/01a_grid_generation.py --config 03_configs/01a_grid_generation.yaml --campaign 04_data/campaigns/mi_campana/campaign_config.yaml

# 01b — DOCK6 docking
python 02_scripts/01b_dock6_run.py --config 03_configs/01b_dock6_run.yaml --campaign 04_data/campaigns/mi_campana/campaign_config.yaml
```

### PyCharm Run Configuration

Para cada módulo crear un Run Configuration:

| Campo | Valor |
|---|---|
| Script | `02_scripts/00b_receptor_preparation.py` |
| Parameters | `--config 03_configs/00b_receptor_preparation.yaml --campaign 04_data/campaigns/mi_campana/campaign_config.yaml` |
| Working directory | Raíz del proyecto |
| Python interpreter | `molecular_docking_env` |

## Configuración

Cada campaña se define en `campaign_config.yaml`. Los parámetros clave:

```yaml
campaign_id: "mi_campana"
docking_ph: 7.2

receptor:
  pdb: "receptor/mi_receptor.pdb"
  protonation:
    enabled: true
    tool: "pdb2pqr"       # pdb2pqr | chimerax | obabel

molecules:
  input_file: "molecules/"
  protonation:
    tool: "obabel"        # obabel | dimorphite_dl

grids:
  generate: true
  binding_site:
    method: "reference_ligand"
    reference_mol2: "molecules/ligando_cristalografico.mol2"
    radius: 10.0
```

## Notas técnicas

**DOCK6 y el límite de 80 caracteres.** Los programas Fortran de DOCK6 (sphgen, showbox) truncan paths a ~80 chars. El pipeline usa symlinks y filenames cortos automáticamente — no necesitas preocuparte de esto.

**Receptor preparation.** El módulo 00b usa ChimeraX para generar el mol2 del receptor con Sybyl atom types y AMBER ff14SB charges. Si usas la estrategia `pdb2pqr`, PDB2PQR+PROPKA predicen pKa por residuo antes de que ChimeraX asigne las cargas.

**Coordenadas cristalográficas.** Si tienes un ligando co-cristalizado, asegúrate de usar las coordenadas extraídas del PDB (no regeneradas desde SMILES) para `binding_site.reference_mol2`.

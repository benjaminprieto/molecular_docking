# MOLECULAR_DOCKING

Pipeline de virtual screening con dos engines complementarios: **DOCK6** (físico, basado en grids) y **GNINA** (Vina + CNN scoring). Produce outputs compatibles con **dock2profile**.

---

## Requisitos

| Dependencia | Instalación | Notas |
|---|---|---|
| Python 3.9-3.12 | conda | Vía `environment.yaml` |
| RDKit, OpenBabel, AmberTools | conda | Cheminformática + cargas AM1-BCC |
| scikit-learn | conda | Clustering y hotspots (m05) |
| PDB2PQR | pip | Protonación pH-aware del receptor |
| **DOCK6** | Manual | Licencia académica gratuita de [UCSF](https://dock.compbio.ucsf.edu/DOCK_6/index.htm) — engine 01x |
| **ChimeraX** | Manual | [Descargar](https://www.cgl.ucsf.edu/chimerax/download.html) — necesario para preparar el receptor (00b) |
| **GNINA** | Manual | [github.com/gnina/gnina](https://github.com/gnina/gnina) — engine 02x; CUDA recomendado para GPU |

---

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

Requiere licencia académica gratuita. Tras obtenerla:

```bash
# Compilar e instalar
tar -xzf dock.6.X.tar.gz
cd dock6
./configure gnu
make

# Configurar entorno (en ~/.bashrc o ~/.profile)
export DOCK_HOME=/opt/dock6           # raíz de la instalación
export PATH=$DOCK_HOME/bin:$PATH
```

El pipeline localiza los archivos de parámetros (`vdw_AMBER_parm99.defn`, `flex.defn`, `flex_drive.tbl`) en este orden:

1. `$DOCK_HOME / $DOCK6_HOME / $DOCK_BASE` (env vars)
2. `dirname(which dock6)/../parameters/`
3. Rutas comunes: `/opt/dock6`, `/usr/local/dock6`, `~/dock6`, `~/software/dock6`

### ChimeraX

Descargar e instalar desde [cgl.ucsf.edu/chimerax](https://www.cgl.ucsf.edu/chimerax/download.html). El pipeline busca el binario en `chimerax-daily` o `chimerax` en PATH.

### GNINA

Binario standalone (incluye Vina + CNN scoring entrenado). Soporta CPU y GPU (CUDA):

```bash
# Opción 1 — descargar binario precompilado
wget https://github.com/gnina/gnina/releases/latest/download/gnina
chmod +x gnina && sudo mv gnina /usr/local/bin/

# Opción 2 — compilar desde fuente (ver repo upstream)
```

GPU acelera 02b ~10-50x dependiendo del hardware. Si no hay GPU disponible, GNINA cae a CPU automáticamente.

---

## Estructura del proyecto

```
molecular_docking/
├── 01_src/molecular_docking/       Core modules (lógica, sin CLI)
│   ├── m00_preparation/            Preparación compartida (00a-00d)
│   ├── m01_docking/                DOCK6 engine (01a-01e)
│   ├── m02_gnina/                  GNINA engine (02a-02c)
│   ├── m04_dock6_analysis/         Análisis DOCK6 (04a-04e)
│   ├── m05_gnina_analysis/         Análisis GNINA (05a-05g)
│   └── m07_cross_engine/           Comparación cruzada (07a)
├── 02_scripts/                     CLI wrappers (argparse + YAML → core)
├── 03_configs/                     YAML por módulo (parámetros algorítmicos)
├── 04_data/campaigns/              Campañas (receptor + moléculas + grids)
│   └── <campaign_id>/
│       ├── campaign_config.yaml    Fuente de verdad de la campaña
│       ├── receptor/               PDB del receptor
│       ├── molecules/              Moléculas a dockear (SDF/CSV/SMILES)
│       ├── reference/              Ligando de referencia (cristalográfico)
│       └── grids/                  Pre-existentes (opcional; 01b los genera)
├── 05_results/                     Outputs por campaña/módulo (gitignored)
├── environment.yaml
├── pyproject.toml
├── check_dependencies.sh
└── run_pipeline.sh
```

### Convención de dos capas

Cada módulo tiene exactamente dos archivos:
1. **Core** (`01_src/.../module.py`): lógica pura, recibe paths/dicts, devuelve resultados.
2. **Script** (`02_scripts/module.py`): wrapper CLI que mergea YAML config + campaign_config + overrides y llama al core.

---

## Flujo del pipeline

```
Shared preparation (00x):
  00a  Molecule parser              → unique_molecules.csv + .sdf
  00b  Receptor preparation         → rec_charged.mol2 + rec_noH.pdb
  00c  Ionization profiling         → SDF protonados al pH de docking
  00d  Binding site definition      → rec_noH_site.pdb (opcional)

DOCK6 engine (01x):
  01a  Antechamber preparation      → mol2 con cargas AM1-BCC + Sybyl types
  01b  Grid generation              → DMS, spheres, box, grid.nrg/bmp
  01c  DOCK6 docking                → scored mol2 por molécula
  01d  Footprint re-scoring         → mol2 con descomposición por residuo
  01e  Score collection             → Excel compatible dock2profile

GNINA engine (02x):
  02a  GNINA preparation            → receptor.pdbqt + ligandos PDBQT
  02b  GNINA docking                → poses + CNN affinity
  02c  GNINA score collection       → CSV/Excel

DOCK6 analysis (04x):
  04a  Score ranking                → ranking compuesto multiplicativo
  04b  Footprint analysis           → consistencia residuo-residuo
  04c  Binding modes                → clusters de poses
  04d  Contact mapping              → contactos receptor-ligando
  04e  Campaign report              → Excel + plots resumen

GNINA analysis (05x):
  05a  Parse & fragment             → parseo + fragmentación BRICS
  05b  Fragment clustering          → DBSCAN sobre fragmentos
  05c  Score decomposition          → contribución por fragmento
  05d  Binding site hotspots        → hotspots por DBSCAN espacial
  05e  Structure export             → mejores poses como mol2/pdb
  05f  Contact mapping              → contactos (consume receptor PDB)
  05g  Campaign report              → Excel + plots resumen

Cross-engine (07x, requiere ambos engines):
  07a  Cross-engine comparison      → consensus DOCK6 ↔ GNINA
  07b  Hit export                   → exportación final de hits
```

---

## Uso

### 1. Crear campaña

```bash
mkdir -p 04_data/campaigns/mi_campana/{receptor,molecules,reference,grids,selections}
cp 04_data/campaigns/<campaña_existente>/campaign_config.yaml 04_data/campaigns/mi_campana/
```

Editar `campaign_config.yaml` con receptor, moléculas, pH y sitio de unión.

### 2. Correr pipeline completo

```bash
# Ambos engines, todos los módulos
bash run_pipeline.sh mi_campana

# Solo DOCK6
bash run_pipeline.sh mi_campana dock6

# Solo GNINA
bash run_pipeline.sh mi_campana gnina

# Rango específico (ej: resumir desde 01c hasta 01e)
bash run_pipeline.sh mi_campana dock6 01c 01e
```

### 3. Correr módulos individuales

Todos los scripts siguen el mismo patrón:

```bash
python 02_scripts/<modulo>.py \
    --config 03_configs/<modulo>.yaml \
    --campaign 04_data/campaigns/mi_campana/campaign_config.yaml
```

### 4. Background runs (servidor)

```bash
mkdir -p logs

# Pipeline completo
nohup bash run_pipeline.sh mi_campana both > logs/mi_campana.log 2>&1 &

# Monitoreo
jobs -l
tail -f logs/mi_campana.log
```

### PyCharm Run Configurations

| Campo | Valor |
|---|---|
| Script | `02_scripts/00b_receptor_preparation.py` |
| Parameters | `--config 03_configs/00b_receptor_preparation.yaml --campaign 04_data/campaigns/mi_campana/campaign_config.yaml` |
| Working directory | Raíz del proyecto |
| Python interpreter | `molecular_docking_env` |

---

## Configuración de campaña

`campaign_config.yaml` es la fuente de verdad. Ejemplo mínimo:

```yaml
campaign_id: "mi_campana"
description: "Cribado virtual de N compuestos contra <target>"

receptor:
  pdb: "receptor/mi_receptor.pdb"
  protonation:
    enabled: true
    tool: "pdb2pqr"          # pdb2pqr | chimerax | obabel

docking_ph: 7.2

molecules:
  input_file: "molecules/"
  protonation:
    enabled: true
    tool: "obabel"           # obabel | dimorphite_dl

grids:
  generate: true
  binding_site:
    method: "reference_ligand"
    reference_mol2: "reference/ligando_cristal.mol2"
    radius: 6.0
```

Los YAML en `03_configs/` controlan solo parámetros algorítmicos (exhaustiveness, charge method, conformer strategy, etc.).

---

## Notas técnicas

**Límite de 80 caracteres en DOCK6.** Los programas Fortran de DOCK6 (sphgen, showbox) truncan paths a ~80 chars. El pipeline usa symlinks y filenames cortos automáticamente.

**SYBYL atom types obligatorios.** DOCK6 flexible docking requiere SYBYL (no GAFF), si no cae silenciosamente a rigid docking. `01a_antechamber_preparation` lo enforce.

**Cargas AMBER ff14SB en receptor.** El módulo `00b_receptor_preparation` genera el mol2 del receptor con Sybyl atom types y AMBER ff14SB charges vía ChimeraX. Si la estrategia es `pdb2pqr`, PDB2PQR+PROPKA predicen pKa por residuo antes de la asignación de cargas.

**Coordenadas cristalográficas.** Para `binding_site.reference_mol2`, usar coordenadas extraídas del PDB (no regeneradas desde SMILES) para definir el sitio de unión correctamente.

---

## Testing

```bash
pytest                                              # smoke tests completos
pytest tests/test_pipeline.py::TestGridGeneration   # una clase
pytest -k "test_parse_vina_output"                  # un test
```

Los tests no requieren DOCK6/ChimeraX/GNINA instalados; verifican imports y lógica básica.
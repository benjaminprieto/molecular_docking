#!/bin/bash
# =============================================================================
# molecular_docking — Dependency Checker
# =============================================================================
# Verifica que las dependencias del pipeline esten instaladas y accesibles.
#
# Usage:
#   bash check_dependencies.sh
#
# Salida con exit code != 0 si falta algo critico (FAIL). Los WARN no
# bloquean: indican modulos que no podras correr hasta instalar la dep.
# =============================================================================

set -u

echo "============================================================"
echo "  MOLECULAR_DOCKING — Dependency Check"
echo "============================================================"
echo ""

PASS=0; FAIL=0; WARN=0

ok()   { echo "  [OK]   $1"; PASS=$((PASS+1)); }
fail() { echo "  [FAIL] $1"; FAIL=$((FAIL+1)); }
warn() { echo "  [WARN] $1"; WARN=$((WARN+1)); }

check_cmd() {
    local name="$1" cmd="$2" req="$3"
    if eval "$cmd" > /dev/null 2>&1; then
        ok "$name"
    elif [ "$req" = "required" ]; then
        fail "$name: NOT FOUND in PATH"
    else
        warn "$name: not found (optional — module unavailable)"
    fi
}

check_pylib() {
    local lib="$1" req="$2"
    if python3 -c "import $lib" 2>/dev/null; then
        ok "python: $lib"
    elif [ "$req" = "required" ]; then
        fail "python: $lib"
    else
        warn "python: $lib (optional)"
    fi
}

# --- Python ---
echo "--- Python ---"
check_cmd "python3" "python3 --version" "required"
echo ""

# --- Python libraries ---
echo "--- Python Libraries ---"
for lib in rdkit pandas numpy openpyxl yaml sklearn; do
    check_pylib "$lib" "required"
done
echo ""

# --- DOCK6 (engine 01x + analysis 04x) ---
echo "--- DOCK6 (engine 01x) ---"
check_cmd "dock6"           "which dock6"           "required"
check_cmd "grid"            "which grid"            "required"
check_cmd "sphgen"          "which sphgen"          "required"
check_cmd "sphere_selector" "which sphere_selector" "required"
check_cmd "showbox"         "which showbox"         "required"
check_cmd "dms"             "which dms"             "required"

# DOCK6 parameter files (needed by grid_generation + footprint_rescore)
PARAM_FILE="parameters/vdw_AMBER_parm99.defn"
if [ -n "${DOCK_HOME:-}" ] && [ -f "$DOCK_HOME/$PARAM_FILE" ]; then
    ok "DOCK6 parameters (\$DOCK_HOME/$PARAM_FILE)"
elif [ -f "/opt/dock6/$PARAM_FILE" ]; then
    ok "DOCK6 parameters (/opt/dock6/$PARAM_FILE)"
else
    warn "DOCK6 parameters not found — set DOCK_HOME=/path/to/dock6"
fi
echo ""

# --- ChimeraX (00b receptor preparation) ---
echo "--- ChimeraX (00b receptor preparation) ---"
if which chimerax > /dev/null 2>&1; then
    ok "chimerax ($(which chimerax))"
elif which chimerax-daily > /dev/null 2>&1; then
    ok "chimerax-daily ($(which chimerax-daily))"
elif [ -x "/usr/bin/chimerax-daily" ]; then
    ok "chimerax-daily (/usr/bin/chimerax-daily)"
else
    warn "ChimeraX: not found — install from https://www.cgl.ucsf.edu/chimerax/download.html"
fi
echo ""

# --- AmberTools (01a antechamber preparation) ---
echo "--- AmberTools (01a antechamber preparation) ---"
check_cmd "antechamber" "which antechamber" "required"
check_cmd "parmchk2"    "which parmchk2"    "required"
check_cmd "reduce"      "which reduce"      "optional"
echo ""

# --- OpenBabel (CLI used by 00b/00c/01a/06) ---
echo "--- OpenBabel ---"
check_cmd "obabel" "obabel -V" "required"
echo ""

# --- pdb2pqr (default protonator in 00b) ---
echo "--- pdb2pqr (00b receptor protonation) ---"
check_cmd "pdb2pqr" "which pdb2pqr" "required"
echo ""

# --- GNINA (engine 02x) ---
echo "--- GNINA (engine 02x) ---"
if which gnina > /dev/null 2>&1; then
    GNINA_VER=$(gnina --version 2>&1 | head -n1)
    ok "gnina ($(which gnina)) — ${GNINA_VER}"
elif [ -x "/usr/local/bin/gnina" ]; then
    ok "gnina (/usr/local/bin/gnina)"
elif [ -x "$HOME/bin/gnina" ]; then
    ok "gnina ($HOME/bin/gnina)"
else
    warn "gnina: not found — engine 02x will fail. Install from https://github.com/gnina/gnina"
fi
echo ""

# --- Summary ---
echo "============================================================"
echo "  SUMMARY: $PASS passed, $FAIL failed, $WARN warnings"
echo "============================================================"
if [ $FAIL -gt 0 ]; then
    echo "  Fix FAIL items before running the pipeline."
    echo ""
    echo "  Quick fix (Python deps):"
    echo "    conda activate molecular_docking_env"
    echo "    conda env update -f environment.yaml --prune"
    echo "    pip install -e '.[dev]'"
    echo ""
    echo "  External tools (DOCK6, ChimeraX, GNINA): see environment.yaml header."
    exit 1
fi
if [ $WARN -gt 0 ]; then
    echo "  Some optional modules will be unavailable (see WARN above)."
fi
echo "  Ready."
exit 0

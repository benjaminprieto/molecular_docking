#!/bin/bash
# =============================================================================
# molecular_docking - Dependency Checker
# =============================================================================

echo "============================================================"
echo "  MOLECULAR_DOCKING - Dependency Check"
echo "============================================================"
echo ""

PASS=0; FAIL=0; WARN=0

check_cmd() {
    local name="$1" cmd="$2" req="$3"
    if eval "$cmd" > /dev/null 2>&1; then
        echo "  [OK]   $name"
        ((PASS++))
    elif [ "$req" = "required" ]; then
        echo "  [FAIL] $name: NOT FOUND"
        ((FAIL++))
    else
        echo "  [WARN] $name: not found (optional)"
        ((WARN++))
    fi
}

echo "--- Python ---"
check_cmd "python3"     "python3 --version"         "required"

echo ""
echo "--- Python Libraries ---"
for lib in rdkit pandas numpy openpyxl yaml; do
    mod=$lib; [ "$lib" = "yaml" ] && mod="yaml"
    if python3 -c "import $mod" 2>/dev/null; then
        echo "  [OK]   $lib"; ((PASS++))
    else
        echo "  [FAIL] $lib"; ((FAIL++))
    fi
done
for lib in dimorphite_dl pdb2pqr openbabel; do
    if python3 -c "import $lib" 2>/dev/null; then
        echo "  [OK]   $lib (optional)"; ((PASS++))
    else
        echo "  [WARN] $lib (optional)"; ((WARN++))
    fi
done

echo ""
echo "--- DOCK6 ---"
check_cmd "dock6"            "which dock6"            "required"
check_cmd "grid"             "which grid"             "optional"
check_cmd "sphgen"           "which sphgen"           "optional"
check_cmd "sphere_selector"  "which sphere_selector"  "optional"
check_cmd "showbox"          "which showbox"          "optional"

echo ""
echo "--- AmberTools ---"
check_cmd "antechamber"  "which antechamber"  "required"
check_cmd "parmchk2"     "which parmchk2"     "optional"
check_cmd "tleap"        "which tleap"        "optional"
check_cmd "reduce"       "which reduce"       "optional"

echo ""
echo "--- OpenBabel ---"
check_cmd "obabel"  "obabel -V"  "required"

echo ""
echo "--- Optional ---"
check_cmd "chimera"  "which chimera"  "optional"

echo ""
echo "============================================================"
echo "  SUMMARY: $PASS passed, $FAIL failed, $WARN warnings"
echo "============================================================"
[ $FAIL -gt 0 ] && echo "  Fix FAIL items before running pipeline" && exit 1
echo "  Ready!" && exit 0

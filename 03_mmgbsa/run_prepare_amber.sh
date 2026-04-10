#!/bin/bash
#
# run_prepare_amber.sh
#
# Stage 1 of the MM-GBSA pipeline: generate Amber topology/coordinate
# files from a Boltz-2 CIF prediction.
#
# Usage:
#   sbatch run_prepare_amber.sh <boltz_cif> <design_dir>
#
# Example:
#   sbatch run_prepare_amber.sh \
#       $SCRATCH/boltz/outputs/.../design_01_seed123_model_0.cif \
#       $SCRATCH/mmgbsa/design_01
#
# Idempotent: if the expected outputs already exist, this script exits
# successfully without re-running tLeap. To force a redo, delete the
# output files first.

#SBATCH -A scholar
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --job-name prep_amber
#SBATCH --output=%x-%J-%u.out
#SBATCH --error=%x-%J-%u.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bsdrown@purdue.edu

set -euo pipefail

if [[ -z "${1:-}" || -z "${2:-}" ]]; then
    echo "Error: must provide boltz CIF and design dir"
    echo "Usage: sbatch run_prepare_amber.sh <boltz_cif> <design_dir>"
    exit 1
fi

BOLTZ_CIF="$1"
DESIGN_DIR="$2"

if [[ ! -f "$BOLTZ_CIF" ]]; then
    echo "Error: CIF file '$BOLTZ_CIF' not found"
    exit 1
fi

# Idempotence check. If all expected outputs already exist, we're done.
EXPECTED=(
    "$DESIGN_DIR/complex.prmtop"
    "$DESIGN_DIR/complex.inpcrd"
    "$DESIGN_DIR/receptor.prmtop"
    "$DESIGN_DIR/receptor.inpcrd"
    "$DESIGN_DIR/ligand.prmtop"
    "$DESIGN_DIR/ligand.inpcrd"
)
ALL_EXIST=true
for f in "${EXPECTED[@]}"; do
    if [[ ! -s "$f" ]]; then
        ALL_EXIST=false
        break
    fi
done
if [[ "$ALL_EXIST" == "true" ]]; then
    echo "=== Amber prep SKIPPED ==="
    echo "All expected outputs already exist in $DESIGN_DIR"
    echo "To force redo, remove: ${EXPECTED[*]}"
    exit 0
fi

module purge 2>/dev/null || true
ml ambertools/25

# Locate prepare_amber.py
PREP_PY="${SLURM_SUBMIT_DIR:-$(pwd)}/prepare_amber.py"
if [[ ! -f "$PREP_PY" ]]; then
    PREP_PY="./prepare_amber.py"
fi
if [[ ! -f "$PREP_PY" ]]; then
    echo "Error: cannot find prepare_amber.py"
    exit 1
fi

echo "=== Stage 1: Amber Parameter Generation ==="
echo "Date:       $(date)"
echo "Node:       $(hostname)"
echo "CIF:        $BOLTZ_CIF"
echo "Design dir: $DESIGN_DIR"
echo "Script:     $PREP_PY"
echo ""

python "$PREP_PY" --boltz-cif "$BOLTZ_CIF" --output-dir "$DESIGN_DIR"

echo ""
echo "=== Amber prep complete ==="

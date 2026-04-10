#!/bin/bash
#
# run_mmgbsa.sh
#
# SLURM script to run MM-GBSA analysis on a complex trajectory.
# Uses Scholar's ambertools/25 module for MMPBSA.py and the openmm
# conda env for parmed.
#
# Usage:
#   sbatch run_mmgbsa.sh <design_dir>
#
# Example:
#   sbatch run_mmgbsa.sh $SCRATCH/mmgbsa/design_01

#SBATCH -A scholar
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --job-name mmgbsa
#SBATCH --output=%x-%J-%u.out
#SBATCH --error=%x-%J-%u.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bsdrown@purdue.edu

set -euo pipefail

if [[ -z "${1:-}" ]]; then
    echo "Error: must provide a design directory"
    echo "Usage: sbatch run_mmgbsa.sh <design_dir>"
    exit 1
fi

DESIGN_DIR="$1"
if [[ ! -d "$DESIGN_DIR" ]]; then
    echo "Error: directory '$DESIGN_DIR' not found"
    exit 1
fi

# Idempotence: skip if MMPBSA results already exist
if [[ -s "$DESIGN_DIR/mmgbsa.dat" && -s "$DESIGN_DIR/mmgbsa.csv" ]]; then
    # Verify the dat file has a DELTA TOTAL line (complete run)
    if grep -q "DELTA TOTAL" "$DESIGN_DIR/mmgbsa.dat"; then
        echo "=== MM-GBSA SKIPPED ==="
        echo "mmgbsa.dat and mmgbsa.csv already exist in $DESIGN_DIR"
        echo "To force redo, remove those files first."
        exit 0
    else
        echo "mmgbsa.dat exists but looks incomplete (no DELTA TOTAL line)"
        echo "Re-running MMPBSA..."
    fi
fi

# Verify upstream inputs
if [[ ! -s "$DESIGN_DIR/production.nc" ]]; then
    echo "Error: production.nc missing in $DESIGN_DIR"
    echo "Run run_md.sh first."
    exit 1
fi

module purge 2>/dev/null || true
ml conda
source activate "$SCRATCH/conda/envs/openmm"

# Do NOT load ambertools/25 as a module here. Loading it puts its
# Python 3.12 and bundled (outdated) parmed on the front of the path,
# which shadows our conda env's modern parmed and breaks numpy
# compatibility. Instead we call MMPBSA.py via its absolute path.
export MMPBSA_PY="/apps/external/ambertools/25/ambertools25/bin/MMPBSA.py"
if [[ ! -x "$MMPBSA_PY" ]]; then
    echo "Error: MMPBSA.py not found at $MMPBSA_PY"
    exit 1
fi

RUN_PY="${SLURM_SUBMIT_DIR}/run_mmgbsa.py"
if [[ ! -f "$RUN_PY" ]]; then
    RUN_PY="./run_mmgbsa.py"
fi
if [[ ! -f "$RUN_PY" ]]; then
    echo "Error: cannot find run_mmgbsa.py"
    exit 1
fi

echo "=== MM-GBSA ==="
echo "Date:       $(date)"
echo "Node:       $(hostname)"
echo "Design dir: $DESIGN_DIR"
echo "Script:     $RUN_PY"
echo "MMPBSA.py:  $MMPBSA_PY"
echo ""

python "$RUN_PY" --input-dir "$DESIGN_DIR"

echo ""
echo "=== MM-GBSA complete ==="

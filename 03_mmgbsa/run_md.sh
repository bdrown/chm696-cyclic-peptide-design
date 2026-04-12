#!/bin/bash
#
# run_md.sh
#
# SLURM script for running OpenMM MD (heating + production) on a
# minimized MDM2-cyclic peptide complex.
#
# Usage:
#   sbatch run_md.sh <design_dir>
#
# Example:
#   sbatch run_md.sh $SCRATCH/mmgbsa/design_01
#
# The input directory must contain complex.prmtop and complex_min.rst7
# produced by prepare_amber.py and minimize.py.
# Outputs: production.nc, production.rst7, md.log

#SBATCH -A gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --gpus-per-node=1
#SBATCH --constraint=J
#SBATCH --time=00:30:00
#SBATCH --job-name md
#SBATCH --output=%x-%J-%u.out
#SBATCH --error=%x-%J-%u.err

set -euo pipefail

if [[ -z "${1:-}" ]]; then
    echo "Error: must provide a design directory"
    echo "Usage: sbatch run_md.sh <design_dir>"
    exit 1
fi

DESIGN_DIR="$1"
if [[ ! -d "$DESIGN_DIR" ]]; then
    echo "Error: directory '$DESIGN_DIR' not found"
    exit 1
fi

# Idempotence: skip if trajectory already exists. We also check that it
# has the expected ~500 frames (file size > 10 MB is a good proxy for
# the ~15 MB netcdf from a 10 ns run of a ~2000 atom system).
if [[ -s "$DESIGN_DIR/production.nc" && -s "$DESIGN_DIR/production.rst7" ]]; then
    NC_SIZE=$(stat -c %s "$DESIGN_DIR/production.nc")
    if [[ "$NC_SIZE" -gt 1000000 ]]; then
        echo "=== MD SKIPPED ==="
        echo "production.nc (${NC_SIZE} bytes) and production.rst7 "
        echo "already exist in $DESIGN_DIR"
        echo "To force redo, remove those files first."
        exit 0
    else
        echo "production.nc exists but is suspiciously small (${NC_SIZE} bytes)"
        echo "Removing and re-running MD..."
        rm -f "$DESIGN_DIR/production.nc" "$DESIGN_DIR/production.rst7"
    fi
fi

# Verify upstream inputs
if [[ ! -s "$DESIGN_DIR/complex_min.rst7" ]]; then
    echo "Error: complex_min.rst7 missing in $DESIGN_DIR"
    echo "Run minimize.py (or run_minimize.sh) first."
    exit 1
fi

module purge 2>/dev/null || true
ml conda
source activate "$SCRATCH/conda/envs/openmm"

# Locate run_md.py
RUN_MD_PY="${SLURM_SUBMIT_DIR}/run_md.py"
if [[ ! -f "$RUN_MD_PY" ]]; then
    RUN_MD_PY="./run_md.py"
fi
if [[ ! -f "$RUN_MD_PY" ]]; then
    echo "Error: cannot find run_md.py in $SLURM_SUBMIT_DIR or $(pwd)"
    exit 1
fi

echo "=== OpenMM MD ==="
echo "Date:       $(date)"
echo "Node:       $(hostname)"
echo "GPU:        $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "Design dir: $DESIGN_DIR"
echo "Script:     $RUN_MD_PY"
echo ""

python "$RUN_MD_PY" --input-dir "$DESIGN_DIR"

echo ""
echo "=== MD complete ==="
ls -la "$DESIGN_DIR"/production.* 2>/dev/null || echo "(no trajectory found)"

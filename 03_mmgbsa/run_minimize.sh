#!/bin/bash
#
# run_minimize.sh
#
# SLURM script for running OpenMM energy minimization on an
# Amber-prepared MDM2-cyclic peptide complex.
#
# Usage:
#   sbatch run_minimize.sh <design_dir>
#
# Example:
#   sbatch run_minimize.sh $SCRATCH/mmgbsa/design_01
#
# The input directory must contain complex.prmtop and complex.inpcrd
# produced by prepare_amber.py. Outputs (complex_min.pdb,
# complex_min.rst7, minimize.log) are written to the same directory.

#SBATCH -A gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --gpus-per-node=1
#SBATCH --constraint=J
#SBATCH --time=00:15:00
#SBATCH --job-name minimize
#SBATCH --output=%x-%J-%u.out
#SBATCH --error=%x-%J-%u.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${USER}@purdue.edu

set -euo pipefail

if [[ -z "${1:-}" ]]; then
    echo "Error: must provide a design directory containing complex.prmtop"
    echo "Usage: sbatch run_minimize.sh <design_dir>"
    exit 1
fi

DESIGN_DIR="$1"
if [[ ! -d "$DESIGN_DIR" ]]; then
    echo "Error: directory '$DESIGN_DIR' not found"
    exit 1
fi

# Idempotence: skip if already done
if [[ -s "$DESIGN_DIR/complex_min.rst7" && -s "$DESIGN_DIR/complex_min.pdb" ]]; then
    echo "=== Minimization SKIPPED ==="
    echo "complex_min.rst7 and complex_min.pdb already exist in $DESIGN_DIR"
    echo "To force redo, remove those files first."
    exit 0
fi

# Verify upstream inputs exist
if [[ ! -s "$DESIGN_DIR/complex.prmtop" || ! -s "$DESIGN_DIR/complex.inpcrd" ]]; then
    echo "Error: complex.prmtop or complex.inpcrd missing in $DESIGN_DIR"
    echo "Run prepare_amber.py (or run_prepare_amber.sh) first."
    exit 1
fi

module purge 2>/dev/null || true
ml conda
source activate "$SCRATCH/conda/envs/openmm"

# Locate minimize.py. First try the directory the sbatch was launched
# from (the normal case when students submit from the evaluation dir),
# then fall back to the current working directory.
MINIMIZE_PY="${SLURM_SUBMIT_DIR}/minimize.py"
if [[ ! -f "$MINIMIZE_PY" ]]; then
    MINIMIZE_PY="./minimize.py"
fi
if [[ ! -f "$MINIMIZE_PY" ]]; then
    echo "Error: cannot find minimize.py in $SLURM_SUBMIT_DIR or $(pwd)"
    exit 1
fi

echo "=== OpenMM Minimization ==="
echo "Date:       $(date)"
echo "Node:       $(hostname)"
echo "GPU:        $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "Design dir: $DESIGN_DIR"
echo "Script:     $MINIMIZE_PY"
echo ""

python "$MINIMIZE_PY" --input-dir "$DESIGN_DIR"

echo ""
echo "=== Minimization complete ==="
ls -la "$DESIGN_DIR"/complex_min.* 2>/dev/null || echo "(no minimized outputs found)"

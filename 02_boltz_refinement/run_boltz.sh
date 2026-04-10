#!/bin/bash
#
# run_boltz.sh
#
# SLURM submission script for running Boltz-2 on Scholar GPU nodes.
# Targets the A40 48GB GPUs (sub-cluster J) which have sufficient
# VRAM for ternary complex prediction.
#
# Usage:
#   sbatch run_boltz.sh <input.yaml>
#
# Example:
#   sbatch run_boltz.sh mz1_ternary.yaml
#
# If no argument is provided, Boltz will process all YAML files in
# the inputs/ directory.

#SBATCH -A gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gpus-per-node=1
#SBATCH --constraint=J
#SBATCH --time=00:45:00
#SBATCH --job-name boltz_predict
#SBATCH --output=%x-%J-%u.out
#SBATCH --error=%x-%J-%u.err

# -----------------------------------------------------------------------
# Environment setup
# -----------------------------------------------------------------------
ml conda
source activate "$SCRATCH/conda/envs/boltz2"

# Ensure caches use scratch, not home (redundant with conda activate
# hooks but guards against edge cases in SLURM environments)
export PIP_CACHE_DIR="$SCRATCH/conda/pip-cache"
export XDG_CACHE_HOME="$SCRATCH/conda/xdg-cache"
export BOLTZ_CACHE="$SCRATCH/conda/boltz-cache"

# -----------------------------------------------------------------------
# Directories
# -----------------------------------------------------------------------
BOLTZ_WORKDIR="$SCRATCH/boltz"
INPUT_DIR="$BOLTZ_WORKDIR/inputs"
OUTPUT_DIR="$BOLTZ_WORKDIR/outputs"

mkdir -p "$OUTPUT_DIR"

# -----------------------------------------------------------------------
# Determine input
# -----------------------------------------------------------------------
if [[ -n "${1:-}" ]]; then
    # Specific input file provided as argument
    INPUT_PATH="$INPUT_DIR/$1"
    if [[ ! -f "$INPUT_PATH" ]]; then
        echo "Error: input file '$INPUT_PATH' not found."
        exit 1
    fi
else
    # Process all YAML files in the inputs directory
    INPUT_PATH="$INPUT_DIR"
fi

# -----------------------------------------------------------------------
# Report environment
# -----------------------------------------------------------------------
echo "=== Boltz-2 Prediction ==="
echo "Date:      $(date)"
echo "Node:      $(hostname)"
echo "GPU:       $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "Input:     $INPUT_PATH"
echo "Output:    $OUTPUT_DIR"
echo ""

# -----------------------------------------------------------------------
# Run Boltz prediction
# -----------------------------------------------------------------------
# --use_msa_server: generates MSAs using the Boltz MSA server
#                   (requires internet access from compute node)
# --diffusion_samples 5: generate 5 structural models to pick from
#
# If the compute nodes do not have internet access, you will need to
# pre-compute MSAs on a login node or provide them in the YAML file.
# See: https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md

boltz predict "$INPUT_PATH" \
    --out_dir "$OUTPUT_DIR" \
    --use_msa_server \
    --diffusion_samples 5 \
    --num_workers 4

echo ""
echo "=== Prediction complete ==="
echo "Results in: $OUTPUT_DIR"
ls -la "$OUTPUT_DIR"

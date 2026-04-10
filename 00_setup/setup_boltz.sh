#!/bin/bash
#
# setup_boltz.sh
#
# One-time setup of a conda environment for Boltz-2 on Scholar.
# Run this script from a login node (NOT in a SLURM job).
#
# Usage:
#   bash setup_boltz.sh
#
# After setup, activate the environment with:
#   ml conda
#   source activate $SCRATCH/conda/envs/boltz2
#
# The environment, package cache, and pip cache are all placed in
# scratch to avoid exhausting your home directory quota.
# See: https://www.rcac.purdue.edu/knowledge/scholar/run/examples/apps/python/conda

set -euo pipefail

echo "=== Boltz-2 Environment Setup for Scholar ==="

# Load the conda module (required on RCAC systems)
ml conda

# -----------------------------------------------------------------------
# 0. Clean up any stale caches in home directory from failed installs
# -----------------------------------------------------------------------
if [[ -d "$HOME/.cache/pip" ]]; then
    echo "Removing old pip cache from home directory to free quota..."
    rm -rf "$HOME/.cache/pip"
fi

# -----------------------------------------------------------------------
# 1. Define paths — everything goes to scratch
# -----------------------------------------------------------------------
SCRATCH_DIR="${SCRATCH:-/scratch/scholar/$USER}"
ENV_PREFIX="$SCRATCH_DIR/conda/envs/boltz2"
CONDA_PKGS="$SCRATCH_DIR/conda/pkgs"
PIP_CACHE="$SCRATCH_DIR/conda/pip-cache"
XDG_CACHE="$SCRATCH_DIR/conda/xdg-cache"

mkdir -p "$CONDA_PKGS" "$PIP_CACHE" "$XDG_CACHE"

# Redirect pip and general caches BEFORE any pip commands run
export PIP_CACHE_DIR="$PIP_CACHE"
export XDG_CACHE_HOME="$XDG_CACHE"

# Tell conda to store downloaded packages in scratch
conda config --add pkgs_dirs "$CONDA_PKGS"

echo "Environment prefix: $ENV_PREFIX"
echo "Conda package cache: $CONDA_PKGS"
echo "Pip cache: $PIP_CACHE"

# -----------------------------------------------------------------------
# 2. Create a fresh Python 3.11 environment using --prefix
# -----------------------------------------------------------------------
# On RCAC systems, --name puts envs in $HOME/.conda regardless of config.
# Using --prefix places the env exactly where we specify.

if [[ -d "$ENV_PREFIX" ]]; then
    echo ""
    echo "Environment already exists at: $ENV_PREFIX"
    echo "To recreate, run:"
    echo "  rm -rf $ENV_PREFIX"
    echo "  bash setup_boltz.sh"
    exit 1
fi

echo ""
echo "Creating conda environment at $ENV_PREFIX with Python 3.11..."
conda create --prefix="$ENV_PREFIX" python=3.11 -y

# -----------------------------------------------------------------------
# 3. Activate and install Boltz
# -----------------------------------------------------------------------
echo ""
echo "Activating environment and installing Boltz..."

# RCAC recommends 'source activate' with the full prefix path
source activate "$ENV_PREFIX"

echo "Pip cache redirected to: $PIP_CACHE_DIR"

# Install PyTorch built for CUDA 12.6.
# Scholar's A40 nodes run driver 575.57.08 (CUDA 12.9), which supports
# cu126. We pin torch to cu126 (not the default pip torch) because
# Boltz-2 also needs cuEquivariance, whose CUDA ops require
# nvidia-cublas-cu12 >= 12.5. Torch cu124 ships cublas 12.4 (too old),
# while cu126 ships cublas 12.6+ (compatible with both).
echo "Installing PyTorch (CUDA 12.6 build)..."
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu126

# Install Boltz and cuEquivariance CUDA kernels.
# The NVIDIA cuEquivariance packages are hosted on pypi.nvidia.com.
echo "Installing Boltz + cuEquivariance kernels..."
pip install "boltz[cuda]" --extra-index-url https://pypi.nvidia.com

# -----------------------------------------------------------------------
# 3b. Make cache redirects persistent for this environment
# -----------------------------------------------------------------------
ACTIVATE_DIR="$ENV_PREFIX/etc/conda/activate.d"
DEACTIVATE_DIR="$ENV_PREFIX/etc/conda/deactivate.d"
mkdir -p "$ACTIVATE_DIR" "$DEACTIVATE_DIR"

cat > "$ACTIVATE_DIR/cache_dirs.sh" << 'ACTIVATE_EOF'
# Redirect caches to scratch to preserve home directory quota
export PIP_CACHE_DIR="${SCRATCH}/conda/pip-cache"
export XDG_CACHE_HOME="${SCRATCH}/conda/xdg-cache"
export BOLTZ_CACHE="${SCRATCH}/conda/boltz-cache"
mkdir -p "$PIP_CACHE_DIR" "$XDG_CACHE_HOME" "$BOLTZ_CACHE"
ACTIVATE_EOF

cat > "$DEACTIVATE_DIR/cache_dirs.sh" << 'DEACTIVATE_EOF'
unset PIP_CACHE_DIR
unset XDG_CACHE_HOME
unset BOLTZ_CACHE
DEACTIVATE_EOF

echo "Persistent cache redirects installed in conda environment."

# -----------------------------------------------------------------------
# 4. Verify installation
# -----------------------------------------------------------------------
echo ""
echo "Verifying installation..."

# Check that torch still has the CUDA version we pinned
TORCH_CUDA=$(python -c "import torch; print(torch.version.cuda)" 2>/dev/null)
echo "  PyTorch CUDA version: ${TORCH_CUDA:-FAILED}"
if [[ "$TORCH_CUDA" != 12.6* ]]; then
    echo "  WARNING: Expected CUDA 12.6 but got ${TORCH_CUDA}."
    echo "  The boltz install may have pulled in a different torch."
    echo "  Run: pip install torch --index-url https://download.pytorch.org/whl/cu126"
fi

python -c "import boltz; print(f'  Boltz version: {boltz.__version__}')" 2>/dev/null || \
    echo "  Boltz version: could not determine"

boltz predict --help > /dev/null 2>&1 && \
    echo "  boltz predict command: OK" || \
    echo "  WARNING: boltz predict command not found in PATH"

# -----------------------------------------------------------------------
# 5. Create convenience directories
# -----------------------------------------------------------------------
BOLTZ_WORKDIR="$SCRATCH_DIR/boltz"
mkdir -p "$BOLTZ_WORKDIR/inputs" "$BOLTZ_WORKDIR/outputs"

echo ""
echo "=== Setup Complete ==="
echo ""
echo "Working directory: $BOLTZ_WORKDIR"
echo "  inputs/:   place YAML input files here"
echo "  outputs/:  Boltz will write results here"
echo ""
echo "To use Boltz in a terminal session:"
echo "  ml conda"
echo "  source activate $ENV_PREFIX"
echo ""
echo "To use Boltz in a SLURM job, these lines are already in run_boltz.sh."
echo ""
echo "Note: Boltz model weights (~5 GB) will be downloaded automatically"
echo "on the first run. This download happens during the SLURM job, so"
echo "ensure your first job has sufficient walltime."

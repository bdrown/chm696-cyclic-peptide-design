#!/bin/bash
#
# setup_openmm.sh
#
# One-time setup of a conda environment for OpenMM-based MD simulation
# of cyclic peptide-MDM2 complexes. Run once from a Scholar login node.
#
# Usage:
#   bash setup_openmm.sh
#
# After setup, the environment is activated with:
#   ml conda
#   source activate $SCRATCH/conda/envs/openmm
#
# This environment is separate from boltz2 because OpenMM's CUDA stack
# has different requirements than PyTorch's, and mixing them is fragile.
# We use Scholar's ambertools/25 module directly for tleap/MMPBSA.py
# rather than conda-installing AmberTools (which would be slower and
# less well-tested than the RCAC-maintained module).

set -euo pipefail

echo "=== OpenMM Environment Setup for Scholar ==="

ml conda

# -----------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------
SCRATCH_DIR="${SCRATCH:-/scratch/scholar/$USER}"
ENV_PREFIX="$SCRATCH_DIR/conda/envs/openmm"
CONDA_PKGS="$SCRATCH_DIR/conda/pkgs"
PIP_CACHE="$SCRATCH_DIR/conda/pip-cache"

mkdir -p "$CONDA_PKGS" "$PIP_CACHE"

export PIP_CACHE_DIR="$PIP_CACHE"
conda config --add pkgs_dirs "$CONDA_PKGS" 2>/dev/null || true

echo "Environment prefix: $ENV_PREFIX"

if [[ -d "$ENV_PREFIX" ]]; then
    echo ""
    echo "Environment already exists. To recreate:"
    echo "  rm -rf $ENV_PREFIX"
    echo "  bash setup_openmm.sh"
    exit 1
fi

# -----------------------------------------------------------------------
# Create environment with OpenMM (CUDA 12 build) and supporting libs
# -----------------------------------------------------------------------
# parmed: read/write Amber prmtop files, apply modifications
# mdtraj: trajectory analysis, frame extraction for MMPBSA input
# numpy/scipy: general numerics
echo ""
echo "Creating conda environment with OpenMM..."
conda create --prefix="$ENV_PREFIX" -y \
    -c conda-forge \
    python=3.11 \
    "openmm>=8.1" \
    cuda-version=12 \
    parmed \
    mdtraj \
    numpy \
    scipy \
    matplotlib

source activate "$ENV_PREFIX"

# -----------------------------------------------------------------------
# Verify installation
# -----------------------------------------------------------------------
echo ""
echo "Verifying installation..."
python -c "
import openmm
import parmed
import mdtraj
print(f'OpenMM version: {openmm.version.version}')
print(f'parmed version: {parmed.__version__}')
print(f'mdtraj version: {mdtraj.__version__}')

# Check CUDA platform is available
from openmm import Platform
platforms = [Platform.getPlatform(i).getName()
             for i in range(Platform.getNumPlatforms())]
print(f'Available OpenMM platforms: {platforms}')
if 'CUDA' not in platforms:
    print('WARNING: CUDA platform not found. GPU acceleration will not work.')
"

echo ""
echo "=== Setup Complete ==="
echo ""
echo "To activate the environment:"
echo "  ml conda"
echo "  source activate $ENV_PREFIX"
echo ""
echo "To use with AmberTools, load both modules in your job script:"
echo "  ml conda"
echo "  ml ambertools/25"
echo "  source activate $ENV_PREFIX"

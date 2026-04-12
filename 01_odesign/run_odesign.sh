#!/bin/bash
#
# run_odesign.sh
#
# SLURM submission script for running ODesign via an apptainer container
# on Scholar GPU nodes. Generates cyclic peptide binders for a target
# protein specified in an input JSON file.
#
# Usage:
#   sbatch run_odesign.sh <input.json>
#
# Example:
#   sbatch run_odesign.sh mdm2_cyclic.json
#
# Prerequisites:
#   - Apptainer container built at $ODESIGN_SIF (edit path below).
#     The SIF must have the protenix fastfold_layer_norm_cuda extension
#     pre-built during container build (see odesign.def %post section).
#   - Model checkpoints downloaded to $CKPT_DIR
#   - Inference data (components.v20240608.cif[.rdkit_mol.pkl]) in $DATA_DIR

#SBATCH -A gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --gpus-per-node=1
#SBATCH --constraint=J
#SBATCH --time=00:30:00
#SBATCH --job-name odesign
#SBATCH --output=%x-%J-%u.out
#SBATCH --error=%x-%J-%u.err

set -euo pipefail

# Purge all loaded modules so that host compilers, CUDA toolkits, or
# other library paths from Scholar's spack modules cannot leak into
# the apptainer container. The container provides its own toolchain.
module purge 2>/dev/null || true

# -----------------------------------------------------------------------
# Paths — EDIT THESE to match your Scholar setup
# -----------------------------------------------------------------------
ODESIGN_SIF="/class/bsdrown/apps/odesign/odesign.sif"
CKPT_DIR="/class/bsdrown/apps/odesign/ckpt"
DATA_DIR="/class/bsdrown/apps/odesign/data"

WORKDIR="$SCRATCH/odesign"
INPUT_DIR="$WORKDIR/inputs"
OUTPUT_DIR="$WORKDIR/outputs"

# Writable cache directories for Triton autotune and any HuggingFace
# downloads. These need to be on scratch because home quotas are small
# and the container's internal /root is read-only.
CACHE_DIR="$WORKDIR/cache"
mkdir -p "$OUTPUT_DIR" "$CACHE_DIR/triton" "$CACHE_DIR/huggingface"

# -----------------------------------------------------------------------
# Determine input JSON
# -----------------------------------------------------------------------
if [[ -z "${1:-}" ]]; then
    echo "Error: must provide an input JSON filename as argument."
    echo "Usage: sbatch run_odesign.sh <input.json>"
    exit 1
fi

INPUT_JSON="$INPUT_DIR/$1"
if [[ ! -f "$INPUT_JSON" ]]; then
    echo "Error: input file '$INPUT_JSON' not found."
    exit 1
fi

EXP_NAME="${1%.json}"

# -----------------------------------------------------------------------
# ODesign inference parameters
# -----------------------------------------------------------------------
INFER_MODEL="odesign_base_prot_flex"   # flexible-receptor protein/peptide design
DESIGN_MODALITY="protein"              # cyclic peptide counts as protein modality
SEEDS="[42,123,456]"                   # three seeds for diversity
N_SAMPLE=5                             # 5 samples per seed = 15 total designs
USE_MSA="false"
NUM_WORKERS=0                          # avoid DataLoader worker OOM

# -----------------------------------------------------------------------
# Report environment
# -----------------------------------------------------------------------
echo "=== ODesign Cyclic Peptide Design ==="
echo "Date:      $(date)"
echo "Node:      $(hostname)"
echo "GPU:       $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "Container: $ODESIGN_SIF"
echo "Input:     $INPUT_JSON"
echo "Output:    $OUTPUT_DIR/$EXP_NAME"
echo "Model:     $INFER_MODEL"
echo ""

# -----------------------------------------------------------------------
# Run ODesign inside the apptainer container
# -----------------------------------------------------------------------
# Bind mounts:
#   /app/ODesign/ckpt     — model checkpoints (read-only)
#   /app/ODesign/data     — inference data components (read-only)
#   /app/ODesign/outputs  — where CIF files are written
#   /app/ODesign/inputs   — input JSON location (read-only)
#   /cache                — writable cache for Triton and HuggingFace
#
# --nv enables NVIDIA GPU passthrough
# --cleanenv strips host environment (prevents host compiler/lib leakage)
# --writable-tmpfs lets the container write to /tmp without modifying the SIF

# Use APPTAINERENV_ prefix to pass env vars into the container.
export APPTAINERENV_TRITON_CACHE_DIR=/cache/triton
export APPTAINERENV_HF_HOME=/cache/huggingface

# Write the inference command to a temporary shell script rather than
# inlining it into bash -c. This avoids fragile quote-escaping through
# multiple shell expansion layers.
RUN_SCRIPT="$WORKDIR/run_inference_$SLURM_JOB_ID.sh"
cat > "$RUN_SCRIPT" << SCRIPT_EOF
#!/bin/bash
set -e
cd /app/ODesign

python ./scripts/inference.py \\
    exp=train_${INFER_MODEL} \\
    data_root_dir=./data \\
    ckpt_root_dir=./ckpt \\
    exp.infer_model_name=${INFER_MODEL} \\
    exp.design_modality=${DESIGN_MODALITY} \\
    exp.input_json_path=./inputs/$1 \\
    exp.exp_name=${EXP_NAME} \\
    exp.seeds='${SEEDS}' \\
    exp.model.sample_diffusion.N_sample=${N_SAMPLE} \\
    exp.use_msa=${USE_MSA} \\
    exp.num_workers=${NUM_WORKERS}
SCRIPT_EOF
chmod +x "$RUN_SCRIPT"

apptainer run --nv \
    --cleanenv \
    --writable-tmpfs \
    -B "$CKPT_DIR:/app/ODesign/ckpt:ro" \
    -B "$DATA_DIR:/app/ODesign/data:ro" \
    -B "$OUTPUT_DIR:/app/ODesign/outputs" \
    -B "$INPUT_DIR:/app/ODesign/inputs:ro" \
    -B "$CACHE_DIR:/cache" \
    -B "$RUN_SCRIPT:/tmp/run_inference.sh:ro" \
    "$ODESIGN_SIF" \
    bash /tmp/run_inference.sh

# Clean up the temporary script
rm -f "$RUN_SCRIPT"

echo ""
echo "=== ODesign complete ==="
echo "Results in: $OUTPUT_DIR/$EXP_NAME"
ls -la "$OUTPUT_DIR/$EXP_NAME" 2>/dev/null || echo "(output directory structure depends on ODesign version)"

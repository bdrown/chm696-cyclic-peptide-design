#!/bin/bash
#
# submit_design.sh
#
# Submit the MM-GBSA pipeline for a single design as a SLURM dependency
# chain, but only submit stages that actually need to run. Each stage's
# outputs are checked before submission; completed stages are skipped
# entirely (no job is created).
#
# Usage:
#   bash submit_design.sh <boltz_cif> <design_dir>
#
# Example:
#   bash submit_design.sh \
#       $SCRATCH/boltz/outputs/.../design_01_seed123_model_0.cif \
#       $SCRATCH/mmgbsa/design_01
#
# Behavior:
#   - All four stages already complete: nothing is submitted, exit 0.
#   - Some stages complete: the first incomplete stage and all downstream
#     stages are submitted, dependency-chained with afterok.
#   - Nothing complete: all four stages are submitted as a chain.
#   - Partially-done outputs (e.g. a truncated production.nc from a
#     killed MD run) are detected by the individual stage scripts when
#     they run, which will clear the bad output and re-run.

set -euo pipefail

if [[ -z "${1:-}" || -z "${2:-}" ]]; then
    echo "Usage: bash submit_design.sh <boltz_cif> <design_dir>"
    exit 1
fi

BOLTZ_CIF="$1"
DESIGN_DIR="$2"

if [[ ! -f "$BOLTZ_CIF" ]]; then
    echo "Error: CIF file '$BOLTZ_CIF' not found"
    exit 1
fi

mkdir -p "$DESIGN_DIR"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DESIGN_NAME="$(basename "$DESIGN_DIR")"

# Stage completion detection. These functions encode the same "done"
# criteria used inside each stage's SLURM wrapper script, but run here
# at submission time so we can skip creating jobs that would be no-ops.
stage1_done() {
    [[ -s "$DESIGN_DIR/complex.prmtop" \
       && -s "$DESIGN_DIR/complex.inpcrd" \
       && -s "$DESIGN_DIR/receptor.prmtop" \
       && -s "$DESIGN_DIR/receptor.inpcrd" \
       && -s "$DESIGN_DIR/ligand.prmtop" \
       && -s "$DESIGN_DIR/ligand.inpcrd" ]]
}

stage2_done() {
    [[ -s "$DESIGN_DIR/complex_min.rst7" \
       && -s "$DESIGN_DIR/complex_min.pdb" ]]
}

stage3_done() {
    if [[ ! -s "$DESIGN_DIR/production.nc" \
          || ! -s "$DESIGN_DIR/production.rst7" ]]; then
        return 1
    fi
    # Verify trajectory is a reasonable size (catches truncated files
    # from killed MD runs). A 10 ns, 500-frame trajectory of a ~2000
    # atom system is typically 15+ MB.
    local size
    size=$(stat -c %s "$DESIGN_DIR/production.nc")
    [[ "$size" -gt 1000000 ]]
}

stage4_done() {
    [[ -s "$DESIGN_DIR/mmgbsa.dat" ]] \
        && grep -q "DELTA TOTAL" "$DESIGN_DIR/mmgbsa.dat" 2>/dev/null
}

echo "=== Submitting pipeline for $DESIGN_NAME ==="
echo "CIF:        $BOLTZ_CIF"
echo "Design dir: $DESIGN_DIR"
echo ""

# Report current state before deciding what to submit
echo "Current state:"
stage1_done && echo "  Stage 1 (prep):     done" || echo "  Stage 1 (prep):     pending"
stage2_done && echo "  Stage 2 (minimize): done" || echo "  Stage 2 (minimize): pending"
stage3_done && echo "  Stage 3 (md):       done" || echo "  Stage 3 (md):       pending"
stage4_done && echo "  Stage 4 (mmgbsa):   done" || echo "  Stage 4 (mmgbsa):   pending"
echo ""

# If everything is done, exit successfully without submitting anything
if stage1_done && stage2_done && stage3_done && stage4_done; then
    echo "All stages complete. Nothing to submit."
    exit 0
fi

# Submit the pipeline starting from the first incomplete stage. The
# dependency chain uses afterok so downstream stages wait for upstream
# success. LAST_JID tracks the most recently submitted job for chaining.
LAST_JID=""

submit_stage() {
    local stage_name="$1"
    local script_name="$2"
    shift 2
    local args=("$@")

    local sbatch_args=(
        --parsable
        --job-name="${stage_name}_${DESIGN_NAME}"
    )
    if [[ -n "$LAST_JID" ]]; then
        sbatch_args+=(--dependency=afterok:$LAST_JID)
    fi

    local jid
    jid=$(sbatch "${sbatch_args[@]}" "$SCRIPT_DIR/$script_name" "${args[@]}")
    if [[ -n "$LAST_JID" ]]; then
        echo "  $stage_name: job $jid (after $LAST_JID)"
    else
        echo "  $stage_name: job $jid"
    fi
    LAST_JID="$jid"
}

echo "Submitting:"

if ! stage1_done; then
    submit_stage "prep" "run_prepare_amber.sh" "$BOLTZ_CIF" "$DESIGN_DIR"
fi

if ! stage2_done; then
    submit_stage "min" "run_minimize.sh" "$DESIGN_DIR"
fi

if ! stage3_done; then
    submit_stage "md" "run_md.sh" "$DESIGN_DIR"
fi

if ! stage4_done; then
    submit_stage "mmgbsa" "run_mmgbsa.sh" "$DESIGN_DIR"
fi

echo ""
echo "Monitor with:"
echo "  squeue -u \$USER"
echo "  bash $SCRIPT_DIR/check_status.sh $DESIGN_DIR"

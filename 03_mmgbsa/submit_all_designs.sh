#!/bin/bash
#
# submit_all_designs.sh
#
# Submit the MM-GBSA pipeline for all designs listed in the Boltz
# evaluation manifest. This is the normal entry point for the
# problem set.
#
# Usage:
#   bash submit_all_designs.sh <manifest_tsv> <mmgbsa_root>
#
# Example:
#   bash submit_all_designs.sh \
#       $SCRATCH/boltz/inputs/odesign_eval_manifest.tsv \
#       $SCRATCH/mmgbsa
#
# The manifest is produced by prepare_boltz_inputs.py. It has columns:
#   design_id  source_cif  peptide_sequence  boltz_yaml
# We use design_id to name the output directory and derive the Boltz
# prediction CIF path from where Boltz wrote its results.

set -euo pipefail

if [[ -z "${1:-}" || -z "${2:-}" ]]; then
    echo "Usage: bash submit_all_designs.sh <manifest_tsv> <mmgbsa_root>"
    echo ""
    echo "  manifest_tsv: output of prepare_boltz_inputs.py"
    echo "                (e.g. \$SCRATCH/boltz/inputs/odesign_eval_manifest.tsv)"
    echo "  mmgbsa_root:  directory where per-design subdirectories will be"
    echo "                created (e.g. \$SCRATCH/mmgbsa)"
    exit 1
fi

MANIFEST="$1"
MMGBSA_ROOT="$2"

if [[ ! -f "$MANIFEST" ]]; then
    echo "Error: manifest '$MANIFEST' not found"
    exit 1
fi

# Boltz output root — configurable but defaults to the conventional
# location set by run_boltz_evaluation.sh
BOLTZ_OUTPUT_ROOT="${BOLTZ_OUTPUT_ROOT:-$SCRATCH/boltz/outputs/odesign_eval}"
BOLTZ_PREDICTIONS="$BOLTZ_OUTPUT_ROOT/boltz_results_odesign_eval/predictions"

if [[ ! -d "$BOLTZ_PREDICTIONS" ]]; then
    echo "Error: Boltz predictions directory not found: $BOLTZ_PREDICTIONS"
    echo "Set BOLTZ_OUTPUT_ROOT in your environment if it's elsewhere."
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$MMGBSA_ROOT"

echo "=== Submitting MM-GBSA pipeline for all designs ==="
echo "Manifest:      $MANIFEST"
echo "Boltz results: $BOLTZ_PREDICTIONS"
echo "Output root:   $MMGBSA_ROOT"
echo ""

# Read manifest, skipping comment lines. Column 1 is design_id.
N_SUBMITTED=0
while IFS=$'\t' read -r DESIGN_ID REST; do
    # Skip comment/header lines
    [[ "$DESIGN_ID" =~ ^# ]] && continue
    [[ -z "$DESIGN_ID" ]] && continue

    # Find the model_0 CIF for this design
    CIF="$BOLTZ_PREDICTIONS/$DESIGN_ID/${DESIGN_ID}_model_0.cif"
    if [[ ! -f "$CIF" ]]; then
        echo "WARNING: Boltz prediction not found for $DESIGN_ID: $CIF"
        echo "  Skipping."
        continue
    fi

    DESIGN_DIR="$MMGBSA_ROOT/$DESIGN_ID"
    echo "--- $DESIGN_ID ---"
    bash "$SCRIPT_DIR/submit_design.sh" "$CIF" "$DESIGN_DIR"
    echo ""
    N_SUBMITTED=$((N_SUBMITTED + 1))
done < "$MANIFEST"

echo "=== Submitted $N_SUBMITTED design pipelines ==="
echo ""
echo "Monitor overall status with:"
echo "  squeue -u \$USER"
echo ""
echo "Check stage completion with:"
echo "  bash $SCRIPT_DIR/check_status.sh $MMGBSA_ROOT"

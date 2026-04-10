#!/bin/bash
#
# check_status.sh
#
# Report per-stage completion status for one or more design directories.
# Accepts either a single design dir or a root containing multiple
# design_* subdirectories.
#
# Status codes per stage:
#   done - stage output files exist and look valid
#   RUN  - a SLURM job for this stage is currently running
#   PEND - a SLURM job for this stage is queued but not yet running
#   FAIL - stage log exists but output files are missing (ran and failed)
#   --   - stage not yet started
#
# Usage:
#   bash check_status.sh <path>

set -euo pipefail

if [[ -z "${1:-}" ]]; then
    echo "Usage: bash check_status.sh <design_dir_or_root>"
    exit 1
fi

TARGET="$1"
if [[ ! -d "$TARGET" ]]; then
    echo "Error: '$TARGET' is not a directory"
    exit 1
fi

# Decide whether we're looking at one design or a root of many
if [[ -f "$TARGET/complex.prmtop" || -f "$TARGET/tleap.log" ]]; then
    DESIGN_DIRS=("$TARGET")
else
    DESIGN_DIRS=()
    while IFS= read -r -d '' d; do
        DESIGN_DIRS+=("$d")
    done < <(find "$TARGET" -maxdepth 1 -type d -name 'design_*' -print0 | sort -z)
fi

if [[ ${#DESIGN_DIRS[@]} -eq 0 ]]; then
    echo "No design directories found in $TARGET"
    exit 1
fi

# Snapshot the queue once up front to avoid spawning squeue per stage
# per design. Format: "JOBNAME STATE" per line.
# Valid STATE values from squeue include PD (pending), R (running),
# CG (completing), CF (configuring), S (suspended).
QUEUE_SNAPSHOT=$(squeue -u "$USER" --format="%j %T" --noheader 2>/dev/null || echo "")

# Look up the state of a specific job name in the snapshot.
# Returns "RUN", "PEND", "OTHER", or "" (not in queue).
queue_state_for_job() {
    local jobname="$1"
    local state
    state=$(echo "$QUEUE_SNAPSHOT" \
        | awk -v name="$jobname" '$1 == name { print $2; exit }')
    case "$state" in
        RUNNING)          echo "RUN" ;;
        PENDING)          echo "PEND" ;;
        "")               echo "" ;;
        *)                echo "OTHER" ;;
    esac
}

# Each check_X function takes ($design_dir, $design_name) and returns
# a 6-character status string. Stage naming follows the submit_design.sh
# convention: prep_<name>, min_<name>, md_<name>, mmgbsa_<name>.

check_stage() {
    local design_dir="$1"
    local design_name="$2"
    local job_prefix="$3"
    local done_condition_fn="$4"

    # 1. Is the stage queued or running?
    local job_name="${job_prefix}_${design_name}"
    local qstate
    qstate=$(queue_state_for_job "$job_name")
    if [[ "$qstate" == "RUN" ]]; then
        echo " RUN  "
        return
    fi
    if [[ "$qstate" == "PEND" ]]; then
        echo " PEND "
        return
    fi

    # 2. Is the stage done based on output files?
    if "$done_condition_fn" "$design_dir"; then
        echo " done "
        return
    fi

    # 3. Did it run previously and fail (log exists but outputs don't)?
    case "$job_prefix" in
        prep)
            [[ -f "$design_dir/tleap.log" ]] && echo " FAIL " && return ;;
        min)
            [[ -f "$design_dir/minimize.log" ]] && echo " FAIL " && return ;;
        md)
            [[ -f "$design_dir/md.log" ]] && echo " FAIL " && return ;;
        mmgbsa)
            [[ -f "$design_dir/mmgbsa.log" ]] && echo " FAIL " && return ;;
    esac

    # 4. Never started
    echo "  --  "
}

# Done-condition functions (must match stage output file expectations)
prep_done() {
    [[ -s "$1/complex.prmtop" && -s "$1/receptor.prmtop" && -s "$1/ligand.prmtop" ]]
}

min_done() {
    [[ -s "$1/complex_min.rst7" && -s "$1/complex_min.pdb" ]]
}

md_done() {
    if [[ ! -s "$1/production.nc" || ! -s "$1/production.rst7" ]]; then
        return 1
    fi
    local size
    size=$(stat -c %s "$1/production.nc")
    [[ "$size" -gt 1000000 ]]
}

mmgbsa_done() {
    [[ -s "$1/mmgbsa.dat" ]] && grep -q "DELTA TOTAL" "$1/mmgbsa.dat" 2>/dev/null
}

# Extract DELTA TOTAL from mmgbsa.dat if available
get_delta_total() {
    local d="$1"
    if [[ -s "$d/mmgbsa.dat" ]]; then
        awk '/^DELTA TOTAL/ { printf "%10.2f +/- %-6.2f\n", $3, $5; exit }' \
            "$d/mmgbsa.dat" 2>/dev/null || echo ""
    fi
}

# Print header
printf "%-22s  %-6s  %-6s  %-6s  %-6s  %s\n" \
    "Design" "Prep" "Min" "MD" "MMGB" "DELTA TOTAL (kcal/mol)"
printf "%-22s  %-6s  %-6s  %-6s  %-6s  %s\n" \
    "----------------------" "------" "------" "------" "------" \
    "----------------------"

for d in "${DESIGN_DIRS[@]}"; do
    name=$(basename "$d")
    s1=$(check_stage "$d" "$name" "prep"   prep_done)
    s2=$(check_stage "$d" "$name" "min"    min_done)
    s3=$(check_stage "$d" "$name" "md"     md_done)
    s4=$(check_stage "$d" "$name" "mmgbsa" mmgbsa_done)
    delta=$(get_delta_total "$d")

    printf "%-22s  %s  %s  %s  %s  %s\n" \
        "$name" "$s1" "$s2" "$s3" "$s4" "$delta"
done

echo ""
echo "Legend:  done = complete       RUN  = currently running"
echo "         PEND = queued         FAIL = log exists but no output"
echo "         --   = not yet started"

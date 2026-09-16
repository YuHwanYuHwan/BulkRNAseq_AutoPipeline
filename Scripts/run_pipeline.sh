#!/bin/bash
# run_pipeline.sh <group_dir> [<group_dir> ...]   - the whole thing, start to finish
#   FastQC -> Trimming -> Alignment -> strandedness probe -> HTSeq -> CPM -> MultiQC
#
#   The probe records an unambiguous strandedness itself and the run continues. Only a
#   borderline result stops it (exit 2): fill the value in group.conf and run the same
#   command again - every step skips what it already finished.
#
#   Given several groups it works through them in order, each one written to its own log
#   under logs/. A group that fails or stops at the probe does not hold up the rest; what
#   happened to each is listed at the end.
#
#   Submitted with sbatch on a cluster with more than one node, several groups are dealt out
#   one job per node and each node then works through its share in order. Run with bash it
#   never submits anything.
#
#   NODES  nodes to deal the groups out to. Default: every node sinfo reports.
#SBATCH --job-name=rnaseq
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=logs/rnaseq_%j.out
set -euo pipefail
[ $# -ge 1 ] || { sed -n "2,10p" "$0"; exit 1; }

STEPS=(FastQC Trimming Alignment probe_strandedness ReadCount CalcCPM MultiQC)

# sbatch copies this file to a spool directory, so its own path says nothing about where the
# repository is. The group directory does: it always lives under <repo>/rawData/.
repo_of() { local abs; abs="$(cd "$1" && pwd)"; [ "$abs" != "${abs%%/rawData/*}" ] && echo "${abs%%/rawData/*}"; }

# Every group is checked before the first one starts. These runs take hours, so a path typed
# wrong is worth hearing about now rather than after the ones before it have finished.
BAD=()
for g in "$@"; do
    [ -d "$g" ] || { BAD+=("$g is not a directory"); continue; }
    R="$(repo_of "$g")" || true
    [ -n "${R:-}" ] && [ -d "$R/Scripts" ] || BAD+=("$g is not under a pipeline rawData/ directory")
done
if [ ${#BAD[@]} -gt 0 ]; then
    printf '[ERROR] %s\n' "${BAD[@]}" >&2
    echo "        nothing was run" >&2
    exit 1
fi

# Two nodes finish two groups in the time one node finishes one, and separate nodes beat
# sharing one: counting is indifferent to company, but alignment is limited by memory
# bandwidth rather than cores, and two STAR processes on a node split that bandwidth.
#
# Only under sbatch, and only once. Started with bash this does nothing, so a foreground run
# stays a foreground run; the jobs it submits carry RNASEQ_LANE so they get on with the work
# instead of dealing the groups out again.
REPO="$(cd "$1" && pwd)"; REPO="${REPO%%/rawData/*}"
if [ $# -gt 1 ] && [ -n "${SLURM_JOB_ID:-}" ] && [ -z "${RNASEQ_LANE:-}" ] &&
   command -v sinfo >/dev/null 2>&1; then
    read -ra NODE_LIST <<< "${NODES:-$(sinfo -h -N -o '%N' | sort -u | tr '
' ' ')}"
    if [ ${#NODE_LIST[@]} -gt 1 ]; then
        # Dealt out in turn rather than cut into blocks: neighbouring groups tend to be the
        # ones most alike in size, so taking every Nth keeps the lanes closer in total.
        declare -A LANE=()
        i=0
        for g in "$@"; do
            n="${NODE_LIST[$(( i % ${#NODE_LIST[@]} ))]}"
            LANE[$n]="${LANE[$n]:-}${LANE[$n]:+ }$(cd "$g" && pwd)"
            i=$((i+1))
        done
        cd "$REPO"      # the job log path in the directives is relative to where sbatch ran
        for n in "${NODE_LIST[@]}"; do
            [ -n "${LANE[$n]:-}" ] || continue          # more nodes than groups
            echo "[NODE] $n <- $(tr ' ' '
' <<< "${LANE[$n]}" | sed 's|.*/rawData/||' | paste -sd' ')"
            # Unquoted on purpose: each group must arrive as its own argument.
            sbatch --export=ALL,RNASEQ_LANE=1 -w "$n" Scripts/run_pipeline.sh ${LANE[$n]}
        done
        exit 0
    fi
fi

# One group, start to finish. Each step is checked on its own: called from a context where
# errexit is suspended, a failing step would otherwise let the next one start on its output.
run_group() {   # $1=group_dir  $2=Scripts dir
    local G="$1" S="$2" step rc
    for step in "${STEPS[@]}"; do
        # Capture the status before anything else touches $? - inside `if ! cmd` it would
        # already be the inverted value, which is always 0.
        rc=0
        bash "$S/${step}.sh" "$G" || rc=$?
        [ "$rc" -eq 0 ] || return "$rc"
    done
}

STAMP="${SLURM_JOB_ID:-$(date +%Y%m%d_%H%M%S)}"
HELD=() FAILED=()

for g in "$@"; do
    REPO="$(repo_of "$g")"
    REL="$(cd "$g" && pwd)"; REL="${REL#*/rawData/}"
    LOG="${REPO}/logs/${REL//\//_}_${STAMP}.log"
    mkdir -p "$(dirname "$LOG")"

    echo "[GROUP] $g  -> $LOG"
    # tee, not a plain redirect: the groups run one after another, so the job log stays the
    # readable timeline it always was while each group also keeps its own copy.
    # errexit off around the pipeline, and PIPESTATUS read immediately: any command in
    # between, `|| true` included, would overwrite it with its own status.
    set +e
    run_group "$g" "${REPO}/Scripts" 2>&1 | tee "$LOG"
    rc=${PIPESTATUS[0]}
    set -e

    if [ "$rc" -eq 2 ]; then
        echo "[HOLD] $g - strandedness is yours to call; record it and run this same command again"
        HELD+=("$g")
    elif [ "$rc" -ne 0 ]; then
        echo "[FAIL] $g - rc=$rc, see $LOG" >&2
        FAILED+=("$g")
    fi
done

echo
echo "[DONE] $# group(s): $(( $# - ${#HELD[@]} - ${#FAILED[@]} )) finished, ${#HELD[@]} held, ${#FAILED[@]} failed"
[ ${#HELD[@]}   -eq 0 ] || printf '[HOLD] %s\n' "${HELD[@]}"
[ ${#FAILED[@]} -eq 0 ] || printf '[FAIL] %s\n' "${FAILED[@]}" >&2

# A failure outranks a hold: one wants fixing, the other is waiting for you.
[ ${#FAILED[@]} -eq 0 ] || exit 1
[ ${#HELD[@]}   -eq 0 ] || exit 2

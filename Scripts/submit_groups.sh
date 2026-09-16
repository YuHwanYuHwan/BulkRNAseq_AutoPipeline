#!/bin/bash
# submit_groups.sh <group_dir> [<group_dir> ...]
#   Spread groups over the cluster's nodes: one job per node, each working through its share
#   in order. Two nodes process two groups in the time one node processes one, and keeping
#   them on separate nodes matters because alignment is limited by memory bandwidth rather
#   than by cores.
#
#   bash Scripts/submit_groups.sh rawData/ProjectA/GroupA rawData/ProjectA/GroupB
#   bash Scripts/submit_groups.sh rawData/ProjectA/*/
#
#   NODES        nodes to use, space separated. Default: every node sinfo reports.
#   SBATCH_OPTS  extra sbatch arguments, e.g. SBATCH_OPTS="-c 16 --mem=48G"
set -euo pipefail

SUBMIT="${SUBMIT:-sbatch}"          # overridden by the self-check, which has no scheduler
[ $# -ge 1 ] || { sed -n '2,12p' "$0"; exit 1; }

# Group directories are passed on to sbatch as separate words, so a path with whitespace in it
# would arrive as two groups. Refusing is better than submitting something else than was asked.
# Not GROUPS: bash keeps that name for the caller's group ids, and assigning to it fails
# silently, taking the rest of its line with it.
GROUP_DIRS=()
for g in "$@"; do
    case "$g" in *[[:space:]]*) echo "[ERROR] path contains whitespace: $g" >&2; exit 1 ;; esac
    [ -d "$g" ] || { echo "[ERROR] not a directory: $g" >&2; exit 1; }
    abs="$(cd "$g" && pwd)"
    [ "$abs" != "${abs%%/rawData/*}" ] || { echo "[ERROR] not under a pipeline rawData/: $g" >&2; exit 1; }
    GROUP_DIRS+=("$abs")
done

# run_pipeline.sh writes its job log to a path relative to the submission directory, so submit
# from the repository root whatever directory this was called from. The group paths are already
# absolute by now, so the change of directory does not move them.
REPO="${GROUP_DIRS[0]%%/rawData/*}"
cd "$REPO"

if [ -z "${NODES:-}" ]; then
    command -v sinfo >/dev/null 2>&1 || { echo "[ERROR] no sinfo found - set NODES=\"node01 node02\"" >&2; exit 1; }
    NODES="$(sinfo -h -N -o '%N' | sort -u | tr '\n' ' ')"
fi
read -ra NODE_LIST <<< "$NODES"
[ ${#NODE_LIST[@]} -gt 0 ] || { echo "[ERROR] no nodes to submit to" >&2; exit 1; }

# Round robin rather than a block split: consecutive groups are usually the ones most alike in
# size, so dealing them out one at a time leaves the lanes closer in total than cutting the
# list in half would.
declare -A LANE=()
i=0
for g in "${GROUP_DIRS[@]}"; do
    n="${NODE_LIST[$(( i % ${#NODE_LIST[@]} ))]}"
    LANE[$n]="${LANE[$n]:-}${LANE[$n]:+ }$g"
    i=$((i+1))
done

for n in "${NODE_LIST[@]}"; do
    [ -n "${LANE[$n]:-}" ] || continue          # more nodes than groups
    echo "[NODE] $n <- $(tr ' ' '\n' <<< "${LANE[$n]}" | sed 's|.*/rawData/||' | paste -sd' ')"
    # Unquoted on purpose: each group must reach run_pipeline.sh as its own argument.
    $SUBMIT ${SBATCH_OPTS:-} -w "$n" Scripts/run_pipeline.sh ${LANE[$n]}
done

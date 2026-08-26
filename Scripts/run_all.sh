#!/bin/bash
# run_all.sh <group_dir>   - the whole pipeline in one go
#   FastQC -> Trimming -> Alignment -> strandedness probe -> HTSeq -> CPM -> MultiQC
#
#   Stops between the stages only when the probe cannot call the strandedness with
#   confidence (exit 2). Write the value into group.conf and run run_stage2.sh.
#SBATCH --job-name=rnaseq_all
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --output=logs/all_%j.out
set -euo pipefail
[ $# -eq 1 ] || { sed -n "2,5p" "$0"; exit 1; }
GROUP_ABS="$(cd "$1" && pwd)"
S="${GROUP_ABS%%/rawData/*}/Scripts"
[ -d "$S" ] || { echo "[ERROR] $1 is not under a pipeline rawData/ directory" >&2; exit 1; }

# Capture the status before anything else touches $? - inside `if ! cmd` it would already
# be the inverted value, which is always 0.
rc=0
bash "$S/run_stage1.sh" "$1" || rc=$?
if [ "$rc" -eq 2 ]; then
    echo "[HOLD] strandedness needs your call - see above, then: bash Scripts/run_stage2.sh $1"
    exit 2
elif [ "$rc" -ne 0 ]; then
    exit "$rc"
fi

bash "$S/run_stage2.sh" "$1"

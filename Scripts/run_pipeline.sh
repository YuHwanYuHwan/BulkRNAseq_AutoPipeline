#!/bin/bash
# run_pipeline.sh <group_dir>   - the whole thing, start to finish
#   FastQC -> Trimming -> Alignment -> strandedness probe -> HTSeq -> CPM -> MultiQC
#
#   The probe records an unambiguous strandedness itself and the run continues. Only a
#   borderline result stops it (exit 2): fill the value in group.conf and run the same
#   command again - every step skips what it already finished.
#SBATCH --job-name=rnaseq
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --output=logs/rnaseq_%j.out
set -euo pipefail
[ $# -eq 1 ] || { sed -n "2,7p" "$0"; exit 1; }

# sbatch copies this file to a spool directory, so its own path says nothing about where the
# repository is. The group directory does: it always lives under <repo>/rawData/.
GROUP_ABS="$(cd "$1" && pwd)"
S="${GROUP_ABS%%/rawData/*}/Scripts"
[ -d "$S" ] || { echo "[ERROR] $1 is not under a pipeline rawData/ directory" >&2; exit 1; }

bash "$S/FastQC.sh"    "$1"
bash "$S/Trimming.sh"  "$1"
bash "$S/Alignment.sh" "$1"

# Capture the status before anything else touches $? - inside `if ! cmd` it would already be
# the inverted value, which is always 0.
rc=0
bash "$S/probe_strandedness.sh" "$1" || rc=$?
if [ "$rc" -eq 2 ]; then
    echo "[HOLD] strandedness is yours to call - record it, then run this same command again"
    exit 2
elif [ "$rc" -ne 0 ]; then
    exit "$rc"
fi

bash "$S/ReadCount.sh" "$1"
bash "$S/CalcCPM.sh"   "$1"
bash "$S/MultiQC.sh"   "$1"

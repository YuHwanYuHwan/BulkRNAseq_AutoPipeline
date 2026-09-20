#!/bin/bash
# CalcCPM.sh <group_dir>
#   TMM-normalized CPM from the group count matrix.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"

[ $# -eq 1 ] || { sed -n '2,3p' "$0"; exit 1; }
init_group "$1"

MATRIX="${OUT_DIR}/${GROUP}_count_matrix.tsv"
[ -s "$MATRIX" ] || { echo "[ERROR] run ReadCount.sh first" >&2; exit 1; }

Rscript "$(dirname "${BASH_SOURCE[0]}")/CalcCPM.R" \
    "$MATRIX" "${OUT_DIR}/${GROUP}_logCPM.tsv"

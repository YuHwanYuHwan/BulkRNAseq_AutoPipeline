#!/bin/bash
# RefIndexing.sh <species> <overhang>
#   Build a STAR index for <species> at <overhang> if it does not exist yet.
#   Indexes are kept forever and shared across projects - build cost is paid once.
#
#   bash RefIndexing.sh Homo_sapiens 149
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"

[ $# -eq 2 ] || { sed -n '2,6p' "$0"; exit 1; }
SPECIES="$1"; OVERHANG="$2"

SP_DIR="${REF_ROOT}/${SPECIES}"
IDX_DIR="${SP_DIR}/index/overhang${OVERHANG}"

[ -d "$SP_DIR" ] || { echo "[ERROR] $SP_DIR not found. Download the genome first (see README)." >&2; exit 1; }
FA="$(ls "$SP_DIR"/*.dna.*.fa 2>/dev/null | head -1)"
GTF="$(ls "$SP_DIR"/*.gtf 2>/dev/null | head -1)"
[ -n "$FA" ]  || { echo "[ERROR] no genome FASTA in $SP_DIR" >&2; exit 1; }
[ -n "$GTF" ] || { echo "[ERROR] no GTF in $SP_DIR" >&2; exit 1; }

if [ -f "${IDX_DIR}/SAindex" ]; then
    echo "[IDX ] reuse $IDX_DIR"
    exit 0
fi

# An index the caller cannot read looks exactly like an index that is not there. Shared
# indexes are owned by whoever built them, so this is the ordinary case on a class machine,
# and rebuilding 29 GB per user is the wrong answer to it.
if [ -e "$IDX_DIR" ] && [ ! -r "$IDX_DIR" ]; then
    echo "[ERROR] $IDX_DIR exists but is not readable by $(id -un)." >&2
    echo "        Ask whoever owns it for: chmod -R a+rX $REF_ROOT" >&2
    exit 1
fi

# Guard: a partial index from a killed run would be silently reused as valid.
rm -rf "$IDX_DIR"
mkdir -p "$IDX_DIR"

SA_NBASES="$(sa_index_nbases "$(stat -c %s "$FA")")"

echo "[IDX ] building $SPECIES overhang=$OVERHANG genomeSAindexNbases=$SA_NBASES"
echo "       FASTA $(basename "$FA")"
echo "       GTF   $(basename "$GTF")"
STAR --runMode genomeGenerate \
    --genomeDir "$IDX_DIR" \
    --genomeFastaFiles "$FA" \
    --sjdbGTFfile "$GTF" \
    --sjdbOverhang "$OVERHANG" \
    --genomeSAindexNbases "$SA_NBASES" \
    --runThreadN "$THREADS"

[ -f "${IDX_DIR}/SAindex" ] || { echo "[ERROR] index build failed" >&2; rm -rf "$IDX_DIR"; exit 1; }

echo "[DONE] $IDX_DIR"

#!/bin/bash
# ReadCount.sh <group_dir>
#   HTSeq gene counts per sample, then a merged count matrix for the group.
#   This is the first step of stage 2: it refuses to run until strandedness is set,
#   because a wrong value silently halves the counts instead of raising an error.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"

[ $# -eq 1 ] || { sed -n '2,5p' "$0"; exit 1; }
init_group "$1"

: "${species:?group.conf must define species}"
case "${strandedness:-}" in
    no|yes|reverse) ;;
    *) echo "[ERROR] group.conf strandedness must be no|yes|reverse (got '${strandedness:-empty}')" >&2
       echo "        run: bash Scripts/probe_strandedness.sh ${GROUP_DIR}" >&2; exit 1 ;;
esac

ALIGN="${PROC_DIR}/Alignment_result"
OUT="${PROC_DIR}/HTseqCount_result"
mkdir -p "$OUT"
GTF="$(ls "${REF_ROOT}/${species}"/*.gtf 2>/dev/null | head -1 || true)"
[ -n "$GTF" ] || { echo "[ERROR] no GTF for $species in ${REF_ROOT}/${species}" >&2; exit 1; }

# Counting a sample already recorded under a different strandedness would be skipped by its
# .done marker, quietly mixing two settings in one matrix. Changing the value starts over.
MODE="${OUT}/.strandedness"
if [ -f "$MODE" ] && [ "$(cat "$MODE")" != "$strandedness" ]; then
    echo "[HTSEQ] strandedness changed: $(cat "$MODE") -> $strandedness, recounting the group"
    rm -f "$OUT"/.*.done "$OUT"/*.gene.counts
fi
printf '%s' "$strandedness" > "$MODE"

echo "[HTSEQ] strandedness=$strandedness  gtf=$(basename "$GTF")  parallel=$THREADS"

# htseq-count is single-threaded and CPU-bound - measured at ~23k read pairs per second
# whether one or four run at once. Samples are independent, so they go side by side.
count_one() {   # $1 = sample
    local sample="$1" cnt="${OUT}/${1}.gene.counts" bam="${ALIGN}/${1}/${1}Aligned.sortedByCoord.out.bam"
    [ -s "$bam" ] || { echo "[ERROR] missing BAM for $sample" >&2; return 1; }
    echo "[CNT ] $sample"
    htseq-count -r pos -s "$strandedness" "$bam" "$GTF" > "$cnt"
    # Wrong strandedness produces no error, only a collapsed assignment rate. Warn loudly.
    awk -F'\t' -v s="$sample" '
        /^__no_feature/ { nf=$2 } { t+=$2 }
        END { if (t>0 && 100*nf/t > 50)
                 printf "[WARN] %s: __no_feature %.1f%% - strandedness may be wrong\n", s, 100*nf/t }' "$cnt"
    mark_done "$OUT" "$sample"
}

n=0 skip=0
SAMPLES=()
while IFS=$'\t' read -r sample _ _; do
    SAMPLES+=("$sample")
    if is_done "$OUT" "$sample"; then echo "[SKIP] $sample"; skip=$((skip+1)); continue; fi
    while [ "$(jobs -rp | wc -l)" -ge "$THREADS" ]; do wait -n; done
    count_one "$sample" &
    n=$((n+1))
done < <(list_samples)
wait

# ── merge into one matrix ────────────────────────────────────────────────────
MATRIX="${OUT_DIR}/${GROUP}_count_matrix.tsv"

# Header starts with a leading "Gene" and uses leading tabs, never a trailing one -
# a trailing tab makes R read.table() add an empty column and turn everything character.
{ printf 'Gene'; printf '\t%s' "${SAMPLES[@]}"; printf '\n'; } > "$MATRIX"

# Gene column from the first sample, then one count column per sample.
CMD="paste <(grep -v '^__' ${OUT}/${SAMPLES[0]}.gene.counts | cut -f1)"
for s in "${SAMPLES[@]}"; do
    CMD="$CMD <(grep -v '^__' ${OUT}/${s}.gene.counts | cut -f2)"
done
eval "$CMD" >> "$MATRIX"

# The counter reports every sample in the group, not only the ones processed this time -
# a resumed run that skips all of them still has to read as success.
SKIPMSG=""
[ "$skip" -gt 0 ] && SKIPMSG="  ($skip already done)"
echo "[DONE] HTSeq $((n+skip)) samples$SKIPMSG"
echo "       matrix: $MATRIX  ($(( $(wc -l < "$MATRIX") - 1 )) genes x ${#SAMPLES[@]} samples)"

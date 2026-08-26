#!/bin/bash
# probe_strandedness.sh <group_dir>
#   Run HTSeq with -s reverse on ONE sample and read __no_feature.
#
#   A single reverse run separates all three cases:
#     __no_feature low   (~10-20%)  -> reverse
#     __no_feature ~50%             -> no        (unstranded: half the reads sit on the other strand)
#     __no_feature very high (80%+) -> yes       (forward: nearly nothing is assigned)
#
#   An unambiguous result is written into group.conf and the run continues. A result that
#   falls between those cases stops and waits for a person: a wrong strandedness raises no
#   error, it just quietly deflates every count, so a coin flip at the boundary is worse
#   than an interruption.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"

[ $# -eq 1 ] || { sed -n '2,9p' "$0"; exit 1; }
init_group "$1"

: "${species:?group.conf must define species}"
ALIGN="${PROC_DIR}/Alignment_result"
GTF="$(ls "${REF_ROOT}/${species}"/*.gtf 2>/dev/null | head -1)"
[ -n "$GTF" ] || { echo "[ERROR] no GTF for $species" >&2; exit 1; }

read -r SAMPLE _ _ < <(list_samples)
BAM="${ALIGN}/${SAMPLE}/${SAMPLE}Aligned.sortedByCoord.out.bam"
[ -s "$BAM" ] || { echo "[ERROR] no BAM for $SAMPLE - run Alignment.sh first" >&2; exit 1; }

echo "[PROBE] sample=$SAMPLE  -s reverse"
TMP="${PROC_DIR}/.strand_probe.counts"
htseq-count -r pos -s reverse "$BAM" "$GTF" > "$TMP"

REPORT="$(strand_report "$TMP")"
rm -f "$TMP"
echo "$REPORT"
VERDICT="$(awk '/likely strandedness/ { print $5 }' <<< "$REPORT")"

# Bands leave gaps on purpose. The old 35/65 cut points touched, so 34% and 36% got
# different answers off a difference that means nothing. Data sitting in a gap - a
# rRNA-depleted library, a degraded sample - is exactly what a person should look at.
FRAC=$(awk '/__no_feature/ { gsub(/[()%]/,"",$4); print $4 }' <<< "$REPORT")
AUTO=$(awk -v r="$FRAC" 'BEGIN {
    if      (r < 25)            print "reverse"
    else if (r >= 40 && r <= 60) print "no"
    else if (r > 75)            print "yes"
}')

# PROBE_AUTO=0 keeps the decision manual - the teaching setup uses that, because reading
# the number is the point of the exercise.
[ "${PROBE_AUTO:-1}" = 1 ] || AUTO=""

if [ -n "$AUTO" ]; then
    sed -i "s/^strandedness.*/strandedness = $AUTO/" "${GROUP_DIR}/group.conf"
    echo "  __no_feature ${FRAC}% is unambiguous -> strandedness = $AUTO, written to group.conf"
    echo
else
    cat <<MSG

  __no_feature ${FRAC}% falls between the three cases, so nothing was written.
  Read the numbers above and decide:

      sed -i 's/^strandedness.*/strandedness = ${VERDICT}/' ${GROUP_DIR}/group.conf

  Then run: bash Scripts/run_stage2.sh ${GROUP_DIR}
MSG
    exit 2        # run_all.sh stops here; a plain stage-1 run just ends
fi

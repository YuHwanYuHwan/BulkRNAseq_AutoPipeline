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
FRAC=$(awk -F'	' '/^__no_feature/ { nf=$2 } { t+=$2 } END { printf "%.1f", 100*nf/t }' "$TMP")
rm -f "$TMP"
echo "$REPORT"

AUTO=$(strand_call "$FRAC")

# Keep the measurement, not just its conclusion. MultiQC.sh puts it in the report, so a matrix
# read later still shows what the strandedness was decided on - and the numbers accumulate
# across datasets, which is the only way these bands stop being a guess.
echo "$FRAC" > "${PROC_DIR}/.strandprobe"

# PROBE_AUTO=0 keeps the decision manual - the teaching setup uses that, because reading
# the number is the point of the exercise. Which of the two reasons applies has to survive,
# or a clear-cut number gets reported as an ambiguous one.
MANUAL=0
[ "${PROBE_AUTO:-1}" = 1 ] || { MANUAL=1; AUTO=""; }

CONF="${GROUP_DIR}/group.conf"
if [ -n "$AUTO" ]; then
    # The line is normally present and empty. A hand-written conf that left it out would
    # make sed match nothing, and the value would silently never be recorded.
    if grep -q "^strandedness" "$CONF"; then
        sed -i "s/^strandedness.*/strandedness = $AUTO/" "$CONF"
    else
        echo "strandedness = $AUTO" >> "$CONF"
    fi
    echo "  __no_feature ${FRAC}% is unambiguous -> strandedness = $AUTO, written to group.conf"
    echo
else
    if [ "$MANUAL" = 1 ]; then
        WHY="automatic recording is off (PROBE_AUTO=0), so nothing was written."
    else
        WHY="__no_feature ${FRAC}% falls between the three cases, so nothing was written."
    fi
    cat <<MSG

  ${WHY}
  This is the call the pipeline will not make for you. Read the numbers above,
  probe another sample if it helps, then record one of:

      sed -i 's/^strandedness.*/strandedness = reverse/' ${CONF}
      sed -i 's/^strandedness.*/strandedness = no/'      ${CONF}
      sed -i 's/^strandedness.*/strandedness = yes/'     ${CONF}

  Then run the pipeline again - finished steps are skipped:

      bash Scripts/run_pipeline.sh ${GROUP_DIR}
MSG
    exit 2        # run_pipeline.sh stops here
fi

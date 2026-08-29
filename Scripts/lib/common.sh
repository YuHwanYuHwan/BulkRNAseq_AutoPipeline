# common.sh - shared helpers sourced by every step script.
# Contract: each step takes ONE group directory. Directory layout is the source of truth.
#   subdirectory present -> directory name = sample_id, FASTQs inside are merged
#   flat FASTQ files     -> filename stem (minus _1/_2) = sample_id

PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
[ -f "${PIPELINE_ROOT}/config.sh" ] && source "${PIPELINE_ROOT}/config.sh"
REF_ROOT="${PIPELINE_ROOT}/reference_Genomes"

# Tool paths: use PATH if available, else *_BIN from config.sh
if [ -n "${CONDA_ENV:-}" ] && ! command -v fastqc >/dev/null 2>&1 && command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate "$CONDA_ENV"
fi
FASTQC_BIN="${FASTQC_BIN:-fastqc}"

# Cores, most authoritative first. Inside a SLURM job the reservation wins over anything in
# config.sh - asking for 32 cores and then using 8 wastes the other 24, and using more than
# reserved fights the cgroup. Outside a job, config.sh wins, which is how you stay polite on a
# shared head node; nproc is the last resort.
THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-$(nproc)}}"

# Every step stamps its own start and end, so a stage log reads as a timeline and a slow step
# is obvious without instrumenting anything. The EXIT trap fires on failure too.
if [ -z "${NO_STEP_LOG:-}" ]; then
    STEP="$(basename "$0" .sh)"
    SECONDS=0
    printf '[%s] ==> %s\n' "$(date '+%F %T')" "$STEP"
    trap 'rc=$?; printf "[%s] <== %s  %02d:%02d:%02d%s\n" "$(date "+%F %T")" "$STEP" \
        $((SECONDS/3600)) $((SECONDS%3600/60)) $((SECONDS%60)) \
        "$(case $rc in 0) ;; 2) echo "  HOLD - waiting for you" ;; *) echo "  FAILED rc=$rc" ;; esac)"' EXIT
fi

# group dir -> project/group relative path. Processed/ and Output/ mirror the same layout.
init_group() {
    GROUP_DIR="$(cd "$1" && pwd)"
    REL="${GROUP_DIR#*/rawData/}"
    [ "$REL" != "$GROUP_DIR" ] || { echo "[ERROR] not under rawData/: $1" >&2; exit 1; }
    PROC_DIR="${PIPELINE_ROOT}/Processed/${REL}"
    OUT_DIR="${PIPELINE_ROOT}/Output/${REL}"
    GROUP="$(basename "$REL")"
    # group.conf is written by hand, so tolerate spaces around '=' and quotes.
    # Sourcing it directly would fail on "species = Homo_sapiens".
    CONF="${GROUP_DIR}/group.conf"
    if [ -f "$CONF" ]; then
        while IFS= read -r line; do
            line="${line%%#*}"
            [[ "$line" == *=* ]] || continue
            local k="${line%%=*}" v="${line#*=}"
            k="$(echo "$k" | tr -d '[:space:]')"
            v="$(echo "$v" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//; s/^"\(.*\)"$/\1/; s/^'"'"'\(.*\)'"'"'$/\1/')"
            [ -n "$k" ] && printf -v "$k" '%s' "$v"
        done < "$CONF"
    fi
    mkdir -p "$PROC_DIR" "$OUT_DIR"
}

# Sample list. One line per sample: sample_id<TAB>R1(comma-joined)<TAB>R2(comma-joined, empty if single-end)
list_samples() {
    local d f stem s
    for d in "$GROUP_DIR"/*/; do
        [ -d "$d" ] || continue
        s="$(basename "$d")"
        printf '%s\t%s\t%s\n' "$s" \
            "$(ls "$d"*_1.fastq* "$d"*_1.fq* 2>/dev/null | paste -sd,)" \
            "$(ls "$d"*_2.fastq* "$d"*_2.fq* 2>/dev/null | paste -sd,)"
    done
    for f in "$GROUP_DIR"/*.fastq.gz "$GROUP_DIR"/*.fastq "$GROUP_DIR"/*.fq.gz "$GROUP_DIR"/*.fq; do
        [ -e "$f" ] || continue
        stem="$(basename "$f")"; stem="${stem%%.*}"
        case "$stem" in
            *_2) continue ;;                  # R2 is handled together with its R1
            # Only a _1 stem can have a mate. Without this branch the substitution below is a
            # no-op and ls returns the file itself, pairing a single-end read with itself.
            *_1) printf '%s\t%s\t%s\n' "${stem%_1}" "$f" "$(ls "${f/_1./_2.}" 2>/dev/null || true)"
                 continue ;;
        esac
        printf '%s\t%s\t\n' "$stem" "$f"
    done
}

is_done()   { [ -f "${1}/.${2}.done" ]; }
mark_done() { touch "${1}/.${2}.done"; }

# Max read length over the first N reads of a FASTQ (line 2 of every 4).
# Trimming makes read length variable, so the MAX must be used - sjdbOverhang has to
# accommodate the longest read that can span a junction.
max_read_length() {   # $1=fastq(.gz)  [$2=reads to scan, default 10000]
    local f="$1" n="${2:-10000}" m
    # awk quits after n reads, which SIGPIPEs zcat. Under pipefail that 141 would propagate
    # and set -e would kill the caller, so the scan is run with pipefail off.
    set +o pipefail
    m=$(zcat -f "$f" 2>/dev/null | awk -v n="$n" 'NR%4==2 { if (length($0)>m) m=length($0); c++ }
                                                  c>=n { exit } END { print m+0 }')
    set -o pipefail
    echo "$m"
}

# Read an HTSeq counts file and report how the reads were assigned. Reporting only: the call
# belongs to probe_strandedness.sh, which owns the thresholds.
strand_report() {   # $1 = *.gene.counts
    awk -F'\t' '
        /^__/ { special[$1]=$2; tot+=$2; next }
        { assigned+=$2; tot+=$2 }
        END {
            nf = special["__no_feature"]+0
            printf "\n  total reads counted : %d\n", tot
            printf "  assigned to genes   : %d (%.1f%%)\n", assigned, 100*assigned/tot
            printf "  __no_feature        : %d (%.1f%%)\n", nf, 100*nf/tot
            printf "  __ambiguous         : %d\n", special["__ambiguous"]+0
        }' "$1"
}

# __no_feature fraction from a -s reverse run -> strandedness, or nothing when the number
# sits between the three cases. The bands leave gaps on purpose: touching cut points would
# make 34% and 36% different answers off a difference that means nothing, and a wrong
# strandedness raises no error - it quietly deflates every count. Near a boundary, stopping
# to ask beats guessing.
strand_call() {   # $1 = __no_feature percentage
    awk -v r="$1" 'BEGIN {
        if      (r < 25)             print "reverse"
        else if (r >= 40 && r <= 60) print "no"
        else if (r > 75)             print "yes"
    }'
}

# Bare version number for each tool. Every tool prints its version differently
# ("multiqc, version 1.35", "R version 4.5.3 (2026-03-11) -- ..."), so the first
# dotted number is taken and the caller decides how to present it.
tool_version() {
    local raw=""
    case "$1" in
        fastqc)   raw=$("$FASTQC_BIN" --version 2>&1) ;;
        cutadapt) raw=$(cutadapt --version 2>&1) ;;
        STAR)     raw=$(STAR --version 2>&1) ;;
        htseq)    raw=$(htseq-count --version 2>&1) ;;
        multiqc)  raw=$(multiqc --version 2>&1) ;;
        R)        raw=$(R --version 2>&1) ;;
        edgeR)    raw=$(Rscript -e 'cat(as.character(packageVersion("edgeR")))' 2>/dev/null) ;;
    esac
    grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?[a-z]?' <<< "$raw" | head -1
}

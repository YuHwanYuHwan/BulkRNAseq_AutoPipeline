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
#
# SLURM_CPUS_ON_NODE sits between them because SLURM_CPUS_PER_TASK is only set for a job that
# asked with --cpus-per-task. Ask with --ntasks instead and it is empty, and nproc would then
# report every core on the machine rather than the few the job was given: a one-core
# reservation quietly running 128 threads.
THREADS="${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-${THREADS:-$(nproc)}}}"

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

    # Two files that resolve to one sample id share one .done marker: the first is processed,
    # the rest are skipped as already finished, and the run reports success having dropped
    # them. Nothing downstream can tell that apart from a group that really held one sample, so
    # it has to stop here - every step script comes through init_group, and this is the only
    # place that can end the step rather than warn from inside a process substitution.
    local dup
    dup="$(list_samples 2>/dev/null | cut -f1 | sort | uniq -d | paste -sd' ')"
    if [ -n "$dup" ]; then
        echo "[ERROR] one sample id comes from more than one place in $1: $dup" >&2
        list_samples 2>/dev/null |
            awk -F'\t' -v d=" $dup " '
                index(d, " " $1 " ") { printf "          %s   R1=%s   R2=%s\n", $1, $2, $3 }' >&2
        echo "        Rename them so each sample has one id, or put one sample's runs in a" >&2
        echo "        folder named after that sample." >&2
        exit 1
    fi
}

# FASTQ filename -> the name with only its extension removed. Cutting at the first dot
# instead turned Sample.L001_1.fastq.gz into "Sample", and so did every other lane of that
# sample: one .done marker between them, all but the first skipped as already finished, and
# the run reporting success having dropped three quarters of the data.
fq_stem() {
    local b="${1##*/}"
    b="${b%.gz}"; b="${b%.fastq}"; b="${b%.fq}"
    printf '%s' "$b"
}

# stem -> "sample_id<TAB>1|2", nothing at all when the name carries no read number.
# Recognised: _1 _2 _R1 _R2, each optionally followed by a block of digits, which is what
# bcl2fastq appends - Sample_S1_L001_R1_001.fastq.gz is the name a sequencing run arrives
# under, and it used to be read as a single-end sample called Sample_S1_L001_R1_001.
#
# Two patterns rather than one with an optional tail: in a single pattern the first group grows
# as far as it can and reads Sample_R1_00 plus read 1 out of Sample_R1_001.
#
# The list is closed on purpose. Pairing any two names that differ by a 1 and a 2 would also
# marry Patient1.fastq.gz to Patient2.fastq.gz - two people, one sample.
read_tag() {
    [[ "$1" =~ ^(.+)_R?([12])$ ]] || [[ "$1" =~ ^(.+)_R?([12])_[0-9]+$ ]] || return 1
    printf '%s\t%s' "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}"
}

# Sample list. One line per sample: sample_id<TAB>R1(comma-joined)<TAB>R2(comma-joined, empty if single-end)
list_samples() {
    local d f s t id tag r1 r2
    local -A R1_OF R2_OF SEEN
    local -a ORDER

    # A subdirectory is one sample, and everything inside it is that sample's runs to merge.
    for d in "$GROUP_DIR"/*/; do
        [ -d "$d" ] || continue
        s="$(basename "$d")"
        r1=""; r2=""
        for f in "$d"*.fastq.gz "$d"*.fastq "$d"*.fq.gz "$d"*.fq; do
            [ -e "$f" ] || continue
            if t="$(read_tag "$(fq_stem "$f")")"; then
                if [ "${t#*$'\t'}" = 2 ]; then r2="${r2:+$r2,}$f"; else r1="${r1:+$r1,}$f"; fi
            else
                # No read number anywhere in the name: a single-end run. A folder of _R1/_R2
                # pairs used to land here for every file, be concatenated into one stream, and
                # have every pair counted twice.
                r1="${r1:+$r1,}$f"
            fi
        done
        [ -n "$r1" ] || { echo "[WARN] $s holds no FASTQ - skipped" >&2; continue; }
        printf '%s\t%s\t%s\n' "$s" "$r1" "$r2"
    done

    # Flat files: the name says which sample and which read. Collected before anything is
    # printed, so an R2 finds its R1 whatever order the shell hands the globs over in.
    for f in "$GROUP_DIR"/*.fastq.gz "$GROUP_DIR"/*.fastq "$GROUP_DIR"/*.fq.gz "$GROUP_DIR"/*.fq; do
        [ -e "$f" ] || continue
        if t="$(read_tag "$(fq_stem "$f")")"; then
            id="${t%$'\t'*}"; tag="${t#*$'\t'}"
        else
            id="$(fq_stem "$f")"; tag=1
        fi
        if [ "$tag" = 2 ]; then R2_OF[$id]="$f"; else R1_OF[$id]="$f"; fi
        [ -n "${SEEN[$id]:-}" ] || { SEEN[$id]=1; ORDER+=("$id"); }
    done

    for id in ${ORDER[@]+"${ORDER[@]}"}; do
        if [ -n "${R1_OF[$id]:-}" ]; then
            printf '%s\t%s\t%s\n' "$id" "${R1_OF[$id]}" "${R2_OF[$id]:-}"
        else
            # Half a pair arrived. Said out loud and kept under its own file name rather than
            # dropped: a sample that simply vanishes from the list leaves nobody anything to
            # look for, and it was vanishing without a word.
            s="$(fq_stem "${R2_OF[$id]}")"
            echo "[WARN] $s: no R1 mate found - treated as single-end" >&2
            printf '%s\t%s\t\n' "$s" "${R2_OF[$id]}"
        fi
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

# STAR's suffix-array index has a size parameter that must come down for a small genome.
# The manual gives min(14, log2(GenomeLength)/2 - 1); the default 14 suits a mammal and on a
# 120 Mb plant genome it wastes memory and draws a warning. The FASTA's size on disk stands in
# for the genome length: newlines and headers inflate it by a couple of percent, which a
# base-2 logarithm halved and floored does not notice.
sa_index_nbases() {   # $1 = genome FASTA size in bytes
    awk -v n="$1" 'BEGIN {
        if (n < 2) n = 2
        v = int(log(n)/log(2)/2 - 1)
        if (v > 14) v = 14
        if (v < 1)  v = 1
        print v
    }'
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

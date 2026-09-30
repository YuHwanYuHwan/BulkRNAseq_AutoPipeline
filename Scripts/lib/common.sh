# common.sh - shared helpers sourced by every step script.
# Contract: each step takes ONE group directory. Directory layout is the source of truth.
#   subdirectory present -> directory name = sample_id, FASTQs inside are merged
#   flat FASTQ files     -> the name minus its extension and its read number = sample_id

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
    # LIST_ONLY means a person is looking at the group rather than a step processing it, and
    # list_samples.sh sets it. It has to print the group before saying what is wrong with it,
    # since these messages are what send you there in the first place.
    [ -z "${LIST_ONLY:-}" ] || return 0

    local out
    # On the error path only, run it again with stderr showing: the happy path must not print
    # every [WARN] twice, once here and once in the step's own run.
    out="$(list_samples 2>/dev/null)" || { list_samples >/dev/null; exit 1; }
    check_sample_ids "$out" "$1" || {
        echo "            bash Scripts/list_samples.sh $1" >&2
        exit 1
    }
}

# Sample ids must be unique: two files that resolve to one id share one .done marker, so the
# first is processed, the rest are skipped as already finished, and the run reports success
# having dropped them. Nothing downstream can tell that apart from a group that really held one
# sample. A function rather than four lines inside init_group so that list_samples.sh reaches
# the same verdict after printing the group, instead of passing it in silence.
check_sample_ids() {   # $1 = list_samples output, $2 = group dir
    local dup
    dup="$(cut -f1 <<< "$1" | sort | uniq -d | paste -sd' ')"
    [ -n "$dup" ] || return 0
    echo "[ERROR] more than one file resolves to each of these sample ids in $2: $dup" >&2
    echo "        Rename them so every sample has an id of its own, or put one sample's runs" >&2
    echo "        in a folder named after it.  What collided:" >&2
    return 1
}

# Only the extension comes off. Cutting at the first dot instead merged every lane of
# Sample.L001_1.fastq.gz into one sample id, and one .done marker between them.
fq_stem() {
    local b="${1##*/}"
    b="${b%.gz}"; b="${b%.fastq}"; b="${b%.fq}"
    printf '%s' "$b"
}

# A stem split into the sample it belongs to and which read of the pair it is.
# Recognised: _1 _2 _R1 _R2, each optionally followed by a block of digits, which is what
# bcl2fastq appends - Sample_S1_L001_R1_001.fastq.gz is the name a sequencing run arrives
# under, and it used to be read as a single-end sample called Sample_S1_L001_R1_001.
#
# Two patterns rather than one with an optional tail: in a single pattern the first group grows
# as far as it can and reads Sample_R1_00 plus read 1 out of Sample_R1_001.
#
# The list is closed on purpose. Pairing any two names that differ by a 1 and a 2 would also
# marry Patient1.fastq.gz to Patient2.fastq.gz - two people, one sample.
# Two files claiming the same reads of the same sample. Chunked deliveries land here too -
# bcl2fastq splits a large read into _R1_001, _R1_002 - and the way out is the same for both:
# a sample whose reads come in several files is what a merge folder is for.
claim_clash() {   # $1 = what is claimed twice, $2 and $3 = the two files
    echo "[ERROR] sample $1 is claimed by two files:" >&2
    echo "        $2" >&2
    echo "        $3" >&2
    echo "        Keep one file per read, or put a sample's several files in a folder named" >&2
    echo "        after that sample, where they are merged." >&2
}

read_tag() {   # stem -> RT_ID and RT_TAG; RT_TAG is empty when the name has no read number
    RT_ID="$1"; RT_TAG=
    [[ "$1" =~ ^(.+)_R?([12])$ ]] || [[ "$1" =~ ^(.+)_R?([12])_[0-9]+$ ]] || return 0
    RT_ID="${BASH_REMATCH[1]}"; RT_TAG="${BASH_REMATCH[2]}"
}

# Sample list. One line per sample: sample_id<TAB>R1(comma-joined)<TAB>R2(comma-joined, empty if single-end)
list_samples() {
    local d f s id r1 r2 bare RT_ID RT_TAG rc=0 aside=0
    local -A R1_OF R2_OF BARE_OF
    local -a ORDER

    # A subdirectory is one sample, and everything inside it is that sample's runs to merge.
    for d in "$GROUP_DIR"/*/; do
        [ -d "$d" ] || continue
        s="$(basename "$d")"
        r1=""; r2=""; bare=""
        for f in "$d"*.fastq.gz "$d"*.fastq "$d"*.fq.gz "$d"*.fq; do
            [ -e "$f" ] || continue
            read_tag "$(fq_stem "$f")"
            case "$RT_TAG" in
                1) r1="${r1:+$r1,}$f" ;;
                2) r2="${r2:+$r2,}$f" ;;
                *) bare="${bare:+$bare,}$f" ;;
            esac
        done
        # What a name with no read number means depends on the rest of the folder, so it is
        # decided only now. On its own it is a single-end run. Beside a pair it is the unmated
        # reads fasterq-dump --split-3 writes out, and concatenating those into R1 would leave
        # R1 longer than R2 - cutadapt stops on that, hours in, with a read-count mismatch.
        [ -n "$r1$r2$bare" ] || { echo "[WARN] $s holds no FASTQ - skipped" >&2; continue; }
        if [ -z "$r1$r2" ]; then
            r1="$bare"
        elif [ -n "$bare" ]; then
            aside=$((aside+1))
            [ -z "${LIST_ONLY:-}" ] || echo "[WARN] $s: unmated reads not used: $bare" >&2
        fi
        # Only read 2 arrived. Saying the folder holds no FASTQ would be untrue and would drop
        # the sample; the reads are there, and what is missing is their mates.
        if [ -z "$r1" ]; then
            echo "[WARN] $s: no R1 - its read 2 files are treated as single-end" >&2
            r1="$r2"; r2=""
        fi
        printf '%s\t%s\t%s\n' "$s" "$r1" "$r2"
    done

    # Flat files: the name says which sample and which read. Collected before anything is
    # printed, so an R2 finds its R1 whatever order the shell hands the globs over in.
    for f in "$GROUP_DIR"/*.fastq.gz "$GROUP_DIR"/*.fastq "$GROUP_DIR"/*.fq.gz "$GROUP_DIR"/*.fq; do
        [ -e "$f" ] || continue
        read_tag "$(fq_stem "$f")"
        [ -n "${R1_OF[$RT_ID]:-}${R2_OF[$RT_ID]:-}${BARE_OF[$RT_ID]:-}" ] || ORDER+=("$RT_ID")
        # One file per sample per read. A second claimant of the same read - D_1.fastq.gz
        # beside D_R1_001.fastq.gz - used to overwrite the first without a word, which loses a
        # file exactly as silently as two samples sharing a .done marker, so it stops the run
        # too. Refused in init_group rather than here, because a step reads this list through a
        # process substitution where an exit is the subshell's and looks like no samples.
        #
        # A name with no read number is held aside instead: whether it is a sample of its own
        # or the unmated half of a --split-3 trio depends on what else turns up.
        case "$RT_TAG" in
            1) if [ -n "${R1_OF[$RT_ID]:-}" ]; then
                   claim_clash "$RT_ID read 1" "${R1_OF[$RT_ID]}" "$f"; rc=1
               else R1_OF[$RT_ID]="$f"; fi ;;
            2) if [ -n "${R2_OF[$RT_ID]:-}" ]; then
                   claim_clash "$RT_ID read 2" "${R2_OF[$RT_ID]}" "$f"; rc=1
               else R2_OF[$RT_ID]="$f"; fi ;;
            *) if [ -n "${BARE_OF[$RT_ID]:-}" ]; then
                   claim_clash "$RT_ID" "${BARE_OF[$RT_ID]}" "$f"; rc=1
               else BARE_OF[$RT_ID]="$f"; fi ;;
        esac
    done

    for id in ${ORDER[@]+"${ORDER[@]}"}; do
        r1="${R1_OF[$id]:-}"; r2="${R2_OF[$id]:-}"
        # SRR001_1, SRR001_2 and a bare SRR001 is one run, not two samples: --split-3 puts the
        # reads whose mate is missing in the third file. They cannot go through cutadapt's
        # paired mode with the rest and every paired-end analysis drops them, but reads going
        # unused is worth a line rather than a silence.
        if [ -n "$r1" ] && [ -n "${BARE_OF[$id]:-}" ]; then
            aside=$((aside+1))
            [ -z "${LIST_ONLY:-}" ] || echo "[WARN] $id: unmated reads not used: ${BARE_OF[$id]}" >&2
        fi
        [ -n "$r1" ] || r1="${BARE_OF[$id]:-}"
        if [ -n "$r1" ]; then
            printf '%s\t%s\t%s\n' "$id" "$r1" "$r2"
        else
            # Half a pair arrived. Said out loud and kept under its own file name rather than
            # dropped: a sample that simply vanishes from the list leaves nobody anything to
            # look for, and it was vanishing without a word.
            s="$(fq_stem "$r2")"
            echo "[WARN] $s: no R1 mate found - treated as single-end" >&2
            printf '%s\t%s\t\n' "$s" "$r2"
        fi
    done

    # Every step calls this, so a line per sample would be the same forty lines eight times over
    # in one group's log. The count survives that; the detail is what list_samples.sh is for.
    [ "$aside" -eq 0 ] || [ -n "${LIST_ONLY:-}" ] ||
        echo "[WARN] $aside sample(s) hold unmated reads that go unused - see: bash Scripts/list_samples.sh $GROUP_DIR" >&2
    return "$rc"
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

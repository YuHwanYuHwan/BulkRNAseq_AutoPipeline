#!/bin/bash
# selfcheck.sh - logic self-check. Runs with no bioinformatics tool installed, so a fresh
# clone can be verified before spending hours on a real dataset. Called by setup.sh.
set -uo pipefail

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
ROOT="$TMP/repo"; mkdir -p "$ROOT/Scripts/lib"
LIB="$(dirname "${BASH_SOURCE[0]}")"
cp "$LIB/common.sh" "$ROOT/Scripts/lib/"
cp "$LIB/../PublicData_download.sh" "$ROOT/Scripts/"
NO_STEP_LOG=1                                 # keep step timestamps out of the report
source "$ROOT/Scripts/lib/common.sh"          # PIPELINE_ROOT now points at the fake repo

# The step scripts the pipeline calls, stubbed out by the tests that drive run_pipeline.sh
STEPS_STUB=(FastQC Trimming Alignment probe_strandedness ReadCount CalcCPM MultiQC)

# Flat file = one sample, subdirectory = one sample whose runs are merged.
t_list_samples() {
    local G="$ROOT/rawData/P/g1"; mkdir -p "$G/GSM_B"
    touch "$G/SRR001_1.fastq.gz" "$G/SRR001_2.fastq.gz" "$G/SRR009.fastq.gz" \
          "$G/GSM_B/SRR002_1.fastq.gz" "$G/GSM_B/SRR002_2.fastq.gz" \
          "$G/GSM_B/SRR003_1.fastq.gz" "$G/GSM_B/SRR003_2.fastq.gz"
    init_group "$G"
    local out; out=$(list_samples | sort)
    [ "$(grep -c . <<< "$out")" -eq 3 ]                                  || { echo "sample count $(grep -c . <<< "$out") != 3"; return 1; }
    grep '^GSM_B' <<< "$out" | grep -q 'SRR002_1.*,.*SRR003_1'           || { echo "GSM_B runs not merged"; return 1; }
    grep '^SRR001' <<< "$out" | grep -q 'SRR001_2.fastq.gz'              || { echo "SRR001 R2 not paired"; return 1; }
    grep '^SRR009' <<< "$out" | awk -F'\t' '$3==""' | grep -q .          || { echo "SRR009 not single-end"; return 1; }
    [ "$PROC_DIR" = "$ROOT/Processed/P/g1" ]                             || { echo "PROC_DIR=$PROC_DIR"; return 1; }
}

# group.conf is hand-written: spaces around '=', quotes, comments and blanks must all parse.
t_groupconf() {
    local G="$ROOT/rawData/P/g2"; mkdir -p "$G"
    cat > "$G/group.conf" <<'CONF'
# comment line
species       = Homo_sapiens
adapter_kit="Illumina_TruSeq"      # trailing comment
strandedness  =

CONF
    init_group "$G"
    [ "${species:-}"     = "Homo_sapiens"    ] || { echo "species='${species:-}'"; return 1; }
    [ "${adapter_kit:-}" = "Illumina_TruSeq" ] || { echo "adapter_kit='${adapter_kit:-}'"; return 1; }
    [ -z "${strandedness:-}" ]                 || { echo "strandedness should be empty"; return 1; }
}

# Trimming makes read length variable - the MAX must win, not the first read or the mean.
t_overhang() {
    local FQ="$TMP/t.fastq" L LZ BIG
    { printf '@r1\n%s\n+\n%s\n' "$(printf 'A%.0s' {1..100})" "$(printf 'I%.0s' {1..100})"
      printf '@r2\n%s\n+\n%s\n' "$(printf 'A%.0s' {1..150})" "$(printf 'I%.0s' {1..150})"
      printf '@r3\n%s\n+\n%s\n' "$(printf 'A%.0s' {1..75})"  "$(printf 'I%.0s' {1..75})"; } > "$FQ"
    gzip -cf "$FQ" > "$FQ.gz"
    L=$(max_read_length "$FQ"); LZ=$(max_read_length "$FQ.gz")
    [ "$L" -eq 150 ] && [ "$LZ" -eq 150 ] || { echo "max_read_length plain=$L gz=$LZ (expected 150)"; return 1; }

    # A real FASTQ is far longer than the scan window, so awk exits early and SIGPIPEs zcat.
    # With pipefail that aborts the whole run - which is invisible on a three-read file.
    BIG="$TMP/big.fastq.gz"
    awk 'BEGIN { for (i=0;i<20000;i++) printf "@r%d\nACGTACGTAC\n+\nIIIIIIIIII\n", i }' | gzip -c > "$BIG"
    L=$(max_read_length "$BIG" 100) || { echo "scan aborted on a file longer than the window"; return 1; }
    [ "$L" -eq 10 ] || { echo "early-exit scan gave $L (expected 10)"; return 1; }
}

# runinfo SampleName decides which runs share a merge folder.
t_grouping() {
    local META="$TMP/metadata.csv" acc smp n
    printf 'Run,LibraryLayout,SampleName\nSRR001,PAIRED,GSM_A\nSRR002,PAIRED,GSM_B\nSRR003,PAIRED,GSM_B\n' > "$META"
    local -A SAMPLE_OF MULTI count
    while IFS=$'\t' read -r acc smp; do
        SAMPLE_OF[$acc]="$smp"
        n=$(( ${count[$smp]:-0} + 1 )); count[$smp]=$n
        [ "$n" -gt 1 ] && MULTI[$smp]=1
    done < <(awk -F, 'NR==1 { for (i=1;i<=NF;i++) if ($i=="SampleName") c=i; next }
                      c { print $1 "\t" $c }' "$META")
    [ "${SAMPLE_OF[SRR003]:-}" = "GSM_B" ] || { echo "SRR003 -> ${SAMPLE_OF[SRR003]:-}"; return 1; }
    [ -n "${MULTI[GSM_B]:-}" ]             || { echo "GSM_B not flagged multi-run"; return 1; }
    [ -z "${MULTI[GSM_A]:-}" ]             || { echo "GSM_A wrongly flagged multi-run"; return 1; }
}

# A trailing tab in the header makes R read.table() add a column and go character.
t_matrix() {
    local OUT="$TMP/counts" MATRIX="$TMP/m.tsv" CMD s; mkdir -p "$OUT"
    local SAMPLES=(S1 S2)
    printf 'G1\t10\nG2\t20\n__no_feature\t5\n' > "$OUT/S1.gene.counts"
    printf 'G1\t30\nG2\t40\n__no_feature\t7\n' > "$OUT/S2.gene.counts"
    { printf 'Gene'; printf '\t%s' "${SAMPLES[@]}"; printf '\n'; } > "$MATRIX"
    CMD="paste <(grep -v '^__' ${OUT}/${SAMPLES[0]}.gene.counts | cut -f1)"
    for s in "${SAMPLES[@]}"; do CMD="$CMD <(grep -v '^__' ${OUT}/${s}.gene.counts | cut -f2)"; done
    eval "$CMD" >> "$MATRIX"
    [ "$(head -1 "$MATRIX")" = "$(printf 'Gene\tS1\tS2')" ] || { echo "header: $(head -1 "$MATRIX")"; return 1; }
    grep -q '__no_feature' "$MATRIX"                        && { echo "__ row leaked into matrix"; return 1; }
    [ "$(sed -n '2p' "$MATRIX")" = "$(printf 'G1\t10\t30')" ] || { echo "row G1 misaligned"; return 1; }
}

# strand_report from common.sh on synthetic counts: does __no_feature map to the right call?
t_probe_verdict() {
    local v got
    #                  clear cases          gaps that must stop
    for v in "13 reverse" "50 no" "85 yes" "30 -" "68 -"; do
        set -- $v
        got=$(strand_call "$1"); got="${got:--}"
        [ "$got" = "$2" ] || { echo "${1}% -> $got (expected $2)"; return 1; }
    done
}

# Several groups at once means a typo in the last one must surface before the first download,
# not twelve hours into the night. Runs without any tool installed: the refusal is reached
# before prefetch is ever called.
t_multigroup_preflight() {
    local A="$ROOT/rawData/P/ga" B="$ROOT/rawData/P/gb" out
    mkdir -p "$A" "$B"
    printf 'SRR0000001\n' > "$A/accessions.csv"          # B deliberately has no list
    out=$(bash "$ROOT/Scripts/PublicData_download.sh" "$A" "$B" 2>&1) &&
        { echo "accepted a group with no accessions.csv"; return 1; }
    grep -q 'accessions.csv not found' <<< "$out" || { echo "unhelpful message: $out"; return 1; }
    [ ! -e "$A/.SRR0000001.done" ] || { echo "downloaded before checking every group"; return 1; }
    [ ! -e "$A/.runinfo.csv" ]     || { echo "hit the network before checking every group"; return 1; }
}

# Several groups in one run_pipeline call: a group that stops at the probe must not take the
# ones after it down with it, and each must keep its own log. Stub steps stand in for the real
# tools, so this says nothing about the tools and everything about the sequencing.
t_pipeline_multigroup() {
    local R="$TMP/repo2" step rc out
    mkdir -p "$R/Scripts" "$R/rawData/P/a" "$R/rawData/P/b"
    cp "$LIB/../run_pipeline.sh" "$R/Scripts/"
    for step in FastQC Trimming Alignment ReadCount CalcCPM MultiQC; do
        printf '#!/bin/bash\necho "ran %s on $1"\n' "$step" > "$R/Scripts/${step}.sh"
    done
    # group a holds at the probe, group b sails through
    printf '#!/bin/bash\necho "probed $1"\ncase "$1" in *P/a) exit 2 ;; esac\n' \
        > "$R/Scripts/probe_strandedness.sh"

    rc=0; out=$(bash "$R/Scripts/run_pipeline.sh" "$R/rawData/P/a" "$R/rawData/P/b" 2>&1) || rc=$?
    [ "$rc" -eq 2 ]                                  || { echo "held run exited $rc (expected 2)"; return 1; }
    grep -q 'ran MultiQC on .*P/b' <<< "$out"        || { echo "group b did not finish after a held"; return 1; }
    grep -q 'ran ReadCount on .*P/a' <<< "$out"      && { echo "group a counted despite the hold"; return 1; }
    [ -s "$R/logs/P_a_"*.log ] && [ -s "$R/logs/P_b_"*.log ] || { echo "per-group logs missing"; return 1; }
    grep -q 'P/b' "$R/logs/P_a_"*.log                && { echo "group b leaked into group a's log"; return 1; }

    # A path that is not under rawData/ has no pipeline to run, and saying so must come before
    # the groups that are fine have already spent hours.
    rc=0; out=$(bash "$R/Scripts/run_pipeline.sh" "$R/rawData/P/a" "$TMP" 2>&1) || rc=$?
    [ "$rc" -eq 1 ]                        || { echo "bad path exited $rc (expected 1)"; return 1; }
    grep -q 'nothing was run' <<< "$out"   || { echo "unhelpful message: $out"; return 1; }
}

# Dealing groups out to nodes. A stub sbatch and sinfo stand in for the scheduler, so this
# runs where there is none and checks who gets submitted where rather than the submission.
t_pipeline_spread() {
    local R="$TMP/repo3" out n step
    mkdir -p "$R/Scripts" "$R/bin"
    cp "$LIB/../run_pipeline.sh" "$R/Scripts/"
    for step in "${STEPS_STUB[@]}"; do
        printf '#!/bin/bash\necho "ran %s on $1"\n' "$step" > "$R/Scripts/${step}.sh"
    done
    printf '#!/bin/bash\nprintf "node01\\nnode02\\n"\n' > "$R/bin/sinfo"
    printf '#!/bin/bash\necho "SBATCH $*"\n'            > "$R/bin/sbatch"
    chmod +x "$R/bin/sinfo" "$R/bin/sbatch"
    for n in a b c d e; do mkdir -p "$R/rawData/P/$n"; done

    # Under sbatch with several groups and more than one node: deal them out, run nothing here.
    out=$(cd "$R" && PATH="$R/bin:$PATH" SLURM_JOB_ID=1 bash Scripts/run_pipeline.sh \
          rawData/P/a rawData/P/b rawData/P/c rawData/P/d rawData/P/e 2>&1) \
        || { echo "spread failed: $out"; return 1; }
    grep -qE 'SBATCH .* -w node01 Scripts/run_pipeline.sh .*/P/a .*/P/c .*/P/e$' <<< "$out" ||
        { echo "node01 lane wrong: $out"; return 1; }
    grep -qE 'SBATCH .* -w node02 Scripts/run_pipeline.sh .*/P/b .*/P/d$' <<< "$out" ||
        { echo "node02 lane wrong: $out"; return 1; }
    [ "$(grep -c '^SBATCH ' <<< "$out")" -eq 2 ] || { echo "expected 2 submissions: $out"; return 1; }
    grep -q 'RNASEQ_LANE=1' <<< "$out"           || { echo "lanes not marked, they will spread again"; return 1; }
    grep -q 'ran FastQC' <<< "$out"              && { echo "the dispatching job also ran the steps"; return 1; }

    # A lane must get on with the work rather than deal the groups out a second time.
    out=$(cd "$R" && PATH="$R/bin:$PATH" SLURM_JOB_ID=2 RNASEQ_LANE=1 \
          bash Scripts/run_pipeline.sh rawData/P/a rawData/P/b 2>&1) || true
    grep -q '^SBATCH ' <<< "$out" && { echo "a lane submitted more jobs"; return 1; }
    grep -q 'ran MultiQC on .*P/b' <<< "$out" || { echo "lane did not run the steps: $out"; return 1; }

    # One group has nothing to spread, and bash rather than sbatch must never submit anything.
    out=$(cd "$R" && PATH="$R/bin:$PATH" SLURM_JOB_ID=3 bash Scripts/run_pipeline.sh rawData/P/a 2>&1) || true
    grep -q '^SBATCH ' <<< "$out" && { echo "single group was spread"; return 1; }
    out=$(cd "$R" && PATH="$R/bin:$PATH" bash Scripts/run_pipeline.sh rawData/P/a rawData/P/b 2>&1) || true
    grep -q '^SBATCH ' <<< "$out" && { echo "a foreground run submitted jobs"; return 1; }
    grep -q 'ran MultiQC on .*P/b' <<< "$out" || { echo "foreground run did not run the steps: $out"; return 1; }
}

PASS=0; FAIL=0
for t in t_list_samples t_groupconf t_overhang t_grouping t_matrix t_probe_verdict \
         t_multigroup_preflight t_pipeline_multigroup t_pipeline_spread; do
    if msg=$("$t" 2>&1); then PASS=$((PASS+1))
    else FAIL=$((FAIL+1)); printf '%s: %s\n' "${t#t_}" "${msg:-failed}" >&2; fi
done
[ "$FAIL" = 0 ] || exit 1
echo "$PASS/$PASS"

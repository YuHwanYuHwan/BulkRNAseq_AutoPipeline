#!/bin/bash
# sra_to_fastq.sh <group_dir>
#   Turn the .sra archives prefetch left in a group into compressed FASTQ.
#
#   Downloading needs the internet and so has to happen where you are logged in. This does
#   not: it reads a local file and writes a local file, which is what a compute node is for.
#   PublicData_download.sh submits it once the downloading is done.
#
#   The directory is the work list, as everywhere else here: every <accession>/ folder holding
#   an archive is one run to convert. Re-running picks up whatever is left.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"

[ $# -eq 1 ] || { sed -n '2,10p' "$0"; exit 1; }
GROUP_DIR="$(cd "$1" && pwd)"

command -v pigz >/dev/null 2>&1 && ZIP="pigz -p $THREADS" || ZIP="gzip"

n=0 skip=0 failed=0
while IFS= read -r d; do
    acc="$(basename "$d")"; dest="$(dirname "$d")"
    # A finished run leaves its .done flag and its FASTQ; the archive is then only taking up
    # room. It is the same size as the reads it holds.
    if is_done "$dest" "$acc"; then
        rm -rf "$d"; skip=$((skip+1)); continue
    fi
    echo "[FQ  ] $acc"
    if fasterq-dump --split-3 --threads "$THREADS" --outdir "$dest" "$d" >/dev/null &&
       $ZIP -f "$dest"/"$acc"*.fastq; then
        rm -rf "$d"
        mark_done "$dest" "$acc"        # only now: the flag means reads on disk, not an archive
        n=$((n+1))
    else
        # Half-written FASTQ would look like a finished run to everything downstream. The
        # archive stays: it is the half that took the network, and a retry can convert it
        # without asking NCBI for the reads again.
        rm -f "$dest/${acc}"*.fastq
        echo "[FAIL] $acc" >&2
        failed=$((failed+1))
    fi
done < <(find "$GROUP_DIR" -mindepth 1 -maxdepth 2 -type d -name '[SED]RR[0-9]*' | sort)

SKIPMSG=""
[ "$skip" -gt 0 ] && SKIPMSG="  ($skip already done)"
echo "[DONE] $n run(s) converted in $GROUP_DIR$SKIPMSG"
[ "$failed" -eq 0 ] || { echo "[FAIL] $failed run(s) left; re-running converts what is missing" >&2; exit 1; }

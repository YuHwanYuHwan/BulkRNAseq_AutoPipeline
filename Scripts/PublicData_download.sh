#!/bin/bash
# PublicData_download.sh <group_dir> [<group_dir> ...]
#                        <group_dir> <SRR ...>
#                        <group_dir> <list file>
#   Downloads raw FASTQ for a list of SRR accessions into the group folder.
#   Runs sharing a SampleName go into a subfolder, which marks them as one sample to merge.
#   Sample metadata is collected afterwards by fetch_metadata.sh.
#
#   Several group directories can be given at once. Each is expected to hold its own
#   accessions.csv, and they are downloaded one after another, so a night of groups is one
#   command instead of one command per group waited out in turn.
#
#   bash PublicData_download.sh rawData/ProjectA/GroupA rawData/ProjectA/GroupB
#   bash PublicData_download.sh rawData/ProjectA/GroupA SRR0000001 SRR0000002
#   bash PublicData_download.sh rawData/ProjectA/GroupA rawData/ProjectA/GroupA/accessions.csv
#
#   The list file has no format: every SRR/ERR/DRR accession found anywhere in it is used,
#   duplicates dropped. A plain list, a comma-separated line, or a whole SraRunTable.csv
#   pasted in all work, so there is nothing to reformat before running this.
set -euo pipefail

RUNINFO_URL="https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc="
ENA_URL="https://www.ebi.ac.uk/ena/portal/api/search"
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"

[ $# -ge 1 ] || { sed -n '2,13p' "$0"; exit 1; }

# Everything here hangs off one request to NCBI, and that is the request that goes away: the
# endpoint has spent whole evenings not answering. ENA holds the same fields under different
# names, so the file is rebuilt from there and everything downstream reads it unchanged.
# Batched because the query goes in the URL, and a few hundred accessions do not fit in one.
ena_runinfo() {   # accessions on stdin, one per line -> runinfo-shaped CSV on stdout
    local batch q
    echo "Run,SampleName,BioSample,BioProject,ScientificName"
    while mapfile -t -n 50 batch && [ "${#batch[@]}" -gt 0 ]; do
        q="$(printf 'run_accession="%s" OR ' "${batch[@]}")"; q="${q% OR }"
        curl -s --max-time 120 --retry 2 -G "$ENA_URL"             --data-urlencode 'result=read_run'             --data-urlencode "query=$q"             --data-urlencode 'fields=run_accession,sample_alias,sample_accession,study_accession,scientific_name'             --data-urlencode 'format=tsv' |
        # An alias is what groups several runs into one sample. Without one the run is its own.
        awk -F'	' 'NR>1 && $1 != "" { print $1 "," ($2 == "" ? $1 : $2) "," $3 "," $4 "," $5 }'
    done
}

# Every accession found anywhere in a file, duplicates dropped. Scraping rather than parsing
# means a run table or a copied web page works without editing.
scrape_accessions() { grep -oE '[SED]RR[0-9]+' "$1" | awk '!seen[$0]++'; }

download_group() {   # $1=group_dir, rest=accessions
    local GROUP_DIR="$1"; shift
    local ACCS=("$@") acc sample dest n
    mkdir -p "$GROUP_DIR"

    # -- 1. runinfo ---------------------------------------------------------
    # runinfo is machinery, not the user's metadata: it exists only to group runs by sample.
    # The conditions (tissue, treatment, donor) are not in it - fetch_metadata.sh gets those.
    # The endpoint takes a comma-separated list, so this is one request rather than one per run.
    # Short deadline, one retry: this endpoint answers in seconds when it answers at all, so
    # waiting longer only postpones asking ENA, which is standing right there with the answer.
    local META="${GROUP_DIR}/.runinfo.csv"
    curl -sf --max-time 45 --retry 1 --retry-delay 5          "${RUNINFO_URL}$(IFS=,; echo "${ACCS[*]}")" > "$META" || true
    if [ ! -s "$META" ]; then
        echo "[WARN] NCBI runinfo did not answer - reading the same fields from ENA" >&2
        printf '%s
' "${ACCS[@]}" | ena_runinfo > "$META"
    fi
    [ "$(wc -l < "$META")" -gt 1 ] || { echo "[ERROR] runinfo lookup failed for $GROUP_DIR" >&2; return 1; }
    echo "[INFO] runinfo: $(( $(wc -l < "$META") - 1 )) runs"

    # -- 1b. group.conf -----------------------------------------------------
    # The species is in runinfo, so for public data there is nothing left for a person to type.
    # Written only when the file is absent - a conf you edited is never overwritten. A group
    # holding more than one organism gets nothing: a mistake to look at, not to guess past.
    local CONF="${GROUP_DIR}/group.conf" SPECIES
    if [ ! -f "$CONF" ]; then
        SPECIES=$(awk -F, 'NR==1 { for (i=1;i<=NF;i++) if ($i=="ScientificName") c=i; next }
                           c { print $c }' "$META" | sort -u | tr ' ' '_')
        if [ "$(wc -l <<< "$SPECIES")" -eq 1 ] && [ -n "$SPECIES" ]; then
            { echo "species      = $SPECIES"; echo "strandedness ="; } > "$CONF"
            echo "[CONF] species = $SPECIES"
            [ -d "${REF_ROOT}/${SPECIES}" ] ||
                echo "[WARN] no ${REF_ROOT}/${SPECIES} yet - download that genome before running the pipeline"
        else
            echo "[WARN] could not settle on one organism - write $CONF yourself"
        fi
    fi

    # -- 2. acc -> SampleName map, one pass. Samples with >1 run go into a subfolder.
    local -A SAMPLE_OF MULTI count
    while IFS=$'\t' read -r acc sample; do
        SAMPLE_OF[$acc]="$sample"
        n=$(( ${count[$sample]:-0} + 1 )); count[$sample]=$n
        if [ "$n" -gt 1 ]; then MULTI[$sample]=1; fi
    done < <(awk -F, '
        NR==1 { for (i=1;i<=NF;i++) if ($i=="SampleName") c=i; next }
        c { print $1 "\t" $c }' "$META")

    # -- 3. download --------------------------------------------------------
    for acc in "${ACCS[@]}"; do
        sample="${SAMPLE_OF[$acc]:-$acc}"
        dest="$GROUP_DIR"
        [ -n "${MULTI[$sample]:-}" ] && dest="${GROUP_DIR}/${sample}"   # multi-run -> merge folder
        mkdir -p "$dest"

        if is_done "$dest" "$acc"; then echo "[SKIP] $acc"; continue; fi
        echo "[GET ] $acc -> $dest"
        # Downloading only: turning the archive into FASTQ is local work, and sra_to_fastq.sh
        # does it on a compute node once everything is here.
        if ! prefetch --output-directory "$dest" "$acc" >/dev/null; then
            # Over a night of downloads one expired link must cost that run, not the other
            # forty. Nothing half-downloaded is left behind, so re-running picks up exactly
            # what failed; the .done flag belongs to the conversion, where the reads appear.
            rm -rf "${dest:?}/${acc}"
            echo "[FAIL] $acc" >&2
            FAILED+=("${GROUP_DIR}/${acc}")
        fi
    done
}

# -- argument dispatch ------------------------------------------------------
# All arguments existing directories -> a list of groups, each reading its own accessions.csv.
# Anything else is the single-group form, where the rest of the line is accessions or a file.
ALL_DIRS=1
for a in "$@"; do [ -d "$a" ] || { ALL_DIRS=0; break; }; done

# Not GROUPS: bash keeps that name for the caller's group ids, and assigning to it fails
# silently - taking the rest of the line down with it.
declare -a GROUP_DIRS=() FAILED=() META_FAILED=() MISSING=()
declare -A ACCS_OF=()
if [ "$ALL_DIRS" = 1 ]; then
    # Every list is read and checked before a single byte is downloaded. A typo in the last
    # group is worth knowing about now, not twelve hours from now.
    for g in "$@"; do
        if [ -s "${g}/accessions.csv" ]; then
            ACCS_OF[$g]="$(scrape_accessions "${g}/accessions.csv" | paste -sd' ')"
            [ -n "${ACCS_OF[$g]}" ] || MISSING+=("${g}/accessions.csv holds no accession")
            GROUP_DIRS+=("$g")
        elif [ -e "${g}/accessions.csv" ]; then
            MISSING+=("${g}/accessions.csv is empty")
        else
            # Naming the near misses saves the next person the minutes it took to notice that
            # accession.csv and accessions.csv are not the same file.
            # || true: finding nothing is the ordinary case, not a reason to stop
            near="$(ls "$g" 2>/dev/null | grep -i 'acc.*\.\(csv\|txt\|tsv\)$' | paste -sd' ' || true)"
            MISSING+=("${g}/accessions.csv not found${near:+ (the folder holds: $near)}")
        fi
    done
    if [ ${#MISSING[@]} -gt 0 ]; then
        printf '[ERROR] %s\n' "${MISSING[@]}" >&2
        echo "        nothing was downloaded" >&2
        exit 1
    fi
else
    [ $# -ge 2 ] || { sed -n '2,13p' "$0"; exit 1; }
    g="$1"; shift
    if [ $# -eq 1 ] && [ -f "$1" ]; then ACCS_OF[$g]="$(scrape_accessions "$1" | paste -sd' ')"
    else ACCS_OF[$g]="$*"; fi
    [ -n "${ACCS_OF[$g]}" ] || { echo "[ERROR] no accession given" >&2; exit 1; }
    GROUP_DIRS+=("$g")
fi

TOTAL=0
for g in "${GROUP_DIRS[@]}"; do
    echo "[GROUP] $g"
    read -ra accs <<< "${ACCS_OF[$g]}"
    TOTAL=$((TOTAL + ${#accs[@]}))
    download_group "$g" "${accs[@]}" || FAILED+=("$g (runinfo)")
done

# Converting is local work: no internet, a lot of CPU, and hours of it. Off to the scheduler
# it goes, one job per group so they run alongside each other rather than one after another.
# --wrap rather than the script itself: sbatch copies a submitted script into a spool
# directory, where lib/common.sh is not next to it.
S="$(dirname "${BASH_SOURCE[0]}")"
SUBMITTED=0
command -v sbatch >/dev/null 2>&1 && HAVE_SBATCH=1 || HAVE_SBATCH=0
for g in "${GROUP_DIRS[@]}"; do
    if [ "$HAVE_SBATCH" = 1 ]; then
        # A refused submission used to pass unremarked, and the run ended saying it was done
        # while the group held archives and nothing had been queued to unpack them. A cluster
        # with no default partition refuses every submission that does not name one, which is
        # what SBATCH_PARTITION in your environment is for.
        if (cd "$PIPELINE_ROOT" && mkdir -p logs &&
            sbatch -J sra2fq -c 8 --mem=16G -o "logs/sra2fq_%j.out"                    --wrap "bash '$S/sra_to_fastq.sh' '$g'"); then
            SUBMITTED=$((SUBMITTED+1))
        else
            echo "[FAIL] could not submit the unpacking job for $g" >&2
            FAILED+=("$g (unpacking not submitted; run Scripts/sra_to_fastq.sh on it)")
        fi
    else
        bash "$S/sra_to_fastq.sh" "$g" || FAILED+=("$g (conversion)")
    fi
done

# Metadata last, not after each group: a GEO hiccup then sits at the end of the log where you
# will see it, instead of scrolled past in the middle of a run that kept going for hours.
for g in "${GROUP_DIRS[@]}"; do
    bash "$S/fetch_metadata.sh" "$g" || META_FAILED+=("$g")
done

echo
echo "[DONE] ${#GROUP_DIRS[@]} group(s), $TOTAL run(s)"
printf '       %s\n' "${GROUP_DIRS[@]}"

if [ ${#FAILED[@]} -gt 0 ] || [ ${#META_FAILED[@]} -gt 0 ]; then
    echo
    [ ${#FAILED[@]}      -eq 0 ] || printf '[FAIL] %s\n' "${FAILED[@]}" >&2
    [ ${#META_FAILED[@]} -eq 0 ] || printf '[FAIL] metadata: %s\n' "${META_FAILED[@]}" >&2
    echo "       Re-run the same command: finished runs are skipped, only these are retried." >&2
    exit 1
fi

[ "$SUBMITTED" -eq 0 ] || cat <<MSG

  $SUBMITTED unpacking job(s) submitted. Until they finish the groups hold archives rather
  than reads, and the pipeline has nothing to work on.

      squeue -u \$USER
      tail -f ${PIPELINE_ROOT}/logs/sra2fq_*.out
MSG

cat <<MSG

  Next steps, once per group:

  1. Look at <group>/metadata.tsv and decide which samples are which.
     The pipeline does not read it; you do, to tell the count-matrix columns apart.

  2. <group>/group.conf is written for you, with the species taken from runinfo.
     Nothing else needs filling in - the probe records the strandedness itself.

  3. sbatch Scripts/run_pipeline.sh ${GROUP_DIRS[0]}
MSG

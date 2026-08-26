#!/bin/bash
# PublicData_download.sh <group_dir> <SRR ...>
#   Downloads raw FASTQ for a list of SRR accessions into the group folder.
#   Runs sharing a SampleName go into a subfolder, which marks them as one sample to merge.
#   Sample metadata is collected afterwards by fetch_metadata.sh.
#
#   bash PublicData_download.sh rawData/ProjectA/GroupA SRR0000001 SRR0000002
#   bash PublicData_download.sh rawData/ProjectA/GroupA rawData/ProjectA/GroupA/accessions.csv
#
#   The list file has no format: every SRR/ERR/DRR accession found anywhere in it is used,
#   duplicates dropped. A plain list, a comma-separated line, or a whole SraRunTable.csv
#   pasted in all work, so there is nothing to reformat before running this.
set -euo pipefail

RUNINFO_URL="https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc="
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"

# gzip is single-threaded and dominates the wall clock on a multi-GB FASTQ. pigz writes the
# same format on N cores; plain gzip stays the fallback.
command -v pigz >/dev/null 2>&1 && ZIP="pigz -p $THREADS" || ZIP="gzip"

[ $# -ge 2 ] || { sed -n '2,7p' "$0"; exit 1; }

GROUP_DIR="$1"; shift
# A single file argument is read as a list of accessions. Scraping them out rather than
# parsing the file means a run table or a copied web page works without editing.
if [ $# -eq 1 ] && [ -f "$1" ]; then
    mapfile -t ACCS < <(grep -oE '[SED]RR[0-9]+' "$1" | awk '!seen[$0]++')
else
    ACCS=("$@")
fi
[ ${#ACCS[@]} -gt 0 ] || { echo "[ERROR] no accession given"; exit 1; }

mkdir -p "$GROUP_DIR"
# runinfo is machinery, not the user's metadata: it exists only to group runs by sample.
# The conditions (tissue, treatment, donor) are not in it - the user fetches those separately.
META="${GROUP_DIR}/.runinfo.csv"

# -- 1. runinfo -------------------------------------------------------------
# The endpoint takes a comma-separated list, so this is one request rather than one per run
curl -sf "${RUNINFO_URL}$(IFS=,; echo "${ACCS[*]}")" > "$META"
[ -s "$META" ] || { echo "[ERROR] runinfo lookup failed"; exit 1; }
echo "[INFO] runinfo: $(( $(wc -l < "$META") - 1 )) runs"
SRP=$(awk -F, 'NR==1 { for (i=1;i<=NF;i++) if ($i=="SRAStudy") c=i; next } c { print $c; exit }' "$META")

# ── 2. acc -> SampleName map, one pass. Samples with >1 run go into a subfolder.
declare -A SAMPLE_OF MULTI count
while IFS=$'	' read -r acc smp; do
    SAMPLE_OF[$acc]="$smp"
    n=$(( ${count[$smp]:-0} + 1 )); count[$smp]=$n
    if [ "$n" -gt 1 ]; then MULTI[$smp]=1; fi
done < <(awk -F, '
    NR==1 { for (i=1;i<=NF;i++) if ($i=="SampleName") c=i; next }
    c { print $1 "	" $c }' "$META")

# -- 3. download ------------------------------------------------------------
for acc in "${ACCS[@]}"; do
    sample="${SAMPLE_OF[$acc]:-$acc}"
    dest="${GROUP_DIR}"
    [ -n "${MULTI[$sample]:-}" ] && dest="${GROUP_DIR}/${sample}"   # multi-run -> merge folder
    mkdir -p "$dest"

    if [ -f "${dest}/.${acc}.done" ]; then
        echo "[SKIP] $acc"
        continue
    fi
    echo "[GET ] $acc -> $dest"
    prefetch --output-directory "$dest" "$acc" >/dev/null
    fasterq-dump --split-3 --threads "$THREADS" --outdir "$dest" "${dest}/${acc}" >/dev/null
    rm -rf "${dest:?}/${acc}"
    $ZIP -f "${dest}"/${acc}*.fastq
    touch "${dest}/.${acc}.done"
done

bash "$(dirname "${BASH_SOURCE[0]}")/fetch_metadata.sh" "$GROUP_DIR" ||     echo "[WARN] metadata fetch failed - rerun Scripts/fetch_metadata.sh $GROUP_DIR"

cat <<MSG

[DONE] ${#ACCS[@]} runs -> $GROUP_DIR

  Next steps:

  1. Look at ${GROUP_DIR}/metadata.tsv and decide which samples are which.
     The pipeline does not read it; you do, to tell the count-matrix columns apart.

  2. Write ${GROUP_DIR}/group.conf

         species      = Homo_sapiens     # must match a folder in reference_Genomes/
         strandedness =                  # leave empty; the probe fills it in

  3. bash Scripts/run_pipeline.sh ${GROUP_DIR}
MSG

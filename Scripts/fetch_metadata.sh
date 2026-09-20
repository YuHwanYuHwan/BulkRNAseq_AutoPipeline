#!/bin/bash
# fetch_metadata.sh <group_dir>
#   Collects the sample metadata for a downloaded group into metadata.tsv.
#
#   runinfo tells you a run's accession and its sample, and nothing about what the sample is.
#   The conditions live in GEO (or, for submissions that never went through GEO, in BioSample),
#   so they are fetched from there and joined on the sample accession.
#
#   Fields are recorded as the submitter wrote them, including the title AND the
#   characteristics. They disagree in real datasets - only you can tell which one is right.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"

GEO="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
EUTILS="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# NCBI answers three requests a second and 429s the fourth, and the BioSample path makes two
# requests per sample. Pausing costs seconds over a run; one refused request used to cost the
# whole table, because curl -sf reports it as a failure and the script stops on those.
eutils() {   # $1 = everything after the endpoint
    local try out
    for try in 1 2 3; do
        sleep 0.4
        out=$(curl -sf --max-time 60 "${EUTILS}/$1") && { printf '%s' "$out"; return 0; }
    done
    return 1
}

[ $# -eq 1 ] || { sed -n '2,3p' "$0"; exit 1; }
GROUP_DIR="$1"
RUNINFO="${GROUP_DIR}/.runinfo.csv"
OUT="${GROUP_DIR}/metadata.tsv"
[ -s "$RUNINFO" ] || { echo "[ERROR] no $RUNINFO - run PublicData_download.sh first" >&2; exit 1; }

# Run, SampleName, BioSample, BioProject - all of it already sits in runinfo
RUNS=$(tr -d '\r' < "$RUNINFO" | awk -F, 'NR==1 { for (i=1;i<=NF;i++) { if ($i=="Run") r=i; if ($i=="SampleName") s=i
                                              if ($i=="BioSample") b=i; if ($i=="BioProject") p=i }
                        next }
                r { print $r "\t" $s "\t" $b "\t" $p }')
FIRST_SAMPLE=$(head -1 <<< "$RUNS" | cut -f2)

RAW="${GROUP_DIR}/.geo_samples.txt"
if [[ "$FIRST_SAMPLE" == GSM* ]]; then
    # One request for the whole series beats one per sample
    # GEO serves this text with CRLF line endings. A carriage return riding along on the
    # series id makes the next URL malformed, and curl fails with a code that says nothing
    # about where it came from.
    SERIES=$(curl -sf --max-time 60 "${GEO}?acc=${FIRST_SAMPLE}&targ=self&form=text&view=brief" |
             tr -d '\r' | awk -F' = ' '/^!Sample_series_id/ { print $2; exit }')
    [ -n "$SERIES" ] || { echo "[ERROR] no series for $FIRST_SAMPLE" >&2; exit 1; }
    echo "[GEO ] $SERIES"
    curl -sf --max-time 120 "${GEO}?acc=${SERIES}&targ=gsm&form=text&view=brief" | tr -d '\r' > "$RAW"
else
    # Not a GEO submission: BioSample carries the same attributes, one fetch per sample.
    # Emitted in the GEO shape so the parser below does not need a second form.
    echo "[BIOS] $FIRST_SAMPLE is not a GSM - reading BioSample instead"
    : > "$RAW"
    missed=0
    while IFS=$'\t' read -r run smp bios proj; do
        [ -n "$bios" ] || continue
        uid=$(eutils "esearch.fcgi?db=biosample&term=${bios}" |
              grep -oE '<Id>[0-9]+' | head -1 | cut -d'>' -f2) || true
        [ -n "$uid" ] || { echo "[WARN] no BioSample record for $bios" >&2; missed=$((missed+1)); continue; }
        eutils "efetch.fcgi?db=biosample&id=${uid}&rettype=full&retmode=xml" |
        awk -v s="$bios" '
            BEGIN { print "^SAMPLE = " s }
            /<Title>/    { t=$0; gsub(/.*<Title>|<\/Title>.*/,"",t); print "!Sample_title = " t }
            /attribute_name=/ {
                k=$0; gsub(/.*attribute_name="/,"",k); gsub(/".*/,"",k)
                v=$0; gsub(/.*>/,"",v); gsub(/<.*/,"",v)
                if (v != "") print "!Sample_characteristics_ch1 = " k ": " v }' >> "$RAW" ||
            { echo "[WARN] could not read the record for $bios" >&2; missed=$((missed+1)); }
    done <<< "$RUNS"
    [ "$missed" -eq 0 ] || echo "[WARN] $missed sample(s) have no attributes below; re-running fills them in" >&2
fi
[ -s "$RAW" ] || { echo "[ERROR] no metadata retrieved" >&2; exit 1; }

# Wide table: one row per run, one column per characteristic key found anywhere in the set.
awk -F'\t' -v raw="$RAW" -v series="${SERIES:--}" '
    BEGIN {
        while ((getline line < raw) > 0) {
            if (line ~ /^\^SAMPLE/)                 { split(line,a," = "); s=a[2]; continue }
            if (line ~ /^!Sample_title/)            { split(line,a," = "); title[s]=a[2]; continue }
            if (line ~ /^!Sample_source_name_ch1/)  { split(line,a," = "); src[s]=a[2]; continue }
            if (line ~ /^!Sample_characteristics/) {
                split(line,a," = "); kv=a[2]
                i=index(kv,":"); if (!i) continue
                k=substr(kv,1,i-1); v=substr(kv,i+2)
                gsub(/^ +| +$/,"",k); gsub(/^ +| +$/,"",v)
                if (!(k in seen)) { seen[k]=++nk; key[nk]=k }
                val[s,k]=v
            }
        }
        printf "Run\tSample\tSeries\tBioProject\tBioSample\tTitle\tSourceName"
        for (j=1;j<=nk;j++) printf "\t%s", key[j]
        print ""
    }
    {
        s=$2
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s", $1, s, series, ($4?$4:"-"), ($3?$3:"-"),
               (title[s]?title[s]:"-"), (src[s]?src[s]:"-")
        for (j=1;j<=nk;j++) printf "\t%s", ((s,key[j]) in val ? val[s,key[j]] : "-")
        print ""
    }' <<< "$RUNS" > "$OUT"

rm -f "$RAW"
echo "[DONE] $OUT   ($(( $(wc -l < "$OUT") - 1 )) runs, $(head -1 "$OUT" | awk -F'\t' '{print NF-7}') attributes)"
echo "       Title and the characteristics are both recorded on purpose. Where they disagree,"
echo "       the submitter made a mistake in one of them and the pipeline cannot tell which."

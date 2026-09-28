#!/usr/bin/env bash
# Download child BMI GWAS sumstats (MoBa 2022 x12 + EGG GCST90002409) to E:/CM/GLP1
# Resumable: safe to re-run; curl -C - continues partial files.
set -u
DEST=/e/CM/GLP1
mkdir -p "$DEST"

MOBA_BASE="https://www.fhi.no/contentassets/380a34f1fc854629a42f7ee1ab2778f3/artikkel-2022"
TPS="birth 6weeks 3months 6months 8months 1year 1.5years 2years 3years 5years 7years 8years"

declare -a URLS FILES
for tp in $TPS; do
  URLS+=("$MOBA_BASE/childhood_bmi_nature_metab_2022_${tp}.gz")
  FILES+=("moba2022_bmi_${tp}.gz")
done
URLS+=("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002409/GCST90002409_buildGRCh37.tsv.gz")
FILES+=("EGG2020_GCST90002409_buildGRCh37.tsv.gz")
URLS+=("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002409/md5sum.txt")
FILES+=("EGG2020_md5sum.txt")
URLS+=("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002409/GCST90002409_buildGRCh37.tsv.gz-meta.yaml")
FILES+=("EGG2020_meta.yaml")

PIDS=()
for i in "${!URLS[@]}"; do
  f="$DEST/${FILES[$i]}"
  curl -fSL --retry 3 --retry-delay 5 --connect-timeout 30 -C - -o "$f" "${URLS[$i]}" \
    >"$f.log" 2>&1 &
  PIDS+=($!)
done

# Wait up to ~260s, then report status; caller re-runs script to continue.
END=$((SECONDS + 260))
while [ $SECONDS -lt $END ]; do
  alive=0
  for p in "${PIDS[@]}"; do kill -0 "$p" 2>/dev/null && alive=$((alive+1)); done
  [ "$alive" -eq 0 ] && break
  sleep 10
done

echo "--- status ---"
for i in "${!FILES[@]}"; do
  f="$DEST/${FILES[$i]}"
  sz=$(stat -c %s "$f" 2>/dev/null || echo 0)
  echo "${FILES[$i]}  ${sz} bytes"
done
# Kill any stragglers (partials persist for resume)
for p in "${PIDS[@]}"; do kill "$p" 2>/dev/null; done
wait 2>/dev/null
exit 0

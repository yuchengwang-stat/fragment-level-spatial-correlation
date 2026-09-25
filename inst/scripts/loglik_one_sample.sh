#!/usr/bin/env bash
# Step 2 for one sample, end to end: keep only the .pat lines that can touch a
# region, then score them.
#
#   loglik_one_sample.sh <regions.rds> <sample.pat.gz> <sample.beta> <out.rds> [threads]
#
# The prefiltered .pat is written next to the output as <out>.pat and removed
# by name when scoring succeeds (set KEEP_PAT=1 to keep it).
set -euo pipefail
if [ $# -lt 4 ]; then
  echo "usage: loglik_one_sample.sh <regions.rds> <sample.pat.gz> <sample.beta> <out.rds> [threads]" >&2
  exit 2
fi
REG=$1; PAT=$2; BETA=$3; OUT=$4; NC=${5:-1}
FC=$(cd "$(dirname "$0")" && pwd)
if [ -s "$OUT" ]; then echo "already done: $OUT"; exit 0; fi
mkdir -p "$(dirname "$OUT")"
IDX="${OUT%.rds}.regions.txt"
SMALL="${OUT%.rds}.pat"

Rscript "$FC/export_region_index.R" "$REG" "$IDX"
t0=$SECONDS
case "$PAT" in
  *.gz) gzip -dc "$PAT" ;;
  *)    cat "$PAT" ;;
esac | awk -v RF="$IDX" -f "$FC/prefilter_pat.awk" > "$SMALL"
echo "  prefilter $((SECONDS - t0))s, $(wc -l < "$SMALL") read patterns kept"

Rscript "$FC/02_loglik.R" --regions "$REG" --pat "$SMALL" --beta "$BETA" \
  --threads "$NC" --out "$OUT"

if [ "${KEEP_PAT:-0}" != "1" ]; then rm -f "$SMALL" "$IDX"; fi

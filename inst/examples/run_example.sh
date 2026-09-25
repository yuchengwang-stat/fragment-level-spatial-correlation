#!/usr/bin/env bash
# The whole pipeline on a small simulated dataset, in about a minute.
#
#   bash run_example.sh [workdir]        (default: ./fragcorr_example)
#
# Two cell types are simulated with known kernels, A: phi = 0.8, s = 100 and
# B: phi = 0.3, s = 300, three samples each.  The last step prints what the
# pipeline recovered next to the truth.
set -euo pipefail
W=${1:-fragcorr_example}
FC=$(Rscript -e 'cat(system.file("scripts", package = "fragcorr"))')
[ -n "$FC" ] || { echo "fragcorr is not installed" >&2; exit 1; }
mkdir -p "$W"

echo "== simulate"
Rscript -e "invisible(fragcorr::simulate_dataset('$W/data'))"

echo "== step 1: regions"
Rscript "$FC/01_regions.R" --cpg "$W/data/cpg_list.csv" --col hg38 \
  --ref "$W/data/cpg_index.rds" --out "$W/regions.rds"

echo "== step 2: log-likelihoods, one sample at a time"
tail -n +2 "$W/data/samples.tsv" | while IFS=$'\t' read -r sample celltype pat beta; do
  bash "$FC/loglik_one_sample.sh" "$W/regions.rds" "$pat" "$beta" "$W/lik/$celltype/$sample.rds" 2
done

echo "== step 3: prior and posterior, one cell type at a time"
for d in "$W"/lik/*/; do
  ct=$(basename "$d")
  Rscript "$FC/03_posterior.R" --lik "$d" --out "$W/post/$ct.rds" --label "$ct"
done

echo "== uniformity of each region's methylation, per cell type"
Rscript "$FC/09_heterogeneity.R" --regions "$W/regions.rds" --sheet "$W/data/samples.tsv" \
  --ref "$W/data/cpg_index.rds" --cpg "$W/data/cpg_list.csv" --col hg38 --out "$W/out/het.tsv"

echo "== step 4: result tables"
Rscript "$FC/04_region_table.R" --regions "$W/regions.rds" --cpg "$W/data/cpg_list.csv" \
  --col hg38 --post "$W/post" --het "$W/out/het.tsv" --out "$W/out/region_ct"

echo "== recovered against the truth"
Rscript -e "
x <- data.table::fread('$W/out/region_ct_long.csv')
print(x[, .(regions = .N, phi_mean = round(mean(phi_mean), 3), s_mean = round(mean(s_mean), 1)),
        by = celltype])
cat('truth: A phi 0.8 s 100 | B phi 0.3 s 300\n')"

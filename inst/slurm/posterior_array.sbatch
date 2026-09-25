#!/bin/bash -l
# Step 3 as a SLURM array: one task per cell type.
#
#   ls $OUT/lik > $OUT/celltypes.txt
#   sbatch --array=1-$(wc -l < $OUT/celltypes.txt) posterior_array.sbatch
#
# Output: $OUT/post/<celltype>.rds and the fitted prior in $OUT/prior/<celltype>.rds.
#SBATCH -J fc_post
#SBATCH -p <partition>
#SBATCH -A <account>
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH -t 01:00:00
#SBATCH -o logs/post_%A_%a.out
set -uo pipefail
OUT=/path/to/output
export R_LIBS_USER=${R_LIBS_USER:-$HOME/R/library}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

FC=$(Rscript -e 'cat(system.file("scripts", package = "fragcorr"))')
CT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$OUT/celltypes.txt")
[ -n "$CT" ] || { echo "no cell type at line $SLURM_ARRAY_TASK_ID"; exit 1; }
mkdir -p "$OUT/post" "$OUT/prior"
Rscript "$FC/03_posterior.R" --lik "$OUT/lik/$CT" --out "$OUT/post/$CT.rds" \
  --prior "$OUT/prior/$CT.rds" --label "$CT"

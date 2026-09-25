#!/bin/bash -l
# Step 2 as a SLURM array: one task per sample.
#
#   sbatch --array=1-<number of samples> loglik_array.sbatch
#
# Set the four paths below.  SHEET is a tab-separated file with a header line and
# columns sample, celltype, pat, beta (absolute paths to the .pat.gz and .beta).
# Output: $OUT/lik/<celltype>/<sample>.rds.  A task whose output exists is skipped,
# so the array can be resubmitted after a partial failure.
#SBATCH -J fc_lik
#SBATCH -p <partition>
#SBATCH -A <account>
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH -t 03:00:00
#SBATCH -o logs/lik_%A_%a.out
set -uo pipefail
SHEET=/path/to/samples.tsv
REGIONS=/path/to/regions.rds
OUT=/path/to/output
export R_LIBS_USER=${R_LIBS_USER:-$HOME/R/library}
# SLURM does not set these; without them BLAS / OpenMP see every core of the node
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
NC=${SLURM_NTASKS_PER_NODE:-1}

FC=$(Rscript -e 'cat(system.file("scripts", package = "fragcorr"))')
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SHEET")
[ -n "$LINE" ] || { echo "no sample at line $((SLURM_ARRAY_TASK_ID + 1)) of $SHEET"; exit 1; }
IFS=$'\t' read -r SAMPLE CT PAT BETA <<< "$LINE"
echo "=== $SLURM_ARRAY_TASK_ID  $CT / $SAMPLE ==="
bash "$FC/loglik_one_sample.sh" "$REGIONS" "$PAT" "$BETA" "$OUT/lik/$CT/$SAMPLE.rds" "$NC"

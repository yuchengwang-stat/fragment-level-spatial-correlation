#!/usr/bin/env Rscript
## Step 4 -- the result tables: per region, its coordinates, the listed CpGs it
## holds, and per cell type the winning (phi, s), the posterior means and the
## region's methylation.
##
##   Rscript 04_region_table.R --regions regions.rds --cpg cpgs.csv --post <dir> --out prefix
##     [--col hg38] [--category ...] [--het het.tsv]
##
## --post is the directory of step-3 files, one per cell type; the file name is
## the cell-type label.  --cpg, --col and --category must match step 1.
## Writes <prefix>_long.csv (one row per region x cell type) and
## <prefix>_wide.csv (one row per region).  --het adds the uniformity columns of
## 09_heterogeneity.R (meth_sd, meth_jump, r2_split, target_lowhigh, layout) to
## the long table.
suppressMessages(library(fragcorr))
a <- fc_args(list(regions = "", cpg = "", post = "", out = "", col = "hg38", category = "",
                  het = ""))
fc_require(a, c("regions", "cpg", "post", "out"),
           "04_region_table.R --regions <rds> --cpg <csv> --post <dir> --out <prefix>")
cats <- if (nzchar(a$category)) strsplit(a$category, ",")[[1]] else NULL
tabs <- region_tables(readRDS(a$regions), read_cpg_list(a$cpg, col = a$col, categories = cats),
                      read_posteriors(a$post),
                      het = if (nzchar(a$het)) data.table::fread(a$het) else NULL)
dir.create(dirname(a$out), showWarnings = FALSE, recursive = TRUE)
fl <- paste0(a$out, "_long.csv"); fw <- paste0(a$out, "_wide.csv")
data.table::fwrite(tabs$long, fl); data.table::fwrite(tabs$wide, fw)
cat("wrote ", fl, " (", nrow(tabs$long), " rows) and ", fw, " (", nrow(tabs$wide), " x ",
    ncol(tabs$wide), ")\n", sep = "")

#!/usr/bin/env Rscript
## How uniform is each region's methylation, per cell type?  Reads the .beta files
## of every sample in the sheet, pools counts per cell type for every CpG inside
## the regions, and writes one row per region x cell type.  See
## heterogeneity_stats() for the columns; pass the output to 04_region_table.R
## with --het.
##
##   Rscript 09_heterogeneity.R --regions regions.rds --sheet samples.tsv
##     --ref cpg_index.rds --cpg cpgs.csv --out het.tsv
##     [--col hg38] [--category ...] [--min-n 20] [--low 0.2] [--high 0.8]
##     [--sd-uniform 0.10] [--r2-cut 0.75]
##
## --sheet has columns sample, celltype, beta (and pat, which is not used here).
## --cpg / --col / --category mark which CpGs are listed (targets); they must match
## step 1.
suppressMessages(library(fragcorr))
a <- fc_args(list(regions = "", sheet = "", ref = "", cpg = "", out = "", col = "hg38",
                  category = "", min_n = "20", low = "0.2", high = "0.8",
                  sd_uniform = "0.10", r2_cut = "0.75"))
fc_require(a, c("regions", "sheet", "ref", "cpg", "out"),
           "09_heterogeneity.R --regions <rds> --sheet <tsv> --ref <cpg index> --cpg <csv> --out <tsv>")
reg <- readRDS(a$regions)
sheet <- data.table::fread(a$sheet)
cats <- if (nzchar(a$category)) strsplit(a$category, ",")[[1]] else NULL
cat(nrow(reg), " regions, ", nrow(sheet), " samples in ", length(unique(sheet$celltype)),
    " groups\n", sep = "")
cnt <- pooled_region_counts(reg, sheet)
tgt <- target_rows(reg, read_cpg_list(a$cpg, col = a$col, categories = cats),
                   load_cpg_index(a$ref))
h <- heterogeneity_stats(cnt$M, cnt$N, cnt$region_id, tgt, min_n = as.integer(a$min_n),
                         low = as.numeric(a$low), high = as.numeric(a$high),
                         sd_uniform = as.numeric(a$sd_uniform), r2_cut = as.numeric(a$r2_cut))
dir.create(dirname(a$out), showWarnings = FALSE, recursive = TRUE)
data.table::fwrite(h, a$out, sep = "\t")
cat("wrote ", a$out, " (", nrow(h), " rows)\n\n", sep = "")
print(h[, .N, by = layout][order(-N)], row.names = FALSE)

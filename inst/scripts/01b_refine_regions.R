#!/usr/bin/env Rscript
## Step 1b -- split regions whose target CpGs sit at opposite methylation
## extremes, so that each region's targets share one methylation state.
##
## Methylation is pooled per cell type from the .beta files.  A run of targets is
## cut when, in at least --min-celltypes cell types, one of its targets is at
## <= --low and another at >= --high (pooled depth >= --min-n).  The split is the
## same for every cell type.  Pieces left with one target are kept and marked in
## the result tables (single_target); every piece records its parent region.
##
##   Rscript 01b_refine_regions.R --regions regions.rds --sheet samples.tsv
##     --ref cpg_index.rds --cpg cpgs.csv --out regions_refined.rds
##     [--col hg38] [--category ...] [--min-celltypes 5] [--low 0.2] [--high 0.8]
##     [--min-n 20] [--pad 250]
##
## --pad should equal step 1's.  Use the refined regions in every later step,
## with the same --cpg / --col / --category.
suppressMessages(library(fragcorr))
a <- fc_args(list(regions = "", sheet = "", ref = "", cpg = "", out = "", col = "hg38",
                  category = "", min_celltypes = "5", low = "0.2", high = "0.8",
                  min_n = "20", pad = "250"))
fc_require(a, c("regions", "sheet", "ref", "cpg", "out"),
           "01b_refine_regions.R --regions <rds> --sheet <tsv> --ref <cpg index> --cpg <csv> --out <rds>")
reg <- readRDS(a$regions)
sheet <- data.table::fread(a$sheet)
cats <- if (nzchar(a$category)) strsplit(a$category, ",")[[1]] else NULL
cpg <- read_cpg_list(a$cpg, col = a$col, categories = cats)
ref <- load_cpg_index(a$ref)
cat(nrow(reg), " regions; pooling ", nrow(sheet), " samples in ",
    length(unique(sheet$celltype)), " cell types\n", sep = "")

cnt <- pooled_region_counts(reg, sheet)
P <- ifelse(cnt$N >= as.integer(a$min_n), cnt$M / cnt$N, NA_real_)
tg <- match_targets(reg, cpg, ref)
cat("  ", nrow(tg), " targets in regions, ", sum(is.na(tg$row)),
    " without an indexed CpG (no methylation; kept as targets)\n", sep = "")

new <- refine_regions(reg, tg, P, min_celltypes = as.integer(a$min_celltypes),
                      low = as.numeric(a$low), high = as.numeric(a$high), pad = as.integer(a$pad))
validate_regions(new, max_width = max(reg$width), min_cpgs = 1L)
new <- index_regions(new, ref)
validate_indexed(new)
if (sum(new$n_cpg) != nrow(tg)) stop("targets were lost or duplicated in the split")

n_split <- sum(table(new$parent_region) > 1L)
cat(sprintf("  %d of %d regions split; %d regions now; %d hold a single target\n",
            n_split, nrow(reg), nrow(new), sum(new$n_cpg == 1L)))
saveRDS(new, a$out)
cat("wrote ", a$out, "\n", sep = "")

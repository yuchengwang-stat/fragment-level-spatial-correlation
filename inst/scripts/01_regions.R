#!/usr/bin/env Rscript
## Step 1 -- regions from a CpG list, indexed against the genome-wide CpG table.
##
## Each listed CpG is padded by --pad bp on both sides and overlapping pads are
## merged.  A merged cluster wider than --max-width is cut, between CpGs, into
## near-equal pieces.  Regions holding fewer than --min-cpgs listed CpGs are
## dropped.  Each region is then given the CpG index range it spans in the
## reference table (row k of that table is row k of every .beta) and the bp gaps
## between consecutive CpGs.  Run once; every later step reads this file.
##
##   Rscript 01_regions.R --cpg cpgs.csv --ref cpg_index.rds --out regions.rds
##     [--col hg38] [--category HighValue,EPIC] [--pad 250] [--max-width 1000]
##     [--min-cpgs 2]
suppressMessages(library(fragcorr))
a <- fc_args(list(cpg = "", ref = "", out = "", col = "hg38", category = "",
                  pad = "250", max_width = "1000", min_cpgs = "2"))
fc_require(a, c("cpg", "ref", "out"), "01_regions.R --cpg <csv> --ref <cpg index> --out <rds>")

cats <- if (nzchar(a$category)) strsplit(a$category, ",")[[1]] else NULL
cpg <- read_cpg_list(a$cpg, col = a$col, categories = cats)
cat(nrow(cpg), " CpGs on ", length(unique(cpg$chr)), " chromosomes\n", sep = "")

pad <- as.integer(a$pad); mw <- as.integer(a$max_width); mc <- as.integer(a$min_cpgs)
cat("building regions: pad +/-", pad, " bp, max width ", mw, " bp, at least ", mc, " listed CpGs\n", sep = "")
reg <- build_regions(cpg, pad = pad, max_width = mw, min_cpgs = mc)
validate_regions(reg, max_width = mw, min_cpgs = mc)
cat("  ", nrow(reg), " regions, ", round(sum(reg$width) / 1e6, 3), " Mb\n", sep = "")

cat("indexing against ", a$ref, "\n", sep = "")
reg <- index_regions(reg, load_cpg_index(a$ref))
validate_indexed(reg)
cat("  ", nrow(reg), " regions hold ", sum(reg$n_cpg_in_region), " CpGs in total\n", sep = "")
saveRDS(reg, a$out)
cat("wrote ", a$out, "\n", sep = "")

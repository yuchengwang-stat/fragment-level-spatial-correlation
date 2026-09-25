#!/usr/bin/env Rscript
## Write the regions' CpG index ranges as a sorted two-column file, which is what
## prefilter_pat.awk reads.
##
##   Rscript export_region_index.R <regions.rds> <out.txt> [chr1,chr2,...]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: export_region_index.R <regions.rds> <out.txt> [chroms]")
reg <- readRDS(args[1])
if (length(args) >= 3 && nzchar(args[3])) reg <- reg[reg$chr %in% strsplit(args[3], ",")[[1]], ]
reg <- reg[order(reg$startCpG), ]
utils::write.table(reg[, c("startCpG", "endCpG")], args[2],
                   sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
cat(nrow(reg), " region index ranges -> ", args[2], "\n", sep = "")

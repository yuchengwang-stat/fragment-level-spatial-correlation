#!/usr/bin/env Rscript
## Correlation against distance, one curve per methylation bin, with a bootstrap
## band across cell types.  See decay_curves() for the definition.
##
##   Rscript 07_decay_curves.R --post <dir of step-3 files> --out curves.rds
##     [--dmax 400] [--dstep 4] [--boot 4000] [--seed 1] [--min-regions 20]
suppressMessages(library(fragcorr))
a <- fc_args(list(post = "", out = "", dmax = "400", dstep = "4", boot = "4000",
                  seed = "1", min_regions = "20"))
fc_require(a, c("post", "out"), "07_decay_curves.R --post <dir> --out <rds>")
x <- decay_curves(read_posteriors(a$post), dmax = as.numeric(a$dmax), dstep = as.numeric(a$dstep),
                  boot = as.integer(a$boot), seed = as.integer(a$seed),
                  min_regions = as.integer(a$min_regions))
dir.create(dirname(a$out), showWarnings = FALSE, recursive = TRUE)
saveRDS(x, a$out)
data.table::fwrite(x$band, sub("[.]rds$", ".tsv", a$out), sep = "\t")
cat("wrote ", a$out, "\n\n", sep = "")
cat("======== per bin ========\n"); print(x$summary, row.names = FALSE)
show <- x$band[d %in% c(0, 50, 100, 200, 400)]
if (nrow(show)) {
  show[, txt := sprintf("%.3f [%.3f, %.3f]", mean, lo, hi)]
  cat("\n======== mean curve, 95% interval across cell types ========\n")
  print(data.table::dcast(show, bin ~ d, value.var = "txt"), row.names = FALSE)
}

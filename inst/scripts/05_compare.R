#!/usr/bin/env Rscript
## Diagnostic -- where is the information, and do cell types agree?
##
##   * A region that is almost fully (un)methylated has little variation to
##     correlate, so its posterior should sit closer to the prior: a larger
##     effective number of patterns.
##   * The same region is fitted separately in each cell type, so a region whose
##     correlation structure is cell-type specific will have different winners.
##
##   Rscript 05_compare.R --post <dir of step-3 files> [--out compare.rds] [--min-ct 5]
suppressMessages(library(fragcorr))
a <- fc_args(list(post = "", out = "", min_ct = "5"))
fc_require(a, "post", "05_compare.R --post <dir>")
post <- read_posteriors(a$post)
cat(length(post), " cell types\n\n", sep = "")
P <- data.table::rbindlist(lapply(names(post), function(ct) {
  p <- data.table::as.data.table(post[[ct]])
  p[, list(celltype = ct, chr, interval_id, bin = as.character(bin), best = best_pattern,
           best_s, best_phi, best_post, eff = eff_patterns)]
}))

cat("======== information by methylation bin ========\n")
cat("(effective patterns: 1 = the data pick one pattern, ", nrow(fc_patterns()),
    " = the data say nothing)\n\n", sep = "")
bt <- P[!is.na(eff), list(regions = .N, eff_median = round(stats::median(eff), 2),
                          winner_post = round(stats::median(best_post), 3)), by = bin][
  order(match(bin, fc_bins()$label))]
print(bt, row.names = FALSE)

cat("\n======== winning decay length, by bin ========\n")
print(data.table::dcast(P[!is.na(best_s)], bin ~ best_s, fun.aggregate = length,
                        value.var = "best_s"), row.names = FALSE)

cat("\n======== do cell types agree on a region? ========\n")
ag <- P[!is.na(best), list(n_ct = .N, n_distinct = data.table::uniqueN(best),
                           top_frac = max(table(best)) / .N), by = list(chr, interval_id)]
ag <- ag[n_ct >= as.integer(a$min_ct)]
if (nrow(ag)) {
  cat(sprintf("  %d regions seen in at least %s cell types\n", nrow(ag), a$min_ct))
  cat(sprintf("  distinct winning patterns per region: median %g, max %d\n",
              stats::median(ag$n_distinct), max(ag$n_distinct)))
  cat(sprintf("  regions where every cell type agrees: %d (%.1f%%)\n",
              sum(ag$n_distinct == 1), 100 * mean(ag$n_distinct == 1)))
} else cat("  no region is seen in enough cell types\n")
if (nzchar(a$out)) { saveRDS(list(long = P, agree = ag, by_bin = bt), a$out)
                     cat("\nwrote ", a$out, "\n", sep = "") }

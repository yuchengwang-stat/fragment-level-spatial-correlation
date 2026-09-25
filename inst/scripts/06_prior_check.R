#!/usr/bin/env Rscript
## Diagnostic -- is disagreement between cell types real, or the prior?
##
## Each cell type fits its own prior, in five methylation bins, possibly from
## very few samples.  With a modest winning posterior the argmax is easy to move,
## so a noisy prior alone could make cell types disagree.  Two controls:
##   shared  refit every cell type against one prior pooled over all of them;
##           if the disagreement collapses, it was the prior
##   deep    keep only cell types with at least --min-samples samples
##
##   Rscript 06_prior_check.R --lik-dir <dir with one subdirectory of step-2 files
##     per cell type> --post-dir <dir of step-3 files> [--min-samples 4] [--out x.rds]
suppressMessages(library(fragcorr))
a <- fc_args(list(lik_dir = "", post_dir = "", min_samples = "4", out = "", min_ct = "5"))
fc_require(a, c("lik_dir", "post_dir"), "06_prior_check.R --lik-dir <dir> --post-dir <dir>")
pat <- fc_patterns(); pn <- pat$name
cts <- list.dirs(a$lik_dir, recursive = FALSE)
nsamp <- vapply(cts, function(d) length(list.files(d, "[.]rds$")), 0L)
cat(length(cts), " cell types; ", sum(nsamp >= as.integer(a$min_samples)),
    " have >= ", a$min_samples, " samples\n\n", sep = "")

lik_of <- lapply(cts, function(d) {
  parts <- lapply(list.files(d, "[.]rds$", full.names = TRUE), readRDS)
  L <- combine_loglik(lapply(parts, `[[`, "loglik"), pn)
  m <- data.table::rbindlist(lapply(parts, function(p) p$loglik[, list(chr, interval_id, M_sum, N_sum)]))[
    , list(M_sum = sum(M_sum), N_sum = sum(N_sum)), by = list(chr, interval_id)]
  L <- merge(L, m, by = c("chr", "interval_id"))
  L[, bin := assign_bin(M_sum / N_sum)][]
})
names(lik_of) <- basename(cts)

cat("fitting the shared prior\n")
pooled <- data.table::rbindlist(lik_of, use.names = TRUE, fill = TRUE)
shared <- fit_priors(pooled, pn, max_fit = 200000L)
winners <- function(L, pri) {
  s <- summarise_posterior(posterior_from_prior(L, pri$pi, pn), pat)
  data.table::data.table(chr = L$chr, interval_id = L$interval_id, best = s$best_pattern,
                         phi = s$best_phi, post = s$best_post)
}
agree <- function(W, label) {
  ag <- W[!is.na(best), list(n_ct = .N, n_distinct = data.table::uniqueN(best),
                             top_frac = max(table(best)) / .N), by = list(chr, interval_id)]
  ag <- ag[n_ct >= as.integer(a$min_ct)]
  if (!nrow(ag)) { cat(sprintf("  %-30s no region seen in enough cell types\n", label)); return(invisible(ag)) }
  cat(sprintf("  %-30s regions %6d | distinct winners median %g | all agree %5.1f%% | majority<half %5.1f%%\n",
              label, nrow(ag), stats::median(ag$n_distinct),
              100 * mean(ag$n_distinct == 1), 100 * mean(ag$top_frac < 0.5)))
  invisible(ag)
}
cat("\n======== disagreement between cell types ========\n")
W_own <- data.table::rbindlist(lapply(names(lik_of), function(nm) {
  p <- data.table::as.data.table(readRDS(file.path(a$post_dir, paste0(nm, ".rds")))$posterior)
  p[, list(chr, interval_id, best = best_pattern, phi = best_phi, post = best_post, ct = nm)]
}))
agree(W_own, "own prior, all cell types")
W_sh <- data.table::rbindlist(lapply(names(lik_of), function(nm) {
  w <- winners(lik_of[[nm]], shared); w[, ct := nm][] }))
agree(W_sh, "shared prior, all cell types")
deep <- basename(names(nsamp)[nsamp >= as.integer(a$min_samples)])
agree(W_own[ct %in% deep], sprintf("own prior, >=%s samples", a$min_samples))
agree(W_sh[ct %in% deep],  sprintf("shared prior, >=%s samples", a$min_samples))
if (nzchar(a$out)) saveRDS(list(shared_prior = shared, own = W_own, shared = W_sh), a$out)

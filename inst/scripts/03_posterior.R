#!/usr/bin/env Rscript
## Step 3 -- a prior per methylation bin, then a posterior over the correlation
## patterns for every region.  Run once per cell type (or any group of samples).
##
## Log-likelihoods add across the group's samples; methylation is pooled by
## summing counts.  Within each methylation bin, the prior is the mixture weight
## of each pattern that maximises the marginal likelihood of the bin's regions
## (mixsqp).  The posterior is prior x likelihood, normalised.
##
##   Rscript 03_posterior.R --lik <dir of step-2 files, or comma-separated files>
##     --out celltype.rds [--prior prior.rds] [--label name] [--min-pi 0] [--max-fit 0]
##
## --prior: if the file exists it is used instead of fitting (so several groups
## can share one prior); if it does not, the fitted prior is written there.
## --max-fit: fit each bin's prior on at most this many regions (0 = all).
suppressMessages(library(fragcorr))
a <- fc_args(list(lik = "", out = "", prior = "", min_pi = "0", label = "", max_fit = "0"))
fc_require(a, c("lik", "out"), "03_posterior.R --lik <dir|files> --out <rds>")

files <- if (dir.exists(a$lik)) list.files(a$lik, "[.]rds$", full.names = TRUE) else
  strsplit(a$lik, ",")[[1]]
if (!length(files)) stop("no likelihood files found at ", a$lik)
cat(length(files), " sample(s)\n", sep = "")

parts <- lapply(files, readRDS)
pat_grid <- parts[[1]]$patterns
pn <- pat_grid$name

lik <- combine_loglik(lapply(parts, `[[`, "loglik"), pn)
meth <- data.table::rbindlist(lapply(parts, function(p)
  p$loglik[, list(chr, interval_id, M_sum, N_sum)]))[
  , list(M_sum = sum(M_sum, na.rm = TRUE), N_sum = sum(N_sum, na.rm = TRUE)),
  by = list(chr, interval_id)]
lik <- merge(lik, meth, by = c("chr", "interval_id"), all.x = TRUE)
lik[, meth := M_sum / N_sum]
lik[, bin := assign_bin(meth)]
cat(nrow(lik), " regions\n", sep = "")

if (nzchar(a$prior) && file.exists(a$prior)) {
  pri <- readRDS(a$prior)
  cat("prior read from ", a$prior, "\n", sep = "")
} else {
  cat("fitting the prior in each methylation bin\n")
  pri <- fit_priors(lik, pn, min_pi = as.numeric(a$min_pi), max_fit = as.integer(a$max_fit))
  if (nzchar(a$prior)) { saveRDS(pri, a$prior); cat("  prior written to ", a$prior, "\n", sep = "") }
}
cat("regions per bin used for the prior:\n"); print(pri$n)

post <- posterior_from_prior(lik, pri$pi, pn)
summ <- summarise_posterior(post, pat_grid)
out <- cbind(lik[, list(chr, interval_id, meth, bin)], summ, data.table::as.data.table(post))
dir.create(dirname(a$out), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(posterior = out, prior = pri, patterns = pat_grid,
             samples = basename(files), label = a$label), a$out)
cat("wrote ", a$out, "\n\n", sep = "")

cat("winning patterns (top 10):\n")
print(utils::head(sort(table(summ$best_pattern), decreasing = TRUE), 10))
cat(sprintf("\nposterior of the winner: median %.3f\n", stats::median(summ$best_post, na.rm = TRUE)))
cat(sprintf("effective number of patterns: median %.2f of %d\n",
            stats::median(summ$eff_patterns, na.rm = TRUE), length(pn)))

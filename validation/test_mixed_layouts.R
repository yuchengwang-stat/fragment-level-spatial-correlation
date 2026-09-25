#!/usr/bin/env Rscript
## Regions whose CpGs do not share one methylation rate.  Same set-up as
## test_bias_sources.R (truth phi 0.6, s 100; 8 CpGs; 60 fragments; marginals
## counted from the reads), but the per-site rate follows a layout.
suppressMessages({ library(fragcorr); library(data.table); library(parallel) })
PHI <- 0.6; S <- 100; NREG <- 300L; NFRAG <- 60L; NCPG <- 8L
pat <- fc_patterns()
layouts <- list(
  "uniform 0.475"            = rep(0.475, 8),
  "uniform 0.05"             = rep(0.05, 8),
  "uniform 0.90"             = rep(0.90, 8),
  "split 0.05 | 0.90"        = c(rep(0.05, 4), rep(0.90, 4)),
  "alternating 0.05 / 0.90"  = rep(c(0.05, 0.90), 4))
sim_region <- function(p) {
  gaps <- sample(20:120, NCPG - 1L, TRUE); pos <- c(0, cumsum(gaps))
  Sig <- PHI * exp(-abs(outer(pos, pos, "-")) / S); diag(Sig) <- 1
  C <- sweep(matrix(rnorm(NFRAG * NCPG), NFRAG, NCPG) %*% chol(Sig), 2, qnorm(p), `<=`)
  s <- apply(C, 1, function(v) paste(ifelse(v, "C", "T"), collapse = ""))
  tb <- table(s); n <- length(tb)
  L <- loglik_all_patterns(names(tb), rep(paste(gaps, collapse = ","), n),
                           rep(paste(colSums(C), collapse = ","), n),
                           rep(paste(rep(NFRAG, NCPG), collapse = ","), n),
                           pat$phi, pat$s, print_every = 0L)
  colSums(L * as.vector(tb), na.rm = TRUE)
}
res <- rbindlist(mclapply(names(layouts), function(nm) {
  set.seed(7)
  LL <- t(vapply(seq_len(NREG), function(i) sim_region(layouts[[nm]]), numeric(nrow(pat))))
  f <- fit_prior_bin(LL)
  Z <- sweep(LL, 2, log(pmax(f$pi, 1e-12)), `+`); Z <- Z - apply(Z, 1, max)
  po <- exp(Z); po <- po / rowSums(po)
  data.table(layout = nm, mean_meth = round(mean(layouts[[nm]]), 3),
             phi = round(mean(po %*% pat$phi), 3), s = round(mean(po %*% pat$s), 1),
             winner_post = round(median(apply(po, 1, max)), 3),
             eff_patterns = round(median(exp(-rowSums(ifelse(po > 0, po * log(po), 0)))), 2))
}, mc.cores = min(5L, parallel::detectCores())))
print(res, row.names = FALSE)
cat("truth: phi 0.6, s 100\n")

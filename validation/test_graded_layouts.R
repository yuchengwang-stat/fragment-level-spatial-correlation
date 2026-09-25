#!/usr/bin/env Rscript
## Regions whose CpGs sit in the same methylation range but not at one rate,
## e.g. 0.95, 0.9, 0.8, 0.7 along a highly methylated region.  Is the kernel still
## recovered?  Same set-up as test_bias_sources.R (8 CpGs, 60 fragments, 300
## regions, marginals counted from the reads), for two true kernels.
##
## Layouts, per range (high 0.70-0.95, hemi 0.35-0.65, low 0.05-0.30):
##   uniform   every CpG at the range's midpoint (the control)
##   random    each CpG drawn independently from the range, new draws per region
##   gradient  rates running evenly from one end of the range to the other
suppressMessages({ library(fragcorr); library(data.table); library(parallel) })
NREG <- 300L; NFRAG <- 60L; NCPG <- 8L
pat <- fc_patterns()
ranges <- list(high = c(0.70, 0.95), hemi = c(0.35, 0.65), low = c(0.05, 0.30))
layouts <- unlist(lapply(names(ranges), function(nm) {
  r <- ranges[[nm]]
  l <- list(function() rep(mean(r), NCPG),
            function() stats::runif(NCPG, r[1], r[2]),
            function() seq(r[2], r[1], length.out = NCPG))
  names(l) <- paste(nm, c("uniform", "random", "gradient"))
  l
}), recursive = FALSE)
truths <- list(c(phi = 0.6, s = 100), c(phi = 0.3, s = 300))

sim_region <- function(p, PHI, S) {
  gaps <- sample(20:120, NCPG - 1L, TRUE); pos <- c(0, cumsum(gaps))
  Sig <- PHI * exp(-abs(outer(pos, pos, "-")) / S); diag(Sig) <- 1
  C <- sweep(matrix(stats::rnorm(NFRAG * NCPG), NFRAG, NCPG) %*% chol(Sig), 2, stats::qnorm(p), `<=`)
  s <- apply(C, 1, function(v) paste(ifelse(v, "C", "T"), collapse = ""))
  tb <- table(s); n <- length(tb)
  L <- loglik_all_patterns(names(tb), rep(paste(gaps, collapse = ","), n),
                           rep(paste(colSums(C), collapse = ","), n),
                           rep(paste(rep(NFRAG, NCPG), collapse = ","), n),
                           pat$phi, pat$s, print_every = 0L)
  colSums(L * as.vector(tb), na.rm = TRUE)
}

jobs <- CJ(layout = names(layouts), truth = seq_along(truths), sorted = FALSE)
res <- rbindlist(mclapply(seq_len(nrow(jobs)), function(j) {
  lay <- jobs$layout[j]; tr <- truths[[jobs$truth[j]]]
  set.seed(100L + j)
  LL <- t(vapply(seq_len(NREG), function(i) sim_region(layouts[[lay]](), tr[["phi"]], tr[["s"]]),
                 numeric(nrow(pat))))
  f <- fit_prior_bin(LL)
  Z <- sweep(LL, 2, log(pmax(f$pi, 1e-12)), `+`); Z <- Z - apply(Z, 1, max)
  po <- exp(Z); po <- po / rowSums(po)
  data.table(truth = sprintf("phi %.1f, s %g", tr[["phi"]], tr[["s"]]), layout = lay,
             phi = round(mean(po %*% pat$phi), 3), s = round(mean(po %*% pat$s), 1),
             winner_post = round(stats::median(apply(po, 1, max)), 3))
}, mc.cores = min(6L, parallel::detectCores())))

for (t in unique(res$truth)) {
  cat("\n======== truth:", t, "========\n")
  print(res[truth == t, .(layout, phi, s, winner_post)], row.names = FALSE)
}

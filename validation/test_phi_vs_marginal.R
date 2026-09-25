#!/usr/bin/env Rscript
## Does the estimated correlation depend on the methylation level it sits on?
##
## phi is defined as a parameter of the latent Gaussian, conditional on each
## site's own marginal rate, so in principle it is free of that rate.  In
## practice a marginal near 0 or 1 makes the joint table nearly degenerate: the
## likelihood over phi flattens and the posterior is carried by the prior.  If
## that bias is strong enough, the monotone fall of phi across methylation bins
## seen in the real data would be an artefact of the estimator rather than a
## property of the genome.
##
## This generates fragments from the model with phi and s FIXED, varying only the
## marginal methylation, and asks what the pipeline returns.  A flat answer means
## the estimator is doing what it claims.
##
## NOTE (2026-09-25): this test draws each site's M / N independently of the
## fragments it simulates, which a real .beta never does -- its counts come from
## the same reads as the .pat.  That alone produces most of the ~0.2 drop in phi
## it reports at p = 0.05 / 0.95.  test_bias_sources.R counts the marginals from
## the simulated reads, as the real data does, and recovers phi within 0.03.
suppressMessages({ library(fragcorr); library(data.table) })

a <- list(regions = "400", frags = "60", depth = "30", phi = "0.6", s = "100", seed = "1")
args <- commandArgs(trailingOnly = TRUE); i <- 1
while (i <= length(args)) {
  k <- sub("^--", "", args[i])
  if (i < length(args) && !startsWith(args[i + 1], "--")) { a[[k]] <- args[i + 1]; i <- i + 2 }
  else { a[[k]] <- "true"; i <- i + 1 }
}
set.seed(as.integer(a$seed))
PHI <- as.numeric(a$phi); S <- as.numeric(a$s)
NREG <- as.integer(a$regions); NFRAG <- as.integer(a$frags); DEPTH <- as.integer(a$depth)
pat <- fc_patterns()

## one region's worth of fragments at a given marginal methylation
simulate_region <- function(p_target, n_cpg = 8L, n_frag = NFRAG, depth = DEPTH) {
  gaps <- sample(20:120, n_cpg - 1L, TRUE)
  pos  <- c(0L, cumsum(gaps))
  Sig  <- PHI * exp(-abs(outer(pos, pos, "-")) / S)
  diag(Sig) <- 1
  L <- tryCatch(chol(Sig), error = function(e) chol(Sig + diag(1e-8, n_cpg)))
  ## every site in the region sits at the same marginal, which is what a
  ## uniformly hyper- or hypo-methylated region looks like
  thr <- qnorm(p_target)
  Z <- matrix(rnorm(n_frag * n_cpg), n_frag, n_cpg) %*% L
  meth <- Z <= thr                       # TRUE = C
  frag <- apply(meth, 1, function(v) paste(ifelse(v, "C", "T"), collapse = ""))
  ## the reference counts the pipeline would read out of a .beta: the same
  ## marginal, at a realistic depth
  N <- rep(depth, n_cpg)
  M <- rbinom(n_cpg, depth, p_target)
  data.table(fragment = frag,
             dist_sub = paste(gaps, collapse = ","),
             M = paste(M, collapse = ","), N = paste(N, collapse = ","),
             count = 1L)
}

levels_p <- c(0.05, 0.15, 0.30, 0.50, 0.70, 0.85, 0.95)
cat(sprintf("truth: phi = %.2f, s = %g   |   %d regions x %d fragments at each level\n\n",
            PHI, S, NREG, NFRAG))

res <- rbindlist(lapply(levels_p, function(p) {
  fr <- rbindlist(lapply(seq_len(NREG), function(k) {
    d <- simulate_region(p); d[, interval_id := k][] }))
  Lmat <- loglik_all_patterns(fr$fragment, fr$dist_sub, fr$M, fr$N,
                              pat$phi, pat$s, print_every = 0L)
  colnames(Lmat) <- pat$name
  agg <- as.data.table(Lmat)[, lapply(.SD, sum, na.rm = TRUE), by = .(interval_id = fr$interval_id)]
  ## fit one prior for this marginal level directly: fit_priors() only knows the
  ## five real methylation bins, and here every region sits at the same level
  M <- as.matrix(agg[, pat$name, with = FALSE])
  f <- fit_prior_bin(M)
  if (is.null(f)) return(data.table(p = p, regions = 0L, phi_weighted = NA_real_,
                                    s_weighted = NA_real_, phi_argmax = NA_real_))
  lp <- pmax(f$pi, 1e-12); lp <- log(lp)
  Z <- sweep(M, 2, lp, `+`); Z <- Z - apply(Z, 1, max)
  E <- exp(Z); po <- E / rowSums(E)
  ok <- rowSums(is.na(po)) == 0
  data.table(p = p,
             regions = sum(ok),
             phi_weighted = round(mean(as.vector(po[ok, ] %*% pat$phi)), 3),
             s_weighted   = round(mean(as.vector(po[ok, ] %*% pat$s)), 1),
             phi_argmax   = round(mean(pat$phi[max.col(po[ok, , drop = FALSE])]), 3))
}))

cat("======== recovered phi as the marginal moves ========\n")
print(res, row.names = FALSE)
cat(sprintf("\n  truth phi = %.2f\n", PHI))
cat(sprintf("  recovered range: %.3f to %.3f   (spread %.3f)\n",
            min(res$phi_weighted), max(res$phi_weighted),
            max(res$phi_weighted) - min(res$phi_weighted)))
d_extreme <- res[p == 0.05, phi_weighted] - res[p == 0.95, phi_weighted]
cat(sprintf("  hypo (p=0.05) minus hyper (p=0.95): %+.3f\n", d_extreme))
cat("\n  For reference, the real data showed hypo minus hyper = +0.32 to +0.48.\n")
cat("  If the number above is near zero the estimator is marginal-free and that\n")
cat("  gap is real; if it is a large positive number the gap is the estimator.\n")

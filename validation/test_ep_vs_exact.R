#!/usr/bin/env Rscript
## Is the alternating-layout bias the EP approximation?  Score the same fragments
## under a few kernels with EP and with a precise orthant probability (mvtnorm,
## Genz-Bretz at 1e-7), using the true marginals as thresholds in both.
suppressMessages({ library(fragcorr); library(mvtnorm); library(parallel) })
PHI <- 0.6; S <- 100; NCPG <- 8L
layouts <- list(alternating = rep(c(0.05, 0.90), 4), split = c(rep(0.05, 4), rep(0.90, 4)),
                uniform = rep(0.475, 8))
ks <- data.frame(phi = c(0.4, 0.6, 0.8, 1), s = 100)
one_layout <- function(nm) {
  set.seed(11); p <- layouts[[nm]]
  out <- matrix(0, 2, nrow(ks), dimnames = list(c("EP", "exact"), paste0(ks$phi, "_", ks$s)))
  for (r in 1:25) {                        # 25 regions x 20 fragments
    gaps <- sample(20:120, NCPG - 1L, TRUE); pos <- c(0, cumsum(gaps))
    Sig0 <- PHI * exp(-abs(outer(pos, pos, "-")) / S); diag(Sig0) <- 1
    C <- sweep(matrix(rnorm(20 * NCPG), 20, NCPG) %*% chol(Sig0), 2, qnorm(p), `<=`)
    fr <- apply(C, 1, function(v) paste(ifelse(v, "C", "T"), collapse = ""))
    M <- paste(as.integer(round(p * 1e6)), collapse = ","); N <- paste(rep(1000000L, NCPG), collapse = ",")
    ep <- loglik_all_patterns(fr, rep(paste(gaps, collapse = ","), 20), rep(M, 20), rep(N, 20),
                              ks$phi, ks$s, print_every = 0L)
    out["EP", ] <- out["EP", ] + colSums(ep)
    for (k in seq_len(nrow(ks))) {
      Sig <- ks$phi[k] * exp(-abs(outer(pos, pos, "-")) / ks$s[k]); diag(Sig) <- 1
      for (f in seq_len(nrow(C))) {
        lo <- ifelse(C[f, ], -Inf, qnorm(p)); hi <- ifelse(C[f, ], qnorm(p), Inf)
        pr <- pmvnorm(lower = lo, upper = hi, sigma = Sig,
                      algorithm = GenzBretz(maxpts = 2e5, abseps = 1e-7, releps = 1e-5))
        out["exact", k] <- out["exact", k] + log(max(pr, 1e-300))
      }
    }
  }
  ## log-likelihood relative to the truth (phi 0.6): > 0 means that kernel is preferred
  rel <- sweep(out, 1, out[, "0.6_100"])
  cat("\n==", nm, "(500 fragments; log-likelihood minus that of the truth, phi 0.6) ==\n")
  print(round(rel, 1))
}
invisible(lapply(names(layouts), one_layout))

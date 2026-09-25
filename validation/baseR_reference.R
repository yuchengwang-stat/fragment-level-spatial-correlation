#!/usr/bin/env Rscript
## The whole method in base R, no packages.
##
## This is the pipeline written out as plainly as it can be, so the logic can be
## read and checked without Rcpp, data.table or mixsqp in the way.  It is slow --
## the real one runs the same model in C++ with an EP approximation to the
## orthant probability -- but every step corresponds one-to-one.
##
##   1. a region is a short stretch of CpGs; a fragment is one read over it,
##      giving C / T at each covered CpG
##   2. each CpG has its own methylation rate p_i, read off the .beta counts
##   3. under a latent Gaussian, C means Z_i <= qnorm(p_i); the correlation of
##      Z is rho(d) = phi * exp(-d/s)
##   4. P(one fragment | phi, s) is an orthant probability of that Gaussian
##   5. multiply over fragments -> region likelihood under each of 21 patterns
##   6. fit mixture weights over patterns (the prior) by EM
##   7. posterior over patterns per region
##
## The only thing the marginals p_i do is set the thresholds.  phi and s never
## touch them.  That is what makes phi a correlation net of methylation rate.

## ---------------------------------------------------------------- the 21 patterns
patterns <- function() {
  phi <- c(rep(c(1, .8, .6, .4, .2, .1), 3), 1, .8, .6)
  s   <- c(rep(30, 6), rep(100, 6), rep(300, 6), rep(1200, 3))
  data.frame(phi = phi, s = s, name = paste0(phi, "_", s), stringsAsFactors = FALSE)
}

## ------------------------------------------------- P(fragment) by GHK simulation
##
## The probability of seeing this exact C/T string is the mass the correlated
## Gaussian puts in a box: coordinate i is below its threshold where a C was
## read, above it where a T was.  There is no closed form past two dimensions, so
## it is simulated.  GHK does it by drawing each coordinate from the truncated
## normal implied by the ones already drawn, and averaging the product of the
## one-dimensional probabilities -- far steadier than counting hits.
frag_prob <- function(ct, pos, p, phi, s, R = 2000L) {
  keep <- ct %in% c("C", "T")          # "." means the read did not cover it
  ct <- ct[keep]; pos <- pos[keep]; p <- p[keep]
  K <- length(ct)
  if (K == 0L) return(1)
  p <- pmin(pmax(p, 1e-4), 1 - 1e-4)   # a rate of exactly 0 or 1 has no threshold
  thr <- qnorm(p)
  lo <- ifelse(ct == "C", -Inf, thr)   # C: Z <= thr,  T: Z > thr
  hi <- ifelse(ct == "C", thr,  Inf)
  if (K == 1L) return(pnorm(hi) - pnorm(lo))

  Sig <- phi * exp(-abs(outer(pos, pos, "-")) / s)
  diag(Sig) <- 1
  L <- t(chol(Sig))                    # lower triangular, Sig = L %*% t(L)

  v <- matrix(0, R, K)                 # the standardised draws
  w <- rep(1, R)                       # running product of the slice widths
  for (j in seq_len(K)) {
    off <- if (j == 1L) 0 else as.vector(v[, seq_len(j - 1), drop = FALSE] %*% L[j, seq_len(j - 1)])
    a <- (lo[j] - off) / L[j, j]
    b <- (hi[j] - off) / L[j, j]
    Fa <- pnorm(a); Fb <- pnorm(b)
    w <- w * (Fb - Fa)
    u <- runif(R, Fa, Fb)              # inverse-cdf draw from the truncated normal
    v[, j] <- qnorm(pmin(pmax(u, 1e-12), 1 - 1e-12))
  }
  mean(w)
}

## ------------------------------------------------ region log-likelihood, 21 columns
region_loglik <- function(frags, pos, p, pat, eps = 1e-12) {
  out <- numeric(nrow(pat))
  for (k in seq_len(nrow(pat))) {
    lp <- 0
    for (f in seq_along(frags)) {
      ct <- strsplit(frags[f], "")[[1]]
      ## a read starting partway in: line it up on the region's CpGs
      lp <- lp + log(max(frag_prob(ct, pos[seq_along(ct)], p[seq_along(ct)],
                                   pat$phi[k], pat$s[k]), eps))
    }
    out[k] <- lp
  }
  out
}

## ------------------------------------------------------- the prior, fitted by EM
##
## mixsqp and this solve the same problem: maximise sum_i log(sum_k pi_k L_ik)
## over the simplex.  EM is three lines and converges from anywhere; mixsqp is
## just faster.  Working in logs keeps a region with tiny likelihoods from
## underflowing to 0/0.
fit_prior <- function(LL, iter = 500L, tol = 1e-8) {
  K <- ncol(LL); pi <- rep(1 / K, K)
  for (it in seq_len(iter)) {
    Z <- sweep(LL, 2, log(pmax(pi, 1e-300)), `+`)
    Z <- Z - apply(Z, 1, max)
    W <- exp(Z); W <- W / rowSums(W)      # responsibilities
    new <- colMeans(W)
    if (max(abs(new - pi)) < tol) { pi <- new; break }
    pi <- new
  }
  pi
}

## ------------------------------------------------------------------- the posterior
posterior <- function(LL, pi) {
  Z <- sweep(LL, 2, log(pmax(pi, 1e-300)), `+`)
  Z <- Z - apply(Z, 1, max)
  E <- exp(Z)
  E / rowSums(E)
}

## =============================================================== demonstration
## Simulate three kinds of region -- hypo, hemi, hyper -- all with the SAME true
## correlation, and check the fitted phi does not follow the methylation level.
if (sys.nframe() == 0L) {
  set.seed(1)
  pat <- patterns()
  TRUE_PHI <- 0.6; TRUE_S <- 100

  simulate <- function(p_level, n_region = 30L, n_frag = 40L, n_cpg = 8L) {
    lapply(seq_len(n_region), function(r) {
      pos <- c(0, cumsum(sample(20:120, n_cpg - 1, TRUE)))
      Sig <- TRUE_PHI * exp(-abs(outer(pos, pos, "-")) / TRUE_S); diag(Sig) <- 1
      Zm <- matrix(rnorm(n_frag * n_cpg), n_frag, n_cpg) %*% chol(Sig)
      ct <- ifelse(Zm <= qnorm(p_level), "C", "T")
      list(frags = apply(ct, 1, paste, collapse = ""),
           pos = pos, p = rep(p_level, n_cpg))
    })
  }

  for (lev in c(0.05, 0.50, 0.95)) {
    regs <- simulate(lev)
    LL <- t(vapply(regs, function(r) region_loglik(r$frags, r$pos, r$p, pat),
                   numeric(nrow(pat))))
    pi <- fit_prior(LL)
    po <- posterior(LL, pi)
    cat(sprintf("methylation %.2f   fitted phi = %.3f   fitted s = %5.1f bp\n",
                lev, mean(po %*% pat$phi), mean(po %*% pat$s)))
  }
  cat(sprintf("\ntruth: phi = %.2f, s = %g\n", TRUE_PHI, TRUE_S))
}

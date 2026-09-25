#!/usr/bin/env Rscript
## Where does the phi bias at extreme methylation come from?
##
## test_phi_vs_marginal.R recovers phi = 0.35-0.38 at p = 0.05 / 0.95 against a
## truth of 0.6.  But it draws each site's reference counts independently of the
## fragments it simulates, and a real .beta is not like that: its M / N are
## counted from the same reads as the .pat.  So the first question is whether the
## bias belongs to the method or to the test.  One thing changes at a time:
##
##   indep   marginals drawn independently of the fragments, depth 30 (the old test)
##   reads   marginals counted from the simulated fragments themselves (real data)
##   oracle  the true marginal, no estimation error
##
## Everything else is the pipeline: the EP engine, a mixsqp prior per level, and
## the posterior over the same 21 patterns.  The prior's own mean is reported as
## well as the average posterior, because with a correct likelihood that mean
## should sit near the truth however weak each region is on its own.
suppressMessages({ library(fragcorr); library(data.table); library(parallel) })

a <- list(regions = "300", frags = "60", ncpg = "8", depth = "30",
          phi = "0.6", s = "100", seed = "1", cores = "4", out = "")
args <- commandArgs(trailingOnly = TRUE); i <- 1
while (i <= length(args)) {
  k <- sub("^--", "", args[i])
  if (i < length(args) && !startsWith(args[i + 1], "--")) { a[[k]] <- args[i + 1]; i <- i + 2 }
  else { a[[k]] <- "true"; i <- i + 1 }
}
PHI <- as.numeric(a$phi); S <- as.numeric(a$s)
NREG <- as.integer(a$regions); NFRAG <- as.integer(a$frags)
NCPG <- as.integer(a$ncpg); DEPTH <- as.integer(a$depth)
pat <- fc_patterns()
levels_p <- c(0.05, 0.15, 0.30, 0.50, 0.70, 0.85, 0.95)
how_all <- c("indep", "reads", "oracle")

sim_region <- function(p) {
  gaps <- sample(20:120, NCPG - 1L, TRUE); pos <- c(0, cumsum(gaps))
  Sig <- PHI * exp(-abs(outer(pos, pos, "-")) / S); diag(Sig) <- 1
  C <- (matrix(rnorm(NFRAG * NCPG), NFRAG, NCPG) %*% chol(Sig)) <= qnorm(p)
  list(gaps = gaps, C = C, p = p)
}

marginals <- function(r, how) switch(how,
  indep  = list(M = rbinom(NCPG, DEPTH, r$p), N = rep(DEPTH, NCPG)),
  reads  = list(M = colSums(r$C), N = rep(NFRAG, NCPG)),
  ## integers, so paste() never writes 1e+06, which the C++ parser reads as 1
  oracle = list(M = rep(as.integer(round(r$p * 1e6)), NCPG), N = rep(1000000L, NCPG)))

## identical fragments carry identical likelihoods, so each distinct one is
## scored once and weighted by how often it occurs
region_loglik <- function(r, mg) {
  s <- apply(r$C, 1, function(v) paste(ifelse(v, "C", "T"), collapse = ""))
  tb <- table(s); n <- length(tb)
  L <- loglik_all_patterns(names(tb), rep(paste(r$gaps, collapse = ","), n),
                           rep(paste(mg$M, collapse = ","), n),
                           rep(paste(mg$N, collapse = ","), n),
                           pat$phi, pat$s, print_every = 0L)
  colSums(L * as.vector(tb), na.rm = TRUE)
}

## the same simulated regions go through every variant, so the variants differ
## only in the marginals
set.seed(as.integer(a$seed))
regs <- lapply(levels_p, function(p) lapply(seq_len(NREG), function(k) sim_region(p)))
names(regs) <- levels_p

jobs <- CJ(p = levels_p, how = how_all)
res <- rbindlist(mclapply(seq_len(nrow(jobs)), function(j) {
  p <- jobs$p[j]; how <- jobs$how[j]
  set.seed(1000L + j)
  LL <- t(vapply(regs[[as.character(p)]], function(r) region_loglik(r, marginals(r, how)),
                 numeric(nrow(pat))))
  f <- fit_prior_bin(LL)
  lp <- log(pmax(f$pi, 1e-12))
  Z <- sweep(LL, 2, lp, `+`); Z <- Z - apply(Z, 1, max); po <- exp(Z); po <- po / rowSums(po)
  data.table(p = p, how = how,
             phi_post = mean(po %*% pat$phi), s_post = mean(po %*% pat$s),
             phi_prior = sum(f$pi * pat$phi),
             phi_argmax = mean(pat$phi[max.col(po)]),
             p_truth_pattern = mean(po[, which(pat$phi == PHI & pat$s == S)]))
}, mc.cores = as.integer(a$cores)))

cat(sprintf("truth: phi = %.2f, s = %g  |  %d regions x %d fragments x %d CpGs per level\n\n",
            PHI, S, NREG, NFRAG, NCPG))
w <- dcast(res, p ~ how, value.var = "phi_post")
setcolorder(w, c("p", how_all))
cat("======== average posterior-mean phi ========\n")
print(w[, (how_all) := lapply(.SD, round, 3), .SDcols = how_all], row.names = FALSE)
cat("\n======== prior mean phi (mixsqp) ========\n")
print(dcast(res, p ~ how, value.var = "phi_prior")[, c("p", how_all), with = FALSE][
  , (how_all) := lapply(.SD, round, 3), .SDcols = how_all], row.names = FALSE)
cat("\n======== average posterior-mean s ========\n")
print(dcast(res, p ~ how, value.var = "s_post")[, c("p", how_all), with = FALSE][
  , (how_all) := lapply(.SD, round, 1), .SDcols = how_all], row.names = FALSE)
if (nzchar(a$out)) fwrite(res, a$out, sep = "\t")

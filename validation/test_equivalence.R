#!/usr/bin/env Rscript
## The rebuilt engine must reproduce the original's numbers, not merely agree
## with them closely.  Both are compiled here and run on the same fragments.
## The original is compiled from validation/reference/; the rebuild is the
## installed package.
suppressMessages({ library(Rcpp); library(fragcorr) })
here <- dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))

cat("compiling the original ...\n")
sourceCpp(file.path(here, "reference", "fragment_lik.cpp"), verbose = FALSE)

## ---- the 21 patterns ---------------------------------------------------------
phi <- c(1, .8, .6, .4, .2, .1,  1, .8, .6, .4, .2, .1,
         1, .8, .6, .4, .2, .1,  1, .8, .6)
s   <- c(rep(30, 6), rep(100, 6), rep(300, 6), rep(1200, 3))
stopifnot(length(phi) == 21L, length(s) == 21L)

## ---- synthetic fragments covering every code path ----------------------------
set.seed(11)
make_frag <- function(L, depth_lo = 1L) {
  ch <- sample(c("C", "T", "."), L, TRUE, prob = c(.45, .45, .10))
  N  <- sample(seq(depth_lo, 60L), L, TRUE)
  M  <- vapply(N, function(nn) sample.int(nn + 1L, 1L) - 1L, 0L)
  d  <- if (L > 1) sample(c(1:50, 60:400, 500:3000), L - 1, TRUE) else integer(0)
  list(fragment = paste(ch, collapse = ""),
       dist_sub = paste(d, collapse = ","),
       M = paste(M, collapse = ","), N = paste(N, collapse = ","))
}
frs <- c(
  lapply(1:60,  function(i) make_frag(sample(1:2, 1))),        # K = 0,1,2
  lapply(1:120, function(i) make_frag(sample(3:8, 1))),        # small K
  lapply(1:40,  function(i) make_frag(sample(9:16, 1))),       # larger K
  lapply(1:20,  function(i) make_frag(sample(2:6, 1), 1L))     # shallow depth
)
FR <- data.frame(
  fragment = vapply(frs, `[[`, "", "fragment"),
  dist_sub = vapply(frs, `[[`, "", "dist_sub"),
  M = vapply(frs, `[[`, "", "M"), N = vapply(frs, `[[`, "", "N"),
  stringsAsFactors = FALSE)
cat(sprintf("%d fragments, lengths %d..%d\n", nrow(FR),
            min(nchar(FR$fragment)), max(nchar(FR$fragment))))

## ---- original: one call per pattern ------------------------------------------
t0 <- Sys.time()
ORIG <- vapply(seq_along(phi), function(i) {
  exp_table <- phi[i] * exp(-(0:7000) / s[i])
  compute_loglik_cpp_coord_string(FR$fragment, FR$dist_sub, FR$M, FR$N,
                                  exp_table, print_every = 0L)
}, numeric(nrow(FR)))
t_orig <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

## ---- rebuild, default settings: must match the original exactly --------------
t0 <- Sys.time()
NEW_ep <- loglik_all_patterns(FR$fragment, FR$dist_sub, FR$M, FR$N, phi, s,
                              print_every = 0L)
t_new_ep <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

## ---- rebuild with the closed-form bivariate shortcut -------------------------
t0 <- Sys.time()
NEW <- loglik_all_patterns(FR$fragment, FR$dist_sub, FR$M, FR$N, phi, s,
                           use_bvn = TRUE, print_every = 0L)
t_new <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

chk <- function(label, A, B) {
  na_same <- identical(is.na(A), is.na(B))
  ok <- !is.na(A) & !is.na(B)
  d <- if (any(ok)) max(abs(A[ok] - B[ok])) else 0
  ident <- na_same && isTRUE(all.equal(A[ok], B[ok], tolerance = 0))
  cat(sprintf("  %-46s %s  (max |diff| %.3g, NA pattern %s)\n", label,
              if (ident) "IDENTICAL" else if (d <= 5e-4) "equal to 3dp" else "DIFFERS",
              d, if (na_same) "same" else "DIFFERENT"))
  list(ident = ident, d = d, na_same = na_same)
}

cat("\n======== agreement ========\n")
r1 <- chk("rebuild, defaults, vs original", NEW_ep, ORIG)
r2 <- chk("rebuild with use_bvn=TRUE vs original", NEW, ORIG)
cat("  (the second is EXPECTED to differ: the EP solver is approximate at K = 2,\n")
cat("   and the closed form is the accurate one.  Off by default.)\n")

## how many fragments land in each branch
Ks <- vapply(seq_len(nrow(FR)), function(i) {
  ch <- strsplit(FR$fragment[i], "")[[1]]
  N <- as.integer(strsplit(FR$N[i], ",")[[1]])
  sum(ch %in% c("C", "T") & N >= 5)
}, 0L)
cat(sprintf("\n  fragments by usable sites K:  K=0 %d, K=1 %d, K=2 %d, K>=3 %d\n",
            sum(Ks == 0), sum(Ks == 1), sum(Ks == 2), sum(Ks >= 3)))

cat("\n======== speed on this set ========\n")
cat(sprintf("  original, 21 passes          %6.2f s\n", t_orig))
cat(sprintf("  rebuild, one pass (EP only)  %6.2f s   %.1fx\n", t_new_ep, t_orig / t_new_ep))
cat(sprintf("  rebuild, one pass (+bvn)     %6.2f s   %.1fx\n", t_new, t_orig / t_new))

stopifnot(r1$na_same, r1$ident)          # the default path must be exact
stopifnot(r2$na_same)                    # the shortcut must not change which are NA
cat("\nEQUIVALENCE PASSED\n")

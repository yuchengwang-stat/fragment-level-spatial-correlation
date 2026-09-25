## Prior over patterns, and the posterior for each region.
##
## Each region carries a log-likelihood under each of the 21 patterns, summed
## over the fragments that fall in it and then over the samples of a cell type.
## The prior is the mixture weight of each pattern, fitted by mixsqp within a
## methylation bin, and the posterior is the normalised product.
##
## The whole point of the pipeline is the last object: for every region, a
## probability distribution over the 21 patterns.


## Sum region log-likelihoods across a set of per-sample tables.
##
## Each input is a data.table with chr, interval_id and one column per pattern.
## Log-likelihoods add across independent samples, so this is a plain sum.
combine_loglik <- function(tabs, pattern_names) {
  stopifnot(length(tabs) > 0L)
  all <- data.table::rbindlist(tabs, use.names = TRUE, fill = TRUE)
  all[, lapply(.SD, sum, na.rm = TRUE), by = .(chr, interval_id), .SDcols = pattern_names]
}

## Fit the mixture weights over patterns within one bin.
##
## mixsqp maximises sum_i log(sum_k pi_k L_ik) over the simplex, which is exactly
## the marginal likelihood of the regions under a shared prior.  It is given the
## log-likelihoods directly (log = TRUE) so nothing has to be exponentiated.
fit_prior_bin <- function(loglik_mat, min_pi = 0) {
  if (!nrow(loglik_mat)) return(NULL)
  ok <- stats::complete.cases(loglik_mat) & is.finite(rowSums(loglik_mat))
  if (!any(ok)) return(NULL)
  fit <- mixsqp::mixsqp(loglik_mat[ok, , drop = FALSE], log = TRUE,
                        control = list(verbose = FALSE))
  pi <- as.numeric(fit$x)
  if (min_pi > 0) { pi <- pmax(pi, min_pi); pi <- pi / sum(pi) }
  list(pi = pi, n = sum(ok))
}

## Fit one prior per bin.  Returns a matrix, bins x patterns.
fit_priors <- function(lik, pattern_names, bin_col = "bin", min_pi = 0, max_fit = 0L) {
  bins <- fc_bins()$label
  out <- matrix(NA_real_, length(bins), length(pattern_names),
                dimnames = list(bins, pattern_names))
  nobs <- stats::setNames(integer(length(bins)), bins)
  for (b in bins) {
    idx <- which(!is.na(lik[[bin_col]]) & as.character(lik[[bin_col]]) == b)
    if (!length(idx)) next
    M <- as.matrix(lik[idx, ..pattern_names])
    M <- M[subsample_rows(nrow(M), max_fit), , drop = FALSE]
    f <- fit_prior_bin(M, min_pi)
    if (is.null(f)) next
    out[b, ] <- f$pi; nobs[b] <- f$n
  }
  list(pi = out, n = nobs)
}

## Posterior over patterns for every region.
##
##   log posterior_k = loglik_k + log(pi_k), normalised across k
##
## Normalisation subtracts the row maximum before exponentiating, so a region
## whose log-likelihoods are all very negative does not underflow to 0/0.
posterior_from_prior <- function(lik, priors, pattern_names, bin_col = "bin") {
  post <- matrix(NA_real_, nrow(lik), length(pattern_names),
                 dimnames = list(NULL, paste0(pattern_names, "_post")))
  for (b in rownames(priors)) {
    if (all(is.na(priors[b, ]))) next
    idx <- which(!is.na(lik[[bin_col]]) & as.character(lik[[bin_col]]) == b)
    if (!length(idx)) next
    lp <- priors[b, ]
    lp[lp <= 0] <- 1e-12
    lp <- log(lp)
    M <- as.matrix(lik[idx, ..pattern_names])
    M <- sweep(M, 2, lp, `+`)
    M <- M - apply(M, 1, max)
    E <- exp(M)
    post[idx, ] <- E / rowSums(E)
  }
  post
}

## The reported summary for a region: which pattern wins, how sure, and what the
## winning pattern's phi and s actually are.
summarise_posterior <- function(post, patterns) {
  ok <- rowSums(is.na(post)) == 0L
  best <- rep(NA_integer_, nrow(post))
  best[ok] <- max.col(post[ok, , drop = FALSE], ties.method = "first")
  data.frame(
    best_pattern = ifelse(is.na(best), NA_character_, patterns$name[best]),
    best_phi     = ifelse(is.na(best), NA_real_,      patterns$phi[best]),
    best_s       = ifelse(is.na(best), NA_real_,      patterns$s[best]),
    best_post    = ifelse(is.na(best), NA_real_, post[cbind(seq_len(nrow(post)), best)]),
    ## how concentrated the posterior is, in units of patterns: 1 means certain,
    ## 21 means no information at all
    eff_patterns = ifelse(ok, exp(-rowSums(ifelse(post > 0, post * log(post), 0))), NA_real_),
    stringsAsFactors = FALSE)
}

## mixsqp gets slow when it is handed a very tall matrix, and the prior only
## needs to be representative, not exhaustive.  Fitting on a random subset of the
## regions in a bin gives the same weights to within sampling error at a fraction
## of the cost; `max_fit` caps it, and 0 means use every region.
subsample_rows <- function(n, max_fit, seed = 1L) {
  if (max_fit <= 0L || n <= max_fit) return(seq_len(n))
  set.seed(seed)
  sort(sample.int(n, max_fit))
}

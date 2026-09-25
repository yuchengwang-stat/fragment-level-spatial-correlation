## Region log-likelihoods for one sample, under every correlation pattern.
##
## This is step 2 of the pipeline as a function.  One pass over the sample's
## fragments produces all pattern columns, and the result is aggregated to region
## level before anything is returned: one sample becomes one small table.
##
##   reg       an indexed region table (index_regions())
##   pat       read patterns from read_pat(), ideally prefiltered to the regions
##   beta      the sample's n x 2 count matrix from read_beta_counts()
##   patterns  the (phi, s) grid; fc_patterns() by default
##   min_n     a site is used only where the sample's .beta depth is at least this
##   use_bvn   closed-form bivariate normal for two-site fragments.  Off by
##             default: the EP solver is what the original method used, and it
##             differs from the exact value at K = 2 (see src/loglik.cpp)
##   threads   forked workers for the per-fragment loop
##
## Returns a data.table with one row per region: chr, interval_id, one column per
## pattern holding the summed log-likelihood, M_sum, N_sum, meth and bin.
sample_loglik <- function(reg, pat, beta, patterns = fc_patterns(), min_n = 5L,
                          use_bvn = FALSE, threads = 1L) {
  if (max(reg$endCpG) > nrow(beta))
    stop("regions index CpG ", max(reg$endCpG), " but the .beta holds ", nrow(beta),
         " -- the region table and this sample use different CpG indices")
  fr <- fragments_in_regions(pat, reg, beta)
  if (!nrow(fr)) stop("no fragment overlaps any region")

  run_slice <- function(ix)
    loglik_all_patterns(fr$fragment[ix], fr$dist_sub[ix], fr$M[ix], fr$N[ix],
                        patterns$phi, patterns$s, min_n = as.integer(min_n),
                        use_bvn = isTRUE(use_bvn), print_every = 0L)
  threads <- as.integer(threads)
  if (threads > 1L && nrow(fr) > 1000L) {
    ## many small slices rather than one per worker: fragments differ hugely in
    ## cost (K >= 3 runs the EP solver, K <= 1 does not), so equal-sized slices
    ## would leave workers idle.  Forked, not threaded, because the EP solver
    ## allocates R objects.
    nchunk <- min(nrow(fr), threads * 8L)
    slices <- split(seq_len(nrow(fr)), cut(seq_len(nrow(fr)), nchunk, labels = FALSE))
    parts <- parallel::mclapply(slices, run_slice, mc.cores = threads, mc.preschedule = FALSE)
    bad <- vapply(parts, function(x) inherits(x, "try-error") || is.null(x), TRUE)
    if (any(bad)) stop(sum(bad), " of ", length(parts), " slices failed in the worker processes")
    L <- do.call(rbind, parts)
  } else {
    L <- run_slice(seq_len(nrow(fr)))
  }
  stopifnot(nrow(L) == nrow(fr))
  colnames(L) <- patterns$name

  ## weight each distinct fragment by the number of reads it stands for, then
  ## add up within a region
  W <- L * fr$count
  agg <- data.table::as.data.table(W)[, lapply(.SD, sum, na.rm = TRUE),
                                      by = list(chr = fr$chr, interval_id = fr$interval_id)]

  ## region methylation, used to choose the prior's bin
  mn <- fr[, list(M_sum = sum(as.numeric(unlist(strsplit(M, ",", fixed = TRUE)))),
                  N_sum = sum(as.numeric(unlist(strsplit(N, ",", fixed = TRUE))))),
           by = list(chr, interval_id)]
  agg <- merge(agg, mn, by = c("chr", "interval_id"), all.x = TRUE)
  agg[, meth := M_sum / N_sum]
  agg[, bin := assign_bin(meth)]
  data.table::setattr(agg, "n_fragments", nrow(fr))
  data.table::setattr(agg, "n_reads", sum(fr$count))
  agg[]
}

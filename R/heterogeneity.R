## How uniform is a region's methylation, per group of samples?
##
## The kernel is estimated conditional on every CpG's own rate, so a region whose
## CpGs sit at different levels is modelled correctly as long as neighbouring
## CpGs are not at opposite extremes.  Two situations still need flagging:
##
##   * A region whose CpGs mix very low and very high rates has a mean methylation
##     that describes none of them, so it lands in the wrong methylation bin.
##   * If the low and high CpGs are interleaved, the reads carry little
##     information about the kernel, and EP tilts phi upward (VALIDATION.md, 5).
##
## These functions measure it from per-CpG counts pooled over a group's samples.
##
## Per region, over CpGs with pooled depth >= min_n:
##   n_cov            CpGs used
##   meth_all         mean of the per-CpG rates, over every CpG in the region
##   meth_all_sd      their SD
##   meth_all_iqr     their interquartile range
##   meth_target, meth_target_sd, meth_target_iqr
##                    the same over the listed (target) CpGs only (SD and IQR are
##                    NA with fewer than two covered targets)
##   meth_range       max - min, over every CpG
##   meth_jump    mean |p_(j+1) - p_j| over consecutive CpGs
##   r2_split     share of the variance explained by the best single cut into two
##                contiguous pieces: 1 = two flat halves, low = interleaved
##   n_low, n_high  CpGs at <= low / >= high
##   target_lowhigh the listed CpGs include one at <= low and one at >= high
##   layout       "uniform"          SD of the rates <= sd_uniform
##                "graded"           not uniform, but not a low/high mixture
##                "mixed_split"      >= 2 CpGs at each extreme, r2_split >= r2_cut
##                "mixed_interleaved" >= 2 CpGs at each extreme, r2_split < r2_cut
##                "too_few_cpgs"     fewer than 3 CpGs with enough depth

## Statistics from per-CpG pooled counts.
##
##   M, N        matrices, one row per CpG (in genome order within each region),
##               one column per group
##   region_id   the region of each row
##   target      logical, whether each row is a listed CpG (optional)
heterogeneity_stats <- function(M, N, region_id, target = NULL, min_n = 20L,
                                low = 0.2, high = 0.8, sd_uniform = 0.10, r2_cut = 0.75) {
  M <- as.matrix(M); N <- as.matrix(N)
  stopifnot(identical(dim(M), dim(N)), length(region_id) == nrow(M))
  if (is.null(colnames(M))) colnames(M) <- paste0("group", seq_len(ncol(M)))
  if (is.null(target)) target <- rep(FALSE, nrow(M))
  P <- ifelse(N >= min_n, M / N, NA_real_)

  one <- function(p, tgt) {
    ok <- !is.na(p); q <- p[ok]; n <- length(q)
    tq <- p[ok & tgt]
    spread <- function(v) if (length(v) >= 2L) c(stats::sd(v), stats::IQR(v)) else c(NA_real_, NA_real_)
    sa <- spread(q); st <- spread(tq)
    tlh <- length(tq) >= 2L && any(tq <= low) && any(tq >= high)
    base <- list(n_cov = n, meth_all = if (n) mean(q) else NA_real_,
                 meth_all_sd = sa[1], meth_all_iqr = sa[2],
                 meth_target = if (length(tq)) mean(tq) else NA_real_,
                 meth_target_sd = st[1], meth_target_iqr = st[2])
    if (n < 3L)
      return(c(base, list(meth_range = NA_real_, meth_jump = NA_real_, r2_split = NA_real_,
                          n_low = NA_integer_, n_high = NA_integer_, target_lowhigh = tlh,
                          enough = FALSE)))
    sst <- sum((q - mean(q))^2)
    cs <- cumsum(q); cs2 <- cumsum(q^2); k <- seq_len(n - 1L)
    sse <- (cs2[k] - cs[k]^2 / k) + ((cs2[n] - cs2[k]) - (cs[n] - cs[k])^2 / (n - k))
    c(base, list(meth_range = max(q) - min(q), meth_jump = mean(abs(diff(q))),
                 r2_split = if (sst > 1e-12) 1 - min(sse) / sst else 1,
                 n_low = sum(q <= low), n_high = sum(q >= high), target_lowhigh = tlh,
                 enough = TRUE))
  }
  out <- data.table::rbindlist(lapply(colnames(M), function(g) {
    d <- data.table::data.table(region_id = region_id, p = P[, g], tgt = target)
    d[, one(p, tgt), by = region_id][, celltype := g]
  }))
  out[, layout := data.table::fifelse(!enough, "too_few_cpgs",
                   data.table::fifelse(n_low >= 2L & n_high >= 2L,
                     data.table::fifelse(r2_split >= r2_cut, "mixed_split", "mixed_interleaved"),
                   data.table::fifelse(meth_all_sd <= sd_uniform, "uniform", "graded")))]
  out[, enough := NULL]
  data.table::setcolorder(out, c("region_id", "celltype"))
  out[]
}

## Pool per-CpG counts over the samples of each group, straight from the .beta
## files, for every CpG inside the regions.
##
##   reg     indexed region table (step 1)
##   sheet   data frame with columns sample, celltype, beta
## Returns list(region_id, cpg_index, M, N) with one row per CpG and one column
## per cell type.
pooled_region_counts <- function(reg, sheet) {
  ix  <- unlist(Map(seq.int, reg$startCpG, reg$endCpG))
  rid <- rep(reg$interval_id, reg$endCpG - reg$startCpG + 1L)
  cts <- sort(unique(sheet$celltype))
  M <- matrix(0L, length(ix), length(cts), dimnames = list(NULL, cts)); N <- M
  for (k in seq_len(nrow(sheet))) {
    f <- sheet$beta[k]
    b <- readBin(f, "integer", n = file.size(f), size = 1L, signed = FALSE)
    if (2L * max(ix) > length(b)) stop(f, " is shorter than the region index")
    ct <- sheet$celltype[k]
    M[, ct] <- M[, ct] + b[2L * ix - 1L]
    N[, ct] <- N[, ct] + b[2L * ix]
  }
  list(region_id = rid, cpg_index = ix, M = M, N = N)
}

## Which rows of pooled_region_counts() are listed CpGs (a listed position at
## the G of an indexed CpG counts for that CpG; see match_targets()).
target_rows <- function(reg, cpg, ref) {
  n <- sum(reg$endCpG - reg$startCpG + 1L)
  seq_len(n) %in% match_targets(reg, cpg, ref)$row
}

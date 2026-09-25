## The correlation-versus-distance curve for each methylation bin, with a band.
##
## phi alone is the curve's height at zero distance and s alone is how fast it
## falls; neither is the curve itself.  The posterior carries both, so the
## summary is
##
##   rho_i(d) = sum_k post_ik * phi_k * exp(-d / s_k)
##
## averaged over the regions of a bin.  That is a mixture of exponentials, not
## the exponential of an average: plugging mean(phi) and mean(s) into a single
## exponential gives a visibly different tail.
##
## The band is a percentile bootstrap across cell types, not across regions.
## Regions within a cell type are not independent replicates, whereas "does this
## hold in every cell type" is the question the band should answer.
##
##   post         a named list of posterior tables (as for region_tables())
##   dmax, dstep  the distance grid, in bp
##   boot, seed   bootstrap replicates and seed
##   min_regions  a bin is averaged only if a cell type has at least this many
##                regions in it
decay_curves <- function(post, dmax = 400, dstep = 4, boot = 4000L, seed = 1L,
                         min_regions = 20L) {
  pat  <- fc_patterns()
  bins <- fc_bins()$label
  D    <- seq(0, dmax, by = dstep)
  pn   <- paste0(pat$name, "_post")
  cts  <- names(post)

  ## B[k, j]: what pattern k contributes to the curve at distance D[j]
  B <- pat$phi * exp(-outer(pat$s, D, function(s, d) d / s))

  curves <- array(NA_real_, c(length(cts), length(bins), length(D)),
                  dimnames = list(cts, bins, D))
  phis <- ss <- ns <- matrix(NA_real_, length(cts), length(bins), dimnames = list(cts, bins))
  for (i in seq_along(cts)) {
    x  <- data.table::as.data.table(post[[i]])
    P  <- as.matrix(x[, pn, with = FALSE])
    ok <- stats::complete.cases(P) & !is.na(x$bin)
    P  <- P[ok, , drop = FALSE]; bb <- as.character(x$bin[ok])
    C  <- P %*% B
    pv <- as.vector(P %*% pat$phi); sv <- as.vector(P %*% pat$s)
    for (b in bins) {
      j <- which(bb == b)
      if (length(j) < min_regions) next
      curves[i, b, ] <- colMeans(C[j, , drop = FALSE])
      phis[i, b] <- mean(pv[j]); ss[i, b] <- mean(sv[j]); ns[i, b] <- length(j)
    }
  }

  set.seed(seed)
  band <- data.table::rbindlist(lapply(bins, function(b) {
    M <- matrix(curves[, b, ], nrow = length(cts))
    M <- M[stats::complete.cases(M), , drop = FALSE]
    if (!nrow(M)) return(NULL)
    bs <- replicate(boot, colMeans(M[sample.int(nrow(M), nrow(M), TRUE), , drop = FALSE]))
    bs <- matrix(bs, nrow = length(D))
    q  <- apply(bs, 1, stats::quantile, c(0.025, 0.975))
    data.table::data.table(bin = b, d = D, mean = colMeans(M),
                           lo = q[1, ], hi = q[2, ], n_celltype = nrow(M))
  }))
  summ <- data.table::data.table(bin = bins,
                                 phi = round(colMeans(phis, na.rm = TRUE), 3),
                                 s = round(colMeans(ss, na.rm = TRUE), 1),
                                 regions = round(colMeans(ns, na.rm = TRUE)),
                                 celltypes = colSums(!is.na(phis)))
  list(band = band, summary = summ, curves = curves, phi = phis, s = ss, n = ns,
       patterns = pat, D = D)
}

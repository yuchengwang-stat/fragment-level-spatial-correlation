## Split regions whose target CpGs sit at opposite methylation extremes.
##
## Regions are defined once, for every cell type.  A region built by padding and
## merging can still contain one target CpG that is unmethylated and another that
## is methylated; its mean methylation then describes neither, and the two
## targets arguably belong to different regions.  refine_regions() cuts such
## regions between targets, using methylation pooled per cell type, so that the
## split is shared by all cell types.
##
## A run of consecutive targets "conflicts" in a cell type when, among its targets
## with pooled depth >= min_n, one is <= low and another >= high.  A run that
## conflicts in at least min_celltypes cell types is cut between two consecutive
## targets, at the cut leaving the fewest conflicts on the two sides (ties: the
## widest bp gap), and each side is checked again.  Requiring several cell types
## keeps a single noisy estimate from triggering a cut.
##
## A piece spans its targets padded by `pad`, never beyond the original region,
## and two neighbouring pieces meet at the midpoint between them.  Pieces left
## with a single target are kept; region_tables() marks them.

## Listed CpGs inside each region, matched to rows of pooled_region_counts().
## A listed position that is the G of an indexed CpG (position - 1) is matched to
## that CpG.  Unmatched targets are kept with row = NA: they still count as
## targets, they just carry no methylation.
match_targets <- function(reg, cpg, ref) {
  ix  <- unlist(Map(seq.int, reg$startCpG, reg$endCpG))
  rid <- rep(reg$interval_id, reg$endCpG - reg$startCpG + 1L)
  key <- paste(ref$chr[ix], ref$pos[ix])
  cp <- unique(data.table::as.data.table(cpg)[, list(chr, pos)])
  r <- match(paste(cp$chr, cp$pos), key)
  g <- match(paste(cp$chr, cp$pos - 1L), key)
  cp[, row := data.table::fifelse(is.na(r), g, r)]
  ## the region of every listed CpG, by coordinates (also for unmatched ones)
  cp[, region_id := NA_integer_]
  for (ch in unique(reg$chr)) {
    rr <- reg[reg$chr == ch, , drop = FALSE]; rr <- rr[order(rr$start), , drop = FALSE]
    w <- which(cp$chr == ch); if (!length(w)) next
    k <- findInterval(cp$pos[w], rr$start); hit <- k > 0L
    hit[hit] <- cp$pos[w][hit] <= rr$end[k[hit]]
    cp[w[hit], region_id := rr$interval_id[k[hit]]]
  }
  cp <- cp[!is.na(region_id)]
  data.table::setorder(cp, region_id, pos)
  cp[]
}

refine_regions <- function(reg, targets, P, min_celltypes = 5L, low = 0.2, high = 0.8,
                           pad = 250L) {
  reg <- data.table::as.data.table(reg)
  conflicts <- function(rows) {
    rows <- rows[!is.na(rows)]
    if (length(rows) < 2L) return(0L)
    X <- P[rows, , drop = FALSE]
    sum(colSums(X <= low, na.rm = TRUE) > 0 & colSums(X >= high, na.rm = TRUE) > 0)
  }
  split_run <- function(i) {                      # i indexes the region's targets
    rows <- tg$row[i]
    if (length(i) < 2L || conflicts(rows) < min_celltypes) return(list(i))
    m <- length(i)
    sc <- vapply(seq_len(m - 1L), function(c) conflicts(rows[1:c]) + conflicts(rows[(c + 1):m]), 0L)
    gap <- diff(tg$pos[i]); best <- which(sc == min(sc)); c <- best[which.max(gap[best])]
    c(split_run(i[1:c]), split_run(i[(c + 1):m]))
  }

  out <- vector("list", nrow(reg))
  for (k in seq_len(nrow(reg))) {
    r <- reg[k]
    tg <- targets[region_id == r$interval_id]
    pieces <- if (nrow(tg) >= 2L) split_run(seq_len(nrow(tg))) else list(seq_len(nrow(tg)))
    if (length(pieces) <= 1L) {
      out[[k]] <- data.table::data.table(chr = r$chr, start = r$start, end = r$end,
                                         n_cpg = r$n_cpg, parent_region = r$interval_id)
      next
    }
    lo <- vapply(pieces, function(i) min(tg$pos[i]), 0); hi <- vapply(pieces, function(i) max(tg$pos[i]), 0)
    st <- pmax(lo - pad, r$start); en <- pmin(hi + pad, r$end)
    mid <- (hi[-length(hi)] + lo[-1]) %/% 2
    en[-length(en)] <- pmin(en[-length(en)], mid); st[-1] <- pmax(st[-1], mid + 1)
    out[[k]] <- data.table::data.table(chr = r$chr, start = as.integer(st), end = as.integer(en),
                                       n_cpg = lengths(pieces), parent_region = r$interval_id)
  }
  res <- data.table::rbindlist(out)
  res[, width := end - start + 1L]
  res[, interval_id := seq_len(.N)]
  as.data.frame(res)
}

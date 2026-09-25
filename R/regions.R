## Regions from a list of CpG positions.
##
## The correlation model is fitted per region, so a region has to be small enough
## that one (phi, s) plausibly describes all of it, and large enough to contain
## several CpGs on the same fragment.  Padding each CpG by a fixed +/- and taking
## the union does both, except where CpGs are dense: there the union runs away
## and produces one enormous region.  So the width is capped and anything over it
## is cut into near-equal pieces.
##
## Everything here is in base-pair coordinates on one genome build.  CpG indices
## (the .beta / .pat row numbers) are attached later, by build_regions_indexed().

## Merge padded CpG positions into regions.
##
##   pos        integer positions, one per CpG, within a single chromosome
##   pad        half-width added to each CpG before merging
##   max_width  a region wider than this is split
##   min_cpgs   regions with fewer CpGs than this are dropped
##
## Returns a data.frame with start, end, n_cpg -- sorted, non-overlapping.
merge_positions <- function(pos, pad = 250L, max_width = 1000L, min_cpgs = 2L) {
  pos <- sort(unique(as.integer(pos)))
  if (!length(pos)) return(data.frame(start = integer(), end = integer(), n_cpg = integer()))

  ## A new cluster starts wherever the padded intervals stop touching.  Two CpGs
  ## more than 2*pad apart can never share a region.
  brk <- c(TRUE, diff(pos) > 2L * pad)
  cl  <- cumsum(brk)

  out <- lapply(split(pos, cl), function(p) split_cluster(p, pad, max_width))
  res <- do.call(rbind, out)
  res <- res[res$n_cpg >= min_cpgs, , drop = FALSE]
  res <- res[order(res$start), , drop = FALSE]
  rownames(res) <- NULL
  res
}

## One cluster of CpGs -> one or more regions, each at most max_width wide.
##
## The cut points are chosen between CpGs, not at arbitrary coordinates, so no
## CpG is orphaned and no region is padded past its own outermost CpG.  A cluster
## that already fits is returned unchanged.
split_cluster <- function(p, pad, max_width) {
  span <- function(v) (max(v) + pad) - (min(v) - pad) + 1L
  if (span(p) <= max_width || length(p) == 1L)
    return(data.frame(start = min(p) - pad, end = max(p) + pad, n_cpg = length(p)))

  ## how many pieces the cluster needs, then cut it into that many contiguous
  ## runs of CpGs of near-equal count
  k <- as.integer(ceiling(span(p) / max_width))
  repeat {
    grp <- cut(seq_along(p), breaks = k, labels = FALSE)
    w <- vapply(split(p, grp), span, 0L)
    if (all(w <= max_width) || k >= length(p)) break
    k <- k + 1L
  }
  parts <- unname(split(p, grp))
  lo <- vapply(parts, min, 0L); hi <- vapply(parts, max, 0L)
  ## Pad outwards, but put the boundary between two pieces at the midpoint of the
  ## gap that separates them, so the pieces abut instead of overlapping.
  starts <- lo - pad; ends <- hi + pad
  if (length(parts) > 1L) {
    mid <- (hi[-length(hi)] + lo[-1]) %/% 2L
    ends[-length(ends)]  <- pmin(ends[-length(ends)],  mid)
    starts[-1]           <- pmax(starts[-1],           mid + 1L)
  }
  data.frame(start = starts, end = ends, n_cpg = vapply(parts, length, 0L))
}

## Build regions for every chromosome in a CpG table.
##
##   cpg  data.frame with `chr` and `pos`
## Returns chr, start, end, n_cpg, width, interval_id (unique, ordered).
build_regions <- function(cpg, pad = 250L, max_width = 1000L, min_cpgs = 2L) {
  stopifnot(all(c("chr", "pos") %in% names(cpg)))
  chrs <- unique(cpg$chr)
  ## order chromosomes naturally (chr1 .. chr22) rather than lexically
  num <- suppressWarnings(as.integer(sub("^chr", "", chrs)))
  chrs <- chrs[order(is.na(num), num, chrs)]

  out <- lapply(chrs, function(ch) {
    r <- merge_positions(cpg$pos[cpg$chr == ch], pad, max_width, min_cpgs)
    if (!nrow(r)) return(NULL)
    cbind(data.frame(chr = ch, stringsAsFactors = FALSE), r)
  })
  res <- do.call(rbind, Filter(Negate(is.null), out))
  if (is.null(res) || !nrow(res))
    stop("no regions survived: check pad / max_width / min_cpgs against the input")
  res$width <- res$end - res$start + 1L
  res$interval_id <- seq_len(nrow(res))
  rownames(res) <- NULL
  res
}

## Check a region table is well formed: sorted, non-overlapping within a
## chromosome, inside the width cap, and carrying at least min_cpgs.
validate_regions <- function(reg, max_width = 1000L, min_cpgs = 2L) {
  stopifnot(all(c("chr", "start", "end", "n_cpg", "width", "interval_id") %in% names(reg)))
  if (any(reg$end < reg$start)) stop("some regions end before they start")
  if (any(reg$width > max_width))
    stop(sum(reg$width > max_width), " regions exceed max_width ", max_width)
  if (any(reg$n_cpg < min_cpgs))
    stop(sum(reg$n_cpg < min_cpgs), " regions hold fewer than ", min_cpgs, " CpGs")
  if (anyDuplicated(reg$interval_id)) stop("interval_id is not unique")
  for (ch in unique(reg$chr)) {
    r <- reg[reg$chr == ch, , drop = FALSE]
    r <- r[order(r$start), , drop = FALSE]
    if (nrow(r) > 1L && any(r$start[-1] <= r$end[-nrow(r)]))
      stop("regions overlap on ", ch)
  }
  invisible(TRUE)
}

## Read the CpG list that defines the regions.
##
## Two layouts are accepted:
##   * a table with `chr` and `pos` columns (1-based position of the CpG's C);
##   * a table with a coordinate column of the form "chrN:pos" or
##     "chrN:start-end", named by `col` (for example "hg38").
## Comma- or tab-separated, plain or gzipped.  A `Category` column, if present,
## is carried along and can be used to keep only some rows.  Files written on
## Windows (CRLF) are handled: every field is trimmed.
read_cpg_list <- function(path, col = "hg38", categories = NULL) {
  d <- as.data.frame(data.table::fread(path, colClasses = "character", showProgress = FALSE))
  names(d) <- trimws(names(d))
  for (j in seq_along(d)) d[[j]] <- trimws(d[[j]])
  if (!is.null(categories)) {
    if (!"Category" %in% names(d)) stop("no Category column to filter on in ", path)
    keep <- d$Category %in% categories
    if (!any(keep)) stop("no rows match categories: ", paste(categories, collapse = ", "))
    d <- d[keep, , drop = FALSE]
  }
  if (all(c("chr", "pos") %in% names(d))) {
    chr <- d$chr; pos <- suppressWarnings(as.integer(d$pos)); loc <- paste0(chr, ":", pos)
  } else {
    if (!col %in% names(d))
      stop("need either chr/pos columns or a coordinate column '", col, "' in ", path,
           " (have: ", paste(names(d), collapse = ", "), ")")
    loc <- d[[col]]
    chr <- sub(":.*$", "", loc)
    rest <- sub("^[^:]*:", "", loc)
    pos <- suppressWarnings(as.integer(sub("-.*$", "", rest)))
    end <- suppressWarnings(as.integer(sub("^.*-", "", rest)))
    if (any(!is.na(pos) & !is.na(end) & pos != end))
      warning(sum(pos != end, na.rm = TRUE), " entries are not single-base; using the start")
  }
  if (anyNA(pos)) stop(sum(is.na(pos)), " positions could not be parsed in ", path)
  data.frame(chr = chr, pos = pos, id = loc,
             category = if ("Category" %in% names(d)) d$Category else NA_character_,
             stringsAsFactors = FALSE)
}

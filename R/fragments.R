## Fragments overlapped onto regions, with the reference counts they need.
##
## A .pat row is one observed read pattern: a start CpG index, a string of C / T
## / . over consecutive CpGs, and how many identical reads it stands for.  What
## the likelihood needs per fragment is the part that falls inside a region, the
## base-pair gaps between those CpGs, and the panel's M / N at each of them.
##
## Two things here are deliberate:
##
##   * Identical (fragment, gaps, M, N) tuples are collapsed and their counts
##     summed.  The log-likelihood depends on nothing else, so duplicates are
##     pure waste -- and reads are highly repetitive, especially short ones.
##   * The result is aggregated to region level before anything is written.  The
##     original kept per-fragment likelihoods for all 21 patterns on disk, which
##     is where its 238 GB of intermediates came from.


## Read a .beta file as an n x 2 integer matrix of (methylated, total).
read_beta_counts <- function(path, n_cpg = NULL) {
  sz <- file.size(path)
  if (is.na(sz)) stop("cannot stat ", path)
  if (sz %% 2L != 0L) stop(path, " has an odd byte count; not a .beta file")
  n <- as.integer(sz %/% 2L)
  if (!is.null(n_cpg) && n_cpg != n)
    stop(path, " holds ", n, " CpGs but ", n_cpg, " were expected")
  raw <- readBin(path, "integer", n = 2L * n, size = 1L, signed = FALSE)
  matrix(raw, ncol = 2L, byrow = TRUE)
}

## Read a .pat.gz, keeping only the columns that matter.
## Columns are chr, startCpG, pattern, count.
read_pat <- function(path, chrom = NULL) {
  ## gzipped files are streamed through gzip so that R.utils is not needed
  d <- if (grepl("[.]gz$", path)) {
    data.table::fread(cmd = paste("gzip -dc", shQuote(path)), sep = "\t", header = FALSE,
                      showProgress = FALSE, select = 1:4,
                      col.names = c("chr", "startCpG", "fragment", "count"))
  } else {
    data.table::fread(path, sep = "\t", header = FALSE, showProgress = FALSE,
                      select = 1:4, col.names = c("chr", "startCpG", "fragment", "count"))
  }
  if (!is.null(chrom)) d <- d[chr %in% chrom]
  d[, startCpG := as.integer(startCpG)]
  d[, count := as.integer(count)]
  d[, endCpG := startCpG + nchar(fragment) - 1L]
  d[]
}

## Overlap fragments onto regions and trim each to the overlapping part.
##
##   pat   from read_pat()
##   reg   an indexed region table (startCpG / endCpG / cpg_gaps / interval_id)
##   beta  the n x 2 count matrix
##
## Returns one row per distinct trimmed fragment within a region.
fragments_in_regions <- function(pat, reg, beta, min_sites = 1L) {
  R <- data.table::as.data.table(reg[, c("chr", "interval_id", "startCpG", "endCpG", "cpg_gaps")])
  data.table::setnames(R, c("startCpG", "endCpG"), c("regStart", "regEnd"))
  P <- data.table::as.data.table(pat)

  data.table::setkey(R, regStart, regEnd)
  ov <- data.table::foverlaps(P, R, by.x = c("startCpG", "endCpG"), nomatch = 0L)
  if (!nrow(ov)) return(ov[0])

  ## the overlapping index range, and where it sits inside the read
  ov[, ov_start := pmax(startCpG, regStart)]
  ov[, ov_end   := pmin(endCpG,   regEnd)]
  ov <- ov[ov_end >= ov_start]
  ov[, trimmed := substr(fragment, ov_start - startCpG + 1L, ov_end - startCpG + 1L)]
  ov[, n_obs := nchar(gsub(".", "", trimmed, fixed = TRUE))]
  ov <- ov[n_obs >= min_sites]
  if (!nrow(ov)) return(ov[0])

  ## collapse identical observations -- the likelihood cannot tell them apart
  fr <- ov[, .(count = sum(count)),
           by = .(chr, interval_id, ov_start, ov_end, trimmed, regStart, cpg_gaps)]

  ## gaps for the trimmed span: the region's gap vector, sliced
  fr[, dist_sub := slice_gaps(cpg_gaps, ov_start - regStart + 1L, ov_end - regStart + 1L)]

  ## reference counts over the same span
  Mv <- as.integer(beta[, 1]); Nv <- as.integer(beta[, 2])
  fr[, M := vapply(seq_len(.N), function(i) paste(Mv[ov_start[i]:ov_end[i]], collapse = ","), "")]
  fr[, N := vapply(seq_len(.N), function(i) paste(Nv[ov_start[i]:ov_end[i]], collapse = ","), "")]

  fr[, c("cpg_gaps", "regStart") := NULL]
  data.table::setnames(fr, "trimmed", "fragment")
  fr[]
}

## Slice a region's comma-separated gap vector down to positions s..e (1-based
## within the region).  A span of k CpGs needs k-1 gaps: those between its own
## consecutive members, i.e. elements s..(e-1) of the region's vector.
slice_gaps <- function(gaps, s, e) {
  vapply(seq_along(gaps), function(i) {
    if (e[i] <= s[i]) return("")
    g <- strsplit(gaps[i], ",", fixed = TRUE)[[1]]
    if (!length(g)) return("")
    lo <- s[i]; hi <- min(e[i] - 1L, length(g))
    if (lo > hi) return("")
    paste(g[lo:hi], collapse = ",")
  }, "")
}

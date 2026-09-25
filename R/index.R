## Attach CpG indices to base-pair regions.
##
## A .pat fragment names its CpGs by row number in the genome-wide CpG list, and
## a .beta file is that same list, so every region has to carry the index range
## it spans as well as its coordinates.  `cpg_gaps` is the base-pair distance
## between consecutive CpGs inside the region, which is what the correlation
## kernel rho(d) = phi * exp(-d/s) is evaluated on.
##
## The reference is ALL_chr_pos.RDS: one row per CpG, in index order, with chr
## and pos.  It is read once and reused for every chromosome.

## Load the genome-wide CpG coordinate table and add the index column.
##
## Row k of the table must be row k of every .beta file, so the table is used in
## stored order and never sorted.  Accepted: an .rds holding a data frame, or a
## text / BED file (optionally gzipped); in either case the first two columns are
## taken as chr and pos unless columns of those names exist.  If the table carries
## its own index column (`idx`, or the third column of a wgbstools CpG.bed), it
## has to equal the row number.
load_cpg_index <- function(path) {
  ref <- if (grepl("[.]rds$", path, ignore.case = TRUE)) readRDS(path) else
    data.table::fread(path, header = "auto", showProgress = FALSE)
  ref <- as.data.frame(ref)
  if (!all(c("chr", "pos") %in% names(ref))) names(ref)[1:2] <- c("chr", "pos")
  own <- if ("idx" %in% names(ref)) ref$idx else if (ncol(ref) >= 3 && is.numeric(ref[[3]])) ref[[3]] else NULL
  if (!is.null(own) && !identical(as.integer(own), seq_len(nrow(ref))))
    stop(path, ": its index column is not the row number, so it is not in .beta order")
  ref <- data.frame(chr = as.character(ref$chr), pos = as.integer(ref$pos),
                    stringsAsFactors = FALSE)
  ref$idx <- seq_len(nrow(ref))
  ref
}

## Add startCpG, endCpG, n_cpg_in_region and cpg_gaps to a region table.
##
## Regions with no CpG in the reference are dropped, with a count reported: that
## happens when the region list and the reference are on different builds, and it
## is better seen than silently carried.
index_regions <- function(reg, ref, verbose = TRUE) {
  stopifnot(all(c("chr", "start", "end") %in% names(reg)))
  out <- list()
  for (ch in unique(reg$chr)) {
    r <- reg[reg$chr == ch, , drop = FALSE]
    p <- ref$pos[ref$chr == ch]
    i0 <- ref$idx[ref$chr == ch]
    if (!length(p)) {
      if (verbose) message("  ", ch, ": absent from the reference, dropped")
      next
    }
    ord <- order(p); p <- p[ord]; i0 <- i0[ord]

    ## first CpG at or after start, last CpG at or before end
    lo <- findInterval(r$start - 1L, p) + 1L
    hi <- findInterval(r$end, p)
    keep <- lo <= hi
    if (verbose && any(!keep))
      message("  ", ch, ": ", sum(!keep), " of ", nrow(r), " regions contain no reference CpG")
    r <- r[keep, , drop = FALSE]; lo <- lo[keep]; hi <- hi[keep]
    if (!nrow(r)) next

    r$startCpG <- i0[lo]
    r$endCpG   <- i0[hi]
    r$n_cpg_in_region <- hi - lo + 1L
    ## gaps between consecutive CpGs, as the original stored them: one comma-
    ## separated string per region, length (n_cpg - 1)
    r$cpg_gaps <- vapply(seq_len(nrow(r)), function(k) {
      if (hi[k] <= lo[k]) return("")
      paste(diff(p[lo[k]:hi[k]]), collapse = ",")
    }, "")
    out[[length(out) + 1L]] <- r
  }
  res <- do.call(rbind, out)
  if (is.null(res) || !nrow(res)) stop("no region could be indexed against the reference")
  res <- res[order(res$startCpG), , drop = FALSE]
  res$interval_id <- seq_len(nrow(res))
  rownames(res) <- NULL
  res
}

## Consistency checks that would otherwise surface much later, as silently wrong
## likelihoods rather than as an error.
validate_indexed <- function(reg) {
  if (any(reg$endCpG < reg$startCpG)) stop("endCpG before startCpG")
  n <- reg$endCpG - reg$startCpG + 1L
  if (!identical(n, reg$n_cpg_in_region)) stop("n_cpg_in_region disagrees with the index range")
  g <- lengths(strsplit(reg$cpg_gaps, ",", fixed = TRUE))
  g[reg$cpg_gaps == ""] <- 0L
  if (any(g != n - 1L))
    stop(sum(g != n - 1L), " regions have the wrong number of gaps for their CpG count")
  ## index ranges must not overlap, or a fragment would be claimed by two regions
  o <- order(reg$startCpG)
  s <- reg$startCpG[o]; e <- reg$endCpG[o]
  if (length(s) > 1L && any(s[-1] <= e[-length(e)]))
    stop("CpG index ranges overlap between regions")
  invisible(TRUE)
}

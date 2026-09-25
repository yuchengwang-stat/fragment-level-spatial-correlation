## The per-region result tables.
##
##   reg   the indexed region table from step 1
##   cpg   the CpG list the regions were built from (read_cpg_list()), with the
##         same filtering
##   post  a named list, one element per cell type, each the `posterior` table
##         written by step 3 (names become the cell-type labels)
##   het   optional, the output of heterogeneity_stats() for the same regions and
##         cell types; its uniformity columns are added to `long`
##
## Returns list(long, wide).
##
##   long  one row per region x cell type:
##           region_id, region, chr, n_list_cpg, list_cpgs, celltype,
##           phi, s, meth, post_winner, phi_mean, s_mean
##   wide  one row per region: coordinates, the list CpGs, and per cell type
##           <ct>_phi, <ct>_s, <ct>_meth
##
## phi and s are the pattern with the largest posterior probability; post_winner
## is that probability; phi_mean and s_mean are posterior expectations over all
## patterns.  meth is the region's methylation in that cell type, the quantity
## that chose its prior bin.  A cell type with no fragment in a region has no row
## in `long` and blank cells in `wide`.
region_tables <- function(reg, cpg, post, het = NULL) {
  if (is.null(names(post)) || any(!nzchar(names(post))))
    stop("post must be a named list: one posterior table per cell type")
  reg <- data.table::as.data.table(reg)
  reg <- reg[, list(chr, interval_id, start, end, width,
                    n_cpg_region = n_cpg_in_region, n_list_expected = n_cpg)]

  ## place the list CpGs in regions (regions do not overlap within a chromosome)
  cp <- unique(data.table::as.data.table(cpg)[, list(chr, pos, category)])
  cp[, interval_id := NA_integer_]
  for (ch in unique(reg$chr)) {
    r <- reg[chr == ch][order(start)]
    w <- which(cp$chr == ch)
    if (!length(w)) next
    k <- findInterval(cp$pos[w], r$start)
    hit <- k > 0L
    hit[hit] <- cp$pos[w][hit] <= r$end[k[hit]]
    cp[w[hit], interval_id := r$interval_id[k[hit]]]
  }
  data.table::setorder(cp, interval_id, pos)
  cl <- cp[!is.na(interval_id), list(n_list_cpg = .N,
                                     list_cpgs = paste0(chr, ":", pos, collapse = ";"),
                                     list_category = paste(category, collapse = ";")),
           by = interval_id]
  reg <- merge(reg, cl, by = "interval_id", all.x = TRUE)
  bad <- reg[is.na(n_list_cpg) | n_list_cpg != n_list_expected]
  if (nrow(bad))
    stop(nrow(bad), " regions hold a different number of list CpGs than they were built",
         " from -- pass the same CpG list, column and category filter as step 1")
  reg[, n_list_expected := NULL]
  if (all(is.na(cp$category))) reg[, list_category := NULL]
  reg[, region := sprintf("%s:%d-%d", chr, start, end)]

  long <- data.table::rbindlist(lapply(names(post), function(ct) {
    d <- data.table::as.data.table(post[[ct]])
    pat <- fc_patterns()
    pm <- as.matrix(d[, paste0(pat$name, "_post"), with = FALSE])
    d[, list(celltype = ct, chr, interval_id, meth, phi = best_phi, s = best_s,
             post_winner = best_post,
             phi_mean = as.vector(pm %*% pat$phi), s_mean = as.vector(pm %*% pat$s))]
  }))
  cts <- sort(unique(long$celltype))

  wide <- data.table::dcast(long, chr + interval_id ~ celltype, value.var = c("phi", "s", "meth"))
  old <- unlist(lapply(cts, function(ct) paste0(c("phi_", "s_", "meth_"), ct)))
  new <- unlist(lapply(cts, function(ct) paste0(ct, c("_phi", "_s", "_meth"))))
  data.table::setnames(wide, old, new)
  for (j in grep("_meth$", new, value = TRUE)) data.table::set(wide, j = j, value = round(wide[[j]], 4))
  lead <- intersect(c("region", "chr", "start", "end", "width", "n_cpg_region", "n_list_cpg",
                      "list_cpgs", "list_category"), names(reg))
  wide <- merge(reg, wide, by = c("chr", "interval_id"), all.x = TRUE)
  data.table::setorder(wide, interval_id)
  data.table::setnames(wide, "interval_id", "region_id")
  data.table::setcolorder(wide, c("region_id", lead, new))

  long <- merge(reg[, list(chr, interval_id, region, n_list_cpg, list_cpgs)], long,
                by = c("chr", "interval_id"))
  data.table::setorder(long, interval_id, celltype)
  data.table::setnames(long, "interval_id", "region_id")
  long[, `:=`(meth = round(meth, 4), post_winner = round(post_winner, 4),
              phi_mean = round(phi_mean, 4), s_mean = round(s_mean, 1))]
  data.table::setcolorder(long, c("region_id", "region", "chr", "n_list_cpg", "list_cpgs",
                                  "celltype", "phi", "s", "meth", "post_winner",
                                  "phi_mean", "s_mean"))
  if (!is.null(het)) {
    h <- data.table::as.data.table(het)[, list(region_id, celltype, meth_sd = round(meth_sd, 4),
                                               meth_jump = round(meth_jump, 4),
                                               r2_split = round(r2_split, 4), target_lowhigh, layout)]
    long <- merge(long, h, by = c("region_id", "celltype"), all.x = TRUE)
    data.table::setorder(long, region_id, celltype)
    data.table::setcolorder(long, c("region_id", "region", "chr", "n_list_cpg", "list_cpgs",
                                    "celltype"))
  }
  list(long = long[], wide = wide[])
}

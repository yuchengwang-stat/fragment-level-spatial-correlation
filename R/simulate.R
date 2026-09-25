## A small synthetic dataset in the file formats the pipeline reads.
##
## It exists so that the whole pipeline can be run end to end in about a minute,
## and so that its output can be checked against a known answer.  Reads are drawn
## from the model the method assumes: within a read, the methylation states of
## the CpGs are thresholded correlated Gaussians with
##
##   rho(d) = phi * exp(-d / s),
##
## and every cell type has its own (phi, s).  The .beta counts are tallied from
## the same reads, as wgbstools does, so the marginals the method plugs in are
## the ones a real .beta would give.
##
## Written into `dir`:
##   cpg_index.rds     chr, pos -- row k is row k of every .beta
##   cpg_list.csv      the CpGs that define regions, as an "hg38" coordinate
##                     column ("chr:pos-pos") and a Category column
##   <sample>.pat.gz   chr, startCpG, pattern, count; sorted by startCpG
##   <sample>.beta     uint8 (M, N) pairs, one per CpG of the index
##   samples.tsv       sample, celltype, pat, beta
## Returns those paths and the truth, invisibly.
simulate_dataset <- function(dir, n_chr = 2L, cpgs_per_chr = 3000L, n_loci = 100L,
                             celltypes = list(A = c(phi = 0.8, s = 100),
                                              B = c(phi = 0.3, s = 300)),
                             samples_per_type = 3L, reads_per_locus = 80L,
                             read_cpgs = c(3L, 9L), background_reads = 3000L,
                             seed = 1L) {
  set.seed(seed)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  ## ---- a genome of CpGs: clustered spacing, like CpG islands in sparse sea ----
  chrs <- paste0("chr", seq_len(n_chr))
  gaps <- function(n) ifelse(stats::runif(n) < 0.5, sample(2:40, n, TRUE), sample(40:400, n, TRUE))
  idx <- data.frame(chr = rep(chrs, each = cpgs_per_chr),
                    pos = unlist(lapply(chrs, function(ch) 10000L + cumsum(gaps(cpgs_per_chr)))),
                    stringsAsFactors = FALSE)
  n <- nrow(idx)
  chr_lo <- match(chrs, idx$chr); chr_hi <- c(chr_lo[-1] - 1L, n)

  ## ---- methylation: every locus sits at one level, sites jitter around it ----
  p <- stats::rbeta(n, 0.5, 0.5)
  ## loci are spaced along the genome, away from chromosome ends
  per_chr <- ceiling(n_loci / n_chr)
  anchors <- unlist(lapply(seq_along(chrs), function(k)
    round(seq(chr_lo[k] + 30L, chr_hi[k] - 60L, length.out = per_chr))))[seq_len(n_loci)]
  levels <- c(0.03, 0.2, 0.5, 0.8, 0.97)
  loci <- lapply(seq_along(anchors), function(i) {
    a <- anchors[i]
    k <- sample(2:5, 1L)                         # list CpGs at this locus
    list_idx <- a + seq_len(k) - 1L
    ## the reads cover the list CpGs and their neighbours within 400 bp
    ch <- idx$chr[a]
    win <- which(idx$chr == ch & idx$pos >= idx$pos[a] - 400L &
                 idx$pos <= idx$pos[a + k - 1L] + 400L)
    list(list_idx = list_idx, win = win, level = sample(levels, 1L))
  })
  for (l in loci) p[l$win] <- stats::plogis(stats::qlogis(l$level) + stats::rnorm(length(l$win), 0, 0.3))

  cat_draw <- sample(c("HighValue", "EPIC"), length(loci), TRUE, prob = c(0.8, 0.2))
  li <- unlist(lapply(loci, `[[`, "list_idx"))
  lc <- rep(cat_draw, vapply(loci, function(l) length(l$list_idx), 0L))
  loc <- paste0(idx$chr[li], ":", idx$pos[li], "-", idx$pos[li])
  utils::write.csv(data.frame(hg38 = loc, Category = lc), file.path(dir, "cpg_list.csv"),
                   row.names = FALSE, quote = FALSE)
  saveRDS(idx, file.path(dir, "cpg_index.rds"))

  ## ---- reads -------------------------------------------------------------------
  one_read <- function(st, len, phi, s, correlated) {
    ch <- idx$chr[st]
    en <- min(st + len - 1L, chr_hi[match(ch, chrs)])
    ii <- st:en
    if (correlated && length(ii) > 1L) {
      Sig <- phi * exp(-abs(outer(idx$pos[ii], idx$pos[ii], "-")) / s); diag(Sig) <- 1
      z <- as.vector(stats::rnorm(length(ii)) %*% chol(Sig))
    } else {
      z <- stats::rnorm(length(ii))
    }
    ct <- ifelse(z <= stats::qnorm(p[ii]), "C", "T")
    ## a few uncalled sites inside the read, as real reads have
    if (length(ct) > 2L) ct[-1L][stats::runif(length(ct) - 1L) < 0.05] <- "."
    list(start = st, pattern = paste(ct, collapse = ""))
  }

  sheet <- list()
  for (ct_name in names(celltypes)) {
    phi <- celltypes[[ct_name]][["phi"]]; s <- celltypes[[ct_name]][["s"]]
    for (r in seq_len(samples_per_type)) {
      smp <- sprintf("%s_%d", ct_name, r)
      reads <- list()
      for (l in loci) {
        st <- sample(l$win, reads_per_locus, TRUE)
        len <- sample(read_cpgs[1]:read_cpgs[2], reads_per_locus, TRUE)
        reads[[length(reads) + 1L]] <- lapply(seq_len(reads_per_locus), function(j)
          one_read(st[j], len[j], phi, s, TRUE))
      }
      bst <- sample.int(n, background_reads, TRUE)
      blen <- sample(read_cpgs[1]:read_cpgs[2], background_reads, TRUE)
      reads[[length(reads) + 1L]] <- lapply(seq_len(background_reads), function(j)
        one_read(bst[j], blen[j], 0, 1, FALSE))
      rd <- data.table::rbindlist(unlist(reads, recursive = FALSE))

      ## .pat: identical reads collapsed into a count, sorted by start CpG
      pt <- rd[, list(count = .N), by = list(start, pattern)]
      data.table::setorder(pt, start, pattern)
      pt[, chr := idx$chr[start]]
      f_pat <- file.path(dir, paste0(smp, ".pat.gz"))
      con <- gzfile(f_pat, "w")
      writeLines(paste(pt$chr, pt$start, pt$pattern, pt$count, sep = "\t"), con)
      close(con)

      ## .beta: tallied from the same reads.  The format holds uint8, so where
      ## N > 255 both counts are scaled down to N = 255 with the ratio kept, as
      ## wgbstools does; capping M and N separately would bias M / N.
      ch <- strsplit(rd$pattern, "", fixed = TRUE)
      site <- rep(rd$start, lengths(ch)) + sequence(lengths(ch)) - 1L
      call <- unlist(ch)
      M <- tabulate(site[call == "C"], n); N <- tabulate(site[call %in% c("C", "T")], n)
      big <- N > 255L
      M[big] <- round(M[big] / N[big] * 255); N[big] <- 255L
      f_beta <- file.path(dir, paste0(smp, ".beta"))
      writeBin(as.integer(as.vector(rbind(M, N))), f_beta, size = 1L)

      sheet[[length(sheet) + 1L]] <- data.frame(sample = smp, celltype = ct_name,
                                                pat = f_pat, beta = f_beta)
    }
  }
  sheet <- do.call(rbind, sheet)
  utils::write.table(sheet, file.path(dir, "samples.tsv"), sep = "\t",
                     quote = FALSE, row.names = FALSE)
  truth <- data.frame(celltype = names(celltypes),
                      phi = vapply(celltypes, `[[`, 0, "phi"),
                      s = vapply(celltypes, `[[`, 0, "s"), row.names = NULL)
  invisible(list(dir = dir, cpg_index = file.path(dir, "cpg_index.rds"),
                 cpg_list = file.path(dir, "cpg_list.csv"),
                 samples = file.path(dir, "samples.tsv"), truth = truth))
}

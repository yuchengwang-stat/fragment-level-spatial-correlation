test_that("the pipeline separates a strongly from a weakly correlated cell type", {
  dir <- file.path(tempdir(), "fragcorr_test")
  d <- simulate_dataset(dir, n_loci = 60L, samples_per_type = 2L, seed = 3L)
  cpg <- read_cpg_list(d$cpg_list, col = "hg38")
  reg <- index_regions(build_regions(cpg), load_cpg_index(d$cpg_index), verbose = FALSE)
  expect_silent(validate_indexed(reg))
  sh <- utils::read.delim(d$samples)
  pn <- fc_patterns()$name
  post <- lapply(split(sh, sh$celltype), function(x) {
    liks <- lapply(seq_len(nrow(x)), function(i)
      sample_loglik(reg, read_pat(x$pat[i]), read_beta_counts(x$beta[i])))
    L <- combine_loglik(liks, pn)
    m <- data.table::rbindlist(lapply(liks, function(l) l[, list(chr, interval_id, M_sum, N_sum)]))[
      , list(M_sum = sum(M_sum), N_sum = sum(N_sum)), by = list(chr, interval_id)]
    L <- merge(L, m, by = c("chr", "interval_id"))
    L[, meth := M_sum / N_sum][, bin := assign_bin(meth)]
    pri <- fit_priors(L, pn)
    po <- posterior_from_prior(L, pri$pi, pn)
    expect_equal(unname(rowSums(po)), rep(1, nrow(po)), tolerance = 1e-8)
    cbind(L[, list(chr, interval_id, meth, bin)], summarise_posterior(po, fc_patterns()),
          data.table::as.data.table(po))
  })
  tb <- region_tables(reg, cpg, post)
  expect_equal(nrow(tb$wide), nrow(reg))
  phi <- tapply(tb$long$phi_mean, tb$long$celltype, mean)
  expect_gt(phi[["A"]], 0.6)          # truth 0.8
  expect_lt(phi[["B"]], 0.5)          # truth 0.3
  cv <- decay_curves(post, boot = 50L, min_regions = 5L)
  expect_true(all(cv$band$lo <= cv$band$mean & cv$band$mean <= cv$band$hi))
})

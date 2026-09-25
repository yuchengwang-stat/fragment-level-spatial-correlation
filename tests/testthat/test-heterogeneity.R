test_that("uniform, graded, split and interleaved regions are told apart", {
  N <- 100L
  p <- list(uniform     = rep(0.5, 8),
            graded      = seq(0.90, 0.40, length.out = 8),   # sd 0.19: not uniform, no extremes
            split       = c(rep(0.05, 4), rep(0.90, 4)),
            interleaved = rep(c(0.05, 0.90), 4))
  M <- matrix(round(unlist(p) * N), ncol = 1, dimnames = list(NULL, "g"))
  Nm <- matrix(N, nrow(M), 1, dimnames = list(NULL, "g"))
  rid <- rep(seq_along(p), each = 8)
  h <- heterogeneity_stats(M, Nm, rid, target = rep(c(TRUE, rep(FALSE, 6), TRUE), 4))
  expect_equal(h$layout, c("uniform", "graded", "mixed_split", "mixed_interleaved"))
  expect_equal(h$r2_split[3], 1)
  expect_lt(h$r2_split[4], 0.2)
  expect_equal(h$target_lowhigh, c(FALSE, FALSE, TRUE, TRUE))   # targets are CpGs 1 and 8
})

test_that("CpGs below the depth threshold are ignored", {
  M <- matrix(c(10L, 10L, 1L, 10L), ncol = 1, dimnames = list(NULL, "g"))
  N <- matrix(c(20L, 20L, 2L, 20L), ncol = 1, dimnames = list(NULL, "g"))
  h <- heterogeneity_stats(M, N, rep(1L, 4), min_n = 20L)
  expect_equal(h$n_cov, 3L)
  expect_equal(h$meth_all_sd, 0)
  expect_equal(h$meth_all_iqr, 0)
  expect_equal(h$meth_all, 0.5)
})

test_that("means and variances are reported for all CpGs and for targets", {
  M <- matrix(c(10L, 90L, 50L, 50L), ncol = 1, dimnames = list(NULL, "g"))
  N <- matrix(100L, 4, 1, dimnames = list(NULL, "g"))
  h <- heterogeneity_stats(M, N, rep(1L, 4), target = c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(h$meth_all, 0.5)
  expect_equal(h$meth_target, 0.5)
  expect_equal(h$meth_target_sd, stats::sd(c(0.1, 0.9)))
  expect_equal(h$meth_all_iqr, stats::IQR(c(0.1, 0.9, 0.5, 0.5)))
  expect_true(h$target_lowhigh)
})

test_that("a region whose targets are low and high in enough groups is cut between them", {
  ## one region, four targets at 100, 150, 600, 650; groups g1..g6
  reg <- data.frame(chr = "chr1", start = 1L, end = 900L, n_cpg = 4L, width = 900L,
                    interval_id = 1L, startCpG = 1L, endCpG = 4L, n_cpg_in_region = 4L)
  targets <- data.table::data.table(chr = "chr1", pos = c(100L, 150L, 600L, 650L),
                                    row = 1:4, region_id = 1L)
  P <- matrix(c(0.05, 0.1, 0.9, 0.95), 4, 6, dimnames = list(NULL, paste0("g", 1:6)))
  r5 <- refine_regions(reg, targets, P, min_celltypes = 5L, pad = 250L)
  expect_equal(nrow(r5), 2L)
  expect_equal(r5$n_cpg, c(2L, 2L))
  expect_equal(r5$end[1] + 1L, r5$start[2])              # the pieces abut
  expect_equal(r5$end[1], (150L + 600L) %/% 2L)          # at the midpoint between them
  expect_true(all(r5$parent_region == 1L))
  ## the same region with the conflict in only two groups is left whole
  P2 <- P; P2[, 3:6] <- 0.5
  r2 <- refine_regions(reg, targets, P2, min_celltypes = 5L)
  expect_equal(nrow(r2), 1L)
  expect_equal(c(r2$start, r2$end), c(1L, 900L))
})

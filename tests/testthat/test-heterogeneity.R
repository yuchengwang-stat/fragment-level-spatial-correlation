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
  expect_equal(h$meth_sd, 0)
})

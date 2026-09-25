test_that("nearby CpGs merge into one region, isolated ones are dropped", {
  r <- merge_positions(c(1000L, 1100L, 1200L, 5000L), pad = 250L, max_width = 1000L, min_cpgs = 2L)
  expect_equal(nrow(r), 1L)
  expect_equal(r$start, 750L)
  expect_equal(r$end, 1450L)
  expect_equal(r$n_cpg, 3L)
})

test_that("a dense cluster is cut into pieces within the width cap, without overlap", {
  pos <- seq(1000L, 5000L, by = 50L)
  r <- merge_positions(pos, pad = 250L, max_width = 1000L, min_cpgs = 2L)
  expect_gt(nrow(r), 1L)
  expect_true(all(r$end - r$start + 1L <= 1000L))
  expect_true(all(r$start[-1] > r$end[-nrow(r)]))
  expect_equal(sum(r$n_cpg), length(pos))
})

test_that("regions are indexed against the CpG table in .beta order", {
  ref <- data.frame(chr = rep(c("chr1", "chr2"), each = 5),
                    pos = c(100L, 150L, 200L, 900L, 950L, 100L, 130L, 160L, 190L, 220L),
                    idx = 1:10)
  cpg <- data.frame(chr = c("chr1", "chr1", "chr2", "chr2"), pos = c(150L, 200L, 130L, 160L))
  reg <- build_regions(cpg, pad = 60L, max_width = 1000L, min_cpgs = 2L)
  ix <- index_regions(reg, ref, verbose = FALSE)
  expect_equal(ix$startCpG, c(1L, 6L))
  expect_equal(ix$endCpG, c(3L, 10L))          # chr2 region is 70-220 bp: CpGs 100..220
  expect_equal(ix$cpg_gaps, c("50,50", "30,30,30,30"))
  expect_silent(validate_indexed(ix))
})

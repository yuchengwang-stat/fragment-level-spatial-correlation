pat <- fc_patterns()

test_that("a one-site fragment scores the same under every pattern", {
  L <- loglik_all_patterns("C", "", "30", "100", pat$phi, pat$s, print_every = 0L)
  expect_equal(length(unique(as.vector(L))), 1L)
  expect_equal(L[1, 1], round(log(0.3), 3))
})

test_that("sites below the depth threshold are ignored", {
  L <- loglik_all_patterns("CT", "50", "1,2", "3,4", pat$phi, pat$s, min_n = 5L, print_every = 0L)
  expect_true(all(is.na(L)))
})

test_that("concordant fragments favour strong correlation, discordant ones weak", {
  frag <- c("CCCC", "CTCT")
  L <- loglik_all_patterns(frag, rep("30,30,30", 2), rep("50,50,50,50", 2),
                           rep("100,100,100,100", 2), pat$phi, pat$s, print_every = 0L)
  strong <- which(pat$phi == 1 & pat$s == 300); weak <- which(pat$phi == 0.1 & pat$s == 30)
  expect_gt(L[1, strong], L[1, weak])
  expect_lt(L[2, strong], L[2, weak])
})

test_that("the closed-form two-site shortcut stays close to EP and keeps NA the same", {
  set.seed(1)
  frag <- replicate(50, paste(sample(c("C", "T"), 2, TRUE), collapse = ""))
  d <- as.character(sample(5:500, 50, TRUE))
  M <- replicate(50, paste(sample(5:95, 2), collapse = ","))
  N <- rep("100,100", 50)
  ep  <- loglik_all_patterns(frag, d, M, N, pat$phi, pat$s, print_every = 0L)
  bvn <- loglik_all_patterns(frag, d, M, N, pat$phi, pat$s, use_bvn = TRUE, print_every = 0L)
  expect_identical(is.na(ep), is.na(bvn))
  expect_lt(max(abs(exp(ep) - exp(bvn)), na.rm = TRUE), 0.05)
})

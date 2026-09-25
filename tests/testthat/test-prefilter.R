test_that("the .pat prefilter reads to the end, so a pipefail pipeline succeeds", {
  skip_on_os("windows")
  if (!nzchar(Sys.which("bash")) || !nzchar(Sys.which("gzip")) || !nzchar(Sys.which("awk")))
    skip("needs bash, gzip and awk")
  d <- tempfile("prefilter"); dir.create(d)
  writeLines(c("10\t12", "20\t25"), file.path(d, "regions.txt"))
  ## a few reads inside the regions, then many past the last one: an awk that stops
  ## early leaves gzip writing into a closed pipe
  pat <- c("chr1\t9\tCCC\t1", "chr1\t21\tTT\t2", "chr1\t30\tCC\t1",
           rep("chr2\t500\tCCCC\t1", 200000))
  con <- gzfile(file.path(d, "x.pat.gz"), "w"); writeLines(pat, con); close(con)
  awk <- system.file("scripts", "prefilter_pat.awk", package = "fragcorr")
  cmd <- sprintf("set -o pipefail; gzip -dc %s | awk -v RF=%s -f %s > %s 2>/dev/null",
                 shQuote(file.path(d, "x.pat.gz")), shQuote(file.path(d, "regions.txt")),
                 shQuote(awk), shQuote(file.path(d, "out.pat")))
  expect_equal(system2("bash", c("-c", shQuote(cmd))), 0L)
  expect_equal(readLines(file.path(d, "out.pat")), c("chr1\t9\tCCC\t1", "chr1\t21\tTT\t2"))
})

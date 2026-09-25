## The correlation patterns.
##
## rho(d) = phi * exp(-d / s): phi is how correlated two CpGs are at zero
## distance, s how fast that decays in base pairs.  The grid is the one the
## original used, and the column names it produced ("1_30", "0.6_1200", ...) are
## kept so results stay comparable.
##
## Note it is not a full crossing: s = 1200 carries only the three strongest phi.
## A long-range pattern that is also weak is indistinguishable from independence
## over a region of a few hundred base pairs, so those cells were left out.

fc_patterns <- function() {
  phi <- c(1, 0.8, 0.6, 0.4, 0.2, 0.1,
           1, 0.8, 0.6, 0.4, 0.2, 0.1,
           1, 0.8, 0.6, 0.4, 0.2, 0.1,
           1, 0.8, 0.6)
  s   <- c(rep(30, 6), rep(100, 6), rep(300, 6), rep(1200, 3))
  stopifnot(length(phi) == length(s))
  data.frame(phi = phi, s = s, name = paste0(phi, "_", s), stringsAsFactors = FALSE)
}

## Methylation bins.  The prior over patterns is estimated separately in each,
## because how methylated a region is changes which correlation structures are
## even possible: a region at 0.99 has almost no variation to correlate.
fc_bins <- function() {
  data.frame(
    label  = c("<0.1", "0.1-0.4", "0.4-0.6", "0.6-0.9", "0.9-1"),
    suffix = c("lt0.1", "0.1_0.4", "0.4_0.6", "0.6_0.9", "0.9_1"),
    stringsAsFactors = FALSE)
}

## Assign regions to bins by mean methylation, using the original's breaks.
assign_bin <- function(meth) {
  b <- fc_bins()
  cut(meth, breaks = c(-Inf, 0.1, 0.4, 0.6, 0.9, 1),
      labels = b$label, include.lowest = TRUE, right = TRUE)
}

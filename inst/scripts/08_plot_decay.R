#!/usr/bin/env Rscript
## Plots for 07_decay_curves.R: correlation against distance, one curve per
## methylation bin with its band (--out), and phi and s by bin (--out2).
## Plain base graphics.
##
##   Rscript 08_plot_decay.R --inp curves.rds --out curves.png [--out2 phi_s.png] [--dmax 400]
##
## Height at d = 0 is phi, the strength.  How fast the line falls is s, the
## range.  Plotting the curve rather than either number is the point: the bins
## do not rank the same way at every distance.

suppressMessages(library(fragcorr))
a <- fc_args(list(inp = "", out = "", out2 = "", dmax = "400"))
fc_require(a, "inp", "08_plot_decay.R --inp <curves.rds> [--out <png>] [--out2 <png>]")
x <- readRDS(a$inp)
band <- as.data.frame(x$band)

bins <- c("<0.1", "0.1-0.4", "0.4-0.6", "0.6-0.9", "0.9-1")
labs <- c("<0.1  hypo", "0.1-0.4", "0.4-0.6  hemi", "0.6-0.9", "0.9-1  hyper")
cols <- colorRampPalette(c("blue", "red"))(5)
ltys <- 1:5

dmax <- as.numeric(a$dmax)
sub  <- band[band$bin %in% bins & band$d <= dmax, ]
get  <- function(b) { z <- sub[sub$bin == b, ]; z[order(z$d), ] }

## ------------------------------------------------------------------ the curves
if (nzchar(a$out)) png(a$out, width = 1400, height = 1000, res = 150)
plot(NA, xlim = c(0, dmax), ylim = c(0, max(sub$hi)),
     xlab = "distance between two CpGs (bp)", ylab = "correlation rho(d)",
     main = "spatial correlation by methylation bin")
for (j in seq_along(bins)) {
  b <- get(bins[j])
  polygon(c(b$d, rev(b$d)), c(b$lo, rev(b$hi)),
          col = adjustcolor(cols[j], alpha.f = 0.2), border = NA)
}
for (j in seq_along(bins)) {
  b <- get(bins[j])
  lines(b$d, b$mean, col = cols[j], lty = ltys[j], lwd = 2)
}
legend("topright", legend = labs, col = cols, lty = ltys, lwd = 2, title = "methylation")
if (nzchar(a$out)) { dev.off(); cat("wrote ", a$out, "\n", sep = "") }

## ------------------------------------------------------- phi and s, bin by bin
phi <- colMeans(x$phi, na.rm = TRUE)[bins]
ss  <- colMeans(x$s,   na.rm = TRUE)[bins]
phi_ci <- apply(x$phi[, bins, drop = FALSE], 2, quantile, c(.025, .975), na.rm = TRUE)
s_ci   <- apply(x$s[,   bins, drop = FALSE], 2, quantile, c(.025, .975), na.rm = TRUE)

if (nzchar(a$out2)) png(a$out2, width = 1400, height = 700, res = 150)
par(mfrow = c(1, 2))
plot(seq_along(bins), phi, type = "b", pch = 19, ylim = c(0, 1), xaxt = "n",
     xlab = "methylation bin", ylab = "phi  (strength at d = 0)", main = "strength")
axis(1, at = seq_along(bins), labels = bins, cex.axis = 0.8)
arrows(seq_along(bins), phi_ci[1, ], seq_along(bins), phi_ci[2, ],
       angle = 90, code = 3, length = 0.05)
plot(seq_along(bins), ss, type = "b", pch = 19, ylim = c(0, max(s_ci)), xaxt = "n",
     xlab = "methylation bin", ylab = "s  (decay length, bp)", main = "range")
axis(1, at = seq_along(bins), labels = bins, cex.axis = 0.8)
arrows(seq_along(bins), s_ci[1, ], seq_along(bins), s_ci[2, ],
       angle = 90, code = 3, length = 0.05)
par(mfrow = c(1, 1))
if (nzchar(a$out2)) { dev.off(); cat("wrote ", a$out2, "\n", sep = "") }

## ----------------------------------------------------------------- the numbers
cat("\n")
for (j in seq_along(bins)) {
  b <- get(bins[j])
  cat(sprintf("%-8s rho(0)=%.3f [%.3f,%.3f]  rho(100)=%.3f [%.3f,%.3f]  rho(400)=%.3f  phi=%.3f  s=%5.1f\n",
              bins[j], b$mean[1], b$lo[1], b$hi[1],
              b$mean[b$d == 100], b$lo[b$d == 100], b$hi[b$d == 100],
              b$mean[b$d == dmax], phi[j], ss[j]))
}
dd <- sort(unique(sub$d))
cross <- function(u, v) { s <- sign(u - v); k <- which(diff(s) != 0)[1]
  if (is.na(k)) NA else round(dd[k] + (dd[k+1]-dd[k]) * abs(u[k]-v[k]) /
                              (abs(u[k]-v[k]) + abs(u[k+1]-v[k+1]))) }
cat("\ncrossings against <0.1 (bp):\n")
u <- get("<0.1")$mean
for (j in 2:5) if (nrow(get(bins[j])) && length(u))
  cat(sprintf("  vs %-8s %s\n", bins[j], cross(u, get(bins[j])$mean)))

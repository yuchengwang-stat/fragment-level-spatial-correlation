#!/usr/bin/env Rscript
## Step 2 -- region log-likelihoods under every correlation pattern, one sample.
##
## Reads the sample's .pat (ideally prefiltered to the regions; see
## loglik_one_sample.sh) and .beta, scores every fragment under all patterns in
## one pass, and writes one row per region.
##
##   Rscript 02_loglik.R --regions regions.rds --pat sample.pat[.gz] --beta sample.beta
##     --out sample.rds [--threads 1] [--min-n 5] [--use-bvn] [--chrom chr1,chr2]
suppressMessages(library(fragcorr))
a <- fc_args(list(regions = "", pat = "", beta = "", out = "", chrom = "",
                  min_n = "5", use_bvn = "false", threads = "1"))
fc_require(a, c("regions", "pat", "beta", "out"),
           "02_loglik.R --regions <rds> --pat <pat> --beta <beta> --out <rds>")
data.table::setDTthreads(as.integer(a$threads))

reg <- readRDS(a$regions)
chrom <- if (nzchar(a$chrom)) strsplit(a$chrom, ",")[[1]] else NULL
if (!is.null(chrom)) reg <- reg[reg$chr %in% chrom, , drop = FALSE]
if (!nrow(reg)) stop("no regions left after the chromosome filter")
cat(nrow(reg), " regions\n", sep = "")

beta <- read_beta_counts(a$beta)
pat  <- read_pat(a$pat, chrom = chrom)
cat("  ", nrow(pat), " read patterns\n", sep = "")

t0 <- Sys.time()
lik <- sample_loglik(reg, pat, beta, min_n = as.integer(a$min_n),
                     use_bvn = identical(tolower(a$use_bvn), "true"),
                     threads = as.integer(a$threads))
el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("  %d distinct fragments (%d reads) scored in %.2f min\n",
            attr(lik, "n_fragments"), attr(lik, "n_reads"), el))

saveRDS(list(loglik = lik, patterns = fc_patterns(), n_fragments = attr(lik, "n_fragments"),
             n_reads = attr(lik, "n_reads"), minutes = el, sample = basename(a$pat)), a$out)
cat("wrote ", a$out, "  (", nrow(lik), " regions)\n", sep = "")
print(table(lik$bin, useNA = "ifany"))

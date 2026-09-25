\name{fragcorr-package}
\alias{fragcorr-package}
\alias{fragcorr}
\alias{fc_patterns}
\alias{fc_bins}
\alias{assign_bin}
\alias{read_cpg_list}
\alias{merge_positions}
\alias{build_regions}
\alias{validate_regions}
\alias{load_cpg_index}
\alias{index_regions}
\alias{validate_indexed}
\alias{read_beta_counts}
\alias{read_pat}
\alias{fragments_in_regions}
\alias{loglik_all_patterns}
\alias{sample_loglik}
\alias{combine_loglik}
\alias{fit_prior_bin}
\alias{fit_priors}
\alias{posterior_from_prior}
\alias{summarise_posterior}
\alias{region_tables}
\alias{decay_curves}
\alias{read_posteriors}
\alias{simulate_dataset}
\alias{fc_args}
\alias{fc_require}
\alias{fc_script}
\docType{package}
\title{Fragment-Level Spatial Correlation of CpG Methylation}
\description{
Estimates, per region and cell type, a posterior over the exponential kernel
\eqn{\rho(d) = \phi \exp(-d/s)} that describes how the methylation states of
CpGs on the same sequenced fragment co-vary with distance, on the latent probit
scale and conditional on each CpG's own methylation rate. The model and the
pipeline are described in the package README and in \file{docs/METHOD.md} of
the source repository; each function's arguments are documented in comments at
its definition.
}
\section{Functions}{
\describe{
  \item{Kernel grid and bins}{\code{fc_patterns()}: the 21 \eqn{(\phi, s)} kernels.
    \code{fc_bins()}, \code{assign_bin(meth)}: the five methylation bins.}
  \item{Regions}{\code{read_cpg_list()}, \code{merge_positions()}, \code{build_regions()},
    \code{validate_regions()}: regions from a CpG list.
    \code{load_cpg_index()}, \code{index_regions()}, \code{validate_indexed()}:
    map regions onto the CpG index (row k = \code{.beta} row k).}
  \item{Reads}{\code{read_beta_counts()}, \code{read_pat()}: wgbstools inputs.
    \code{fragments_in_regions()}: reads trimmed to regions with their counts.}
  \item{Likelihood}{\code{loglik_all_patterns()}: per-fragment log-likelihoods under
    every kernel (C++, EP). \code{sample_loglik()}: region log-likelihoods for one
    sample.}
  \item{Prior and posterior}{\code{combine_loglik()}, \code{fit_prior_bin()},
    \code{fit_priors()}, \code{posterior_from_prior()}, \code{summarise_posterior()}.}
  \item{Results}{\code{region_tables()}: long and wide result tables.
    \code{decay_curves()}: posterior-averaged \eqn{\rho(d)} per methylation bin.
    \code{read_posteriors()}: step-3 files as a named list.}
  \item{Other}{\code{simulate_dataset()}: a synthetic dataset in wgbstools format with
    known kernels. \code{fc_args()}, \code{fc_require()}, \code{fc_script()}:
    helpers for the command-line scripts in \code{system.file("scripts", package = "fragcorr")}.}
}
}
\examples{
fc_patterns()
d <- simulate_dataset(file.path(tempdir(), "fc_demo"), n_loci = 20L, samples_per_type = 1L)
reg <- index_regions(build_regions(read_cpg_list(d$cpg_list)),
                     load_cpg_index(d$cpg_index), verbose = FALSE)
head(reg)
}
\keyword{package}

## fragcorr: fragment-level spatial correlation of CpG methylation.
##
## See README.md for what the method estimates and how to run it, and
## docs/METHOD.md for the model.

## Columns referred to by data.table's non-standard evaluation.  Declared so that
## R CMD check does not report them as undefined globals.
utils::globalVariables(c(
  ".", ".N", ".SD", ":=", "chr", "pos", "idx", "interval_id", "regStart", "regEnd",
  "startCpG", "endCpG", "ov_start", "ov_end", "fragment", "trimmed", "n_obs",
  "count", "cpg_gaps", "dist_sub", "M", "N", "M_sum", "N_sum", "meth", "bin",
  "start", "end", "width", "n_cpg", "n_cpg_in_region", "n_list_expected",
  "n_list_cpg", "list_cpgs", "list_category", "category", "region", "celltype",
  "best_phi", "best_s", "best_post", "post_winner", "phi_mean", "s_mean",
  "pattern", "..pattern_names", "layout", "p", "tgt", "r2_split", "n_low", "n_high",
  "meth_sd", "meth_jump", "meth_range", "n_cov", "target_lowhigh", "region_id", "meth_all", "meth_all_sd", "meth_all_iqr",
  "meth_target", "meth_target_sd", "meth_target_iqr", "enough", "parent_region", "single_target", "row"))

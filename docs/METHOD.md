# Method

This document states the model and each computational step in full. Code
references are to the `fragcorr` package.

## 1. Data

- **Reads.** For each sample, `.pat.gz` holds its bisulfite reads in the
  wgbstools format. A read is a start CpG index `i0`, a string over consecutive
  CpGs `i0, i0+1, ...` with one character per CpG, and a count of identical
  reads. The characters are `C` (methylated), `T` (unmethylated) and `.` (not
  called).
- **Methylation counts.** For each sample, `.beta` holds `(M_j, N_j)` for every
  CpG `j` of the genome-wide index: methylated and total calls at that CpG,
  tallied from the same reads. It is stored as uint8, so depth is capped at
  255.
- **CpG index.** The CpG index table gives the coordinates of every CpG, in
  `.beta` row order.

A fragment here is one read, trimmed to the region it overlaps.

## 2. Regions

A list of CpGs defines where to look (`build_regions()`).

1. Each listed CpG at position `x` becomes the interval `[x - pad, x + pad]`.
2. Overlapping intervals are merged into clusters.
3. A cluster wider than `max_width` is cut into the fewest contiguous runs of
   listed CpGs of near-equal size that each fit. The boundary between two
   pieces is the midpoint of the gap that separates them, so pieces abut and do
   not overlap.
4. Regions with fewer than `min_cpgs` listed CpGs are dropped.

Every region is then mapped onto the CpG index (`index_regions()`). It gets the
index range `[startCpG, endCpG]` of all CpGs inside it, listed or not, and the bp
gaps between consecutive CpGs. The likelihood uses **every CpG in the region**,
not only the listed ones.

Defaults: `pad = 250`, `max_width = 1000`, `min_cpgs = 2`.

**Refinement (`refine_regions()`, step 1b).** Methylation is pooled per cell type
over the `.beta` files. A run of consecutive listed CpGs conflicts in a cell type
when, among its listed CpGs with pooled depth >= 20, one is at <= 0.2 and another
at >= 0.8. A run that conflicts in at least `min_celltypes` cell types (default 5)
is cut between two consecutive listed CpGs. The cut is the one leaving the fewest
conflicts on its two sides, with ties broken by the widest bp gap, and the two
sides are checked recursively. Each piece spans its listed CpGs padded by `pad`,
within the original region, and neighbouring pieces meet at the midpoint between
them. The split does not depend on the cell type being fitted, so all cell types
share one set of regions. Requiring several cell types keeps a single noisy
per-CpG estimate from causing a cut.

## 3. The model for one fragment

Take a fragment in sample `r`. Let `j = 1..K` be its usable CpGs: called `C` or
`T`, with `.beta` depth `N_j >= min_n` (default 5) in that sample. Let `x_j` be
their positions.

**Marginals.**

```
p_j = M_j / N_j,   clamped to [1e-6, 1 - 1e-6]
```

`p_j` is that CpG's methylation rate in sample `r`.

**Latent Gaussian.**

```
Z = (Z_1..Z_K) ~ N(0, R),   R_jj = 1,   R_jk = rho(|x_j - x_k|)
rho(d) = phi * exp(-d / s)          (d capped at 7000 bp)
CpG j reads C  <=>  Z_j <= qnorm(p_j)
```

Each coordinate is standard normal, so `P(CpG j reads C) = p_j` exactly,
whatever `phi` and `s` are. The kernel only governs how calls co-occur. For
distinct positions the exponential kernel is positive definite, so `R` is a
valid correlation matrix for every `phi` in `(0, 1]`.

**Fragment likelihood.** The probability of the observed string is the Gaussian
mass of a box:

```
L(phi, s) = P( Z_j <= u_j for C sites,  Z_j > u_j for T sites ),   u_j = qnorm(p_j)
```

Flipping the sign of each `T` coordinate turns this into one orthant
probability, `P(W <= v)`, where `W ~ N(0, D R D)` and `D` is the diagonal matrix
of signs.

- `K = 1`: `L` is `p_j` or `1 - p_j`, the same for every kernel.
- `K >= 2`: `L` is computed by expectation propagation (EP) for the multivariate
  normal CDF. `--use-bvn` replaces EP by the exact bivariate normal CDF at
  `K = 2`.

Log-likelihoods are rounded to three decimals, as in the original
implementation.

## 4. Region and cell-type likelihood

- **Within a region.** Fragments are treated as independent. The region's
  log-likelihood under kernel `k`, in sample `r`, is

  ```
  l_irk = sum over distinct trimmed fragments f in region i of  count_f * log L_f(k)
  ```

  Reads with identical trimmed strings over the same CpG span share one
  likelihood and are collapsed.
- **Within a cell type.** Samples of the same cell type are treated as
  independent draws of the same kernel, so their log-likelihoods add:

  ```
  l_ik = sum over samples r of l_irk
  ```

- **Region methylation.** The region's methylation in a cell type is

  ```
  m_i = sum M / sum N
  ```

  over the CpG spans of its fragments, pooled over the cell type's samples. It
  is a fragment-weighted average, and is used only to choose the prior bin.

## 5. Kernel grid

The 21 candidate kernels `(phi, s)` are:

- `phi` in {1, 0.8, 0.6, 0.4, 0.2, 0.1}, combined with `s` in {30, 100, 300} bp;
- and `s` = 1200 bp with `phi` in {1, 0.8, 0.6}.

A weak kernel with a very long range is nearly indistinguishable from
independence over a region of a few hundred bp. Those combinations were
therefore left out of the grid.

## 6. Prior: empirical Bayes by methylation bin

**Why the prior depends on methylation.** A region at 0.99 methylation has
almost no variation to correlate, so which kernels are plausible depends on the
level. Regions are therefore binned by `m_i`, with breaks at 0.1, 0.4, 0.6 and
0.9.

**The fit.** Within each bin `b` of each cell type, the prior weights
`pi_b = (pi_b1..pi_b21)` maximise the marginal likelihood of the bin's regions:

```
pi_b = argmax over the simplex of  sum over regions i in b of  log( sum_k pi_k exp(l_ik) )
```

This is the nonparametric maximum likelihood estimate (NPMLE) of the
distribution of kernels in the bin, restricted to the grid. It is solved by
`mixsqp`, which works on the log-likelihood matrix directly. At convergence, the
average of the regions' posteriors equals the fitted prior. The bin-level mean
of `phi` is therefore the NPMLE's mean, and is not biased toward any fixed
prior.

**Options.**

- `--prior`: fix the prior from a file, for example one shared by all cell
  types.
- `--max-fit`: fit on a random subset of regions per bin.

## 7. Posterior and summaries

```
post_ik = pi_bk exp(l_ik) / sum_k' pi_bk' exp(l_ik')        (computed in log space)
```

The reported summaries are:

- the argmax kernel, giving `phi` and `s`;
- its probability, `post_winner`;
- the posterior means, `phi_mean = sum_k post_ik phi_k` and
  `s_mean = sum_k post_ik s_k`;
- the effective number of kernels, `exp(entropy of post_i)`. It is 1 when the
  data pick one kernel and 21 when they carry no information.

## 8. Decay curves

A region's posterior implies a correlation curve:

```
rho_i(d) = sum_k post_ik * phi_k * exp(-d / s_k)
```

This is a mixture of exponentials. It is not the same as plugging `phi_mean` and
`s_mean` into one exponential: the tails differ.

`decay_curves()` averages `rho_i(d)` over the regions of each methylation bin
within each cell type. It then takes a percentile bootstrap across cell types.
The band answers whether a pattern holds across cell types; regions within one
cell type are not independent replicates.

## 9. Uniformity of a region's methylation

`heterogeneity_stats()` and `09_heterogeneity.R` measure how uniform a region's
methylation is, per cell type. For each CpG in the region, counts are pooled
over the cell type's samples, and CpGs with pooled depth below `min_n` (default
20) are dropped. With `p_1..p_n` the rates of the remaining CpGs in genome
order:

- `meth_all`, `meth_all_sd` and `meth_all_iqr` are the mean, SD and
  interquartile range of `p`. `meth_target`, `meth_target_sd` and
  `meth_target_iqr` are the same over the listed CpGs only.
  `meth_range` is `max(p) - min(p)`.
- `meth_jump = mean |p_(j+1) - p_j|`.
- `r2_split = 1 - min_c SSE_c / SST`. Here `SSE_c` is the within-piece sum of
  squares after cutting between CpG `c` and `c+1`, and `SST` is the total sum of
  squares. It is 1 for two flat pieces, and low when high and low CpGs are
  interleaved.
- `layout` is:
  - `mixed_split` or `mixed_interleaved`, if at least two CpGs are at
    `<= 0.2` and two at `>= 0.8` (with `r2_split >= 0.75` for split);
  - otherwise `uniform`, if the SD of `p` is `<= 0.10`;
  - otherwise `graded`.

All cut-offs are arguments. Binomial noise contributes to the SD: at
`p = 0.5` and depth 20, its SD alone is about 0.11.

## 10. Assumptions, and what they imply

- **Marginals are plugged in.** The fit is conditional on each CpG's per-sample
  rate from the `.beta`, which comes from the same reads. This is a two-step
  estimator, like fitting the margins first in a copula model. Simulation shows
  it recovers the kernel across methylation levels from 0.05 to 0.95. Marginals
  estimated independently of the reads would bias `phi` down at extreme levels;
  see VALIDATION.md.
- **One kernel per region.** The kernel is assumed stationary within a region.
  A region spanning two domains gets a compromise.
- **Interleaved extremes carry little information.** A region whose adjacent
  CpGs alternate between very low and very high rates says little about the
  kernel: such pairs are nearly uninformative. In that case, EP's small
  per-fragment error can dominate the likelihood and push `phi` up. A region
  split into a low half and a high half is fine. See VALIDATION.md, section 5.
- **Fragments are independent**, and so are samples within a cell type.
- **Latent scale.** `phi` is a correlation of latent Gaussians, a tetrachoric
  correlation. The implied correlation of the 0/1 calls is smaller, and depends
  on `p_j` and `p_k`.

# Validation

Every script named here is in `validation/` or `tests/`. Each runs in minutes on
a laptop, with `fragcorr` installed. The numbers below come from those runs.

## 1. The likelihood engine reproduces the original implementation

`validation/test_equivalence.R` compiles the original engine
(`validation/reference/fragment_lik.cpp`) next to the package's engine. It then
scores 240 synthetic fragments under all 21 kernels with both.

- The fragment set covers every code path: fragments with K = 0, 1, 2 and >= 3
  usable sites, depths from 1 to 60, and uncalled sites.
- The original computes one kernel per pass. The package computes all kernels in
  one pass, from a single parse of each fragment.

```
rebuild, defaults, vs original                 IDENTICAL  (max |diff| 0, NA pattern same)
rebuild with use_bvn=TRUE vs original          DIFFERS    (max |diff| 0.013, NA pattern same)
```

The second line is expected. `use_bvn` replaces EP by the exact bivariate
normal CDF for two-site fragments, and so changes the numbers. It is off by
default.

## 2. Unit tests

`tests/testthat` checks:

- region merging, splitting and index mapping;
- that a one-site fragment scores the same under every kernel;
- that sites below the depth threshold are dropped;
- that concordant fragments favour strong kernels and discordant ones weak
  kernels;
- that the two-site closed form stays within 0.05 of EP;
- an end-to-end run on simulated data that must separate a strongly from a
  weakly correlated cell type.

## 3. Recovery of known kernels, end to end

`simulate_dataset()` writes a genome of CpGs, a CpG list, and per-sample
`.pat.gz` and `.beta` files. Reads are drawn from the model, and the `.beta`
counts are tallied from the same reads. The data have:

- two cell types with known kernels;
- three samples each;
- 100 loci, at methylation levels 0.03 to 0.97.

The full pipeline recovers:

| reads per locus per sample | A (truth phi 0.8, s 100) | B (truth phi 0.3, s 300) |
|---|---|---|
| 80   | phi 0.76, s 110 | phi 0.31, s 291 |
| 300  | phi 0.75, s 109 | phi 0.29, s 262 |
| 1000 | phi 0.78, s 100 | phi 0.28, s 280 |

The small shortfall for A comes from averaging posteriors over a grid whose
neighbours of 0.8 are 0.6 and 1.

## 4. The estimate does not depend on the methylation level

`validation/test_bias_sources.R` fixes the truth at `phi = 0.6` and `s = 100`.
It simulates 300 regions (8 CpGs, 60 fragments each) at each of seven
methylation levels, and runs them through the engine, the `mixsqp` prior and the
posterior.

The same fragments are scored three times. Only the marginal rates handed to the
likelihood change between runs:

- **indep**: `M ~ Binomial(30, p)`, drawn independently of the fragments;
- **reads**: `M` and `N` counted from the simulated fragments themselves, as a
  real `.beta` is;
- **oracle**: the true `p`.

Average posterior mean of `phi`:

| methylation | indep | reads | oracle |
|---|---|---|---|
| 0.05 | 0.373 | **0.572** | 0.547 |
| 0.15 | 0.529 | **0.591** | 0.591 |
| 0.30 | 0.562 | **0.606** | 0.605 |
| 0.50 | 0.576 | **0.598** | 0.599 |
| 0.70 | 0.574 | **0.607** | 0.604 |
| 0.85 | 0.525 | **0.594** | 0.587 |
| 0.95 | 0.412 | **0.584** | 0.568 |

Average posterior mean of `s` in the reads column: 99 to 116 bp.

**How to read this.**

- **reads matches real data**, and recovers the truth at every level, to within
  0.03.
- **indep shows a spurious drop of about 0.2 at the extremes.** When the plugged-in
  rate disagrees with the reads (at `p = 0.05` and depth 30, `M = 0` one time
  in five), the model explains the mismatch as weaker correlation. An earlier
  version of this check, `validation/test_phi_vs_marginal.R`, used the indep
  scheme and overstated the bias. It is kept, with a note, for the record.
- **oracle is within 0.05 of the truth.** What remains comes from EP and the
  limited information at extreme levels.

`validation/baseR_reference.R` is an independent base-R implementation. It
computes orthant probabilities by GHK simulation instead of EP, and fits the
prior by EM instead of `mixsqp`, using true marginals. It gives `phi` of 0.53,
0.59 and 0.52, and `s` of 153, 100 and 182, at methylation 0.05, 0.5 and 0.95.
That is the same shape as the oracle column, on a smaller sample (30 regions of
40 fragments).

## 5. Regions whose CpGs do not share one methylation rate

`validation/test_mixed_layouts.R` uses the same set-up (truth `phi = 0.6`,
`s = 100`, marginals from the reads). The per-CpG rate follows a layout within
each region:

| layout | mean rate | phi | s |
|---|---|---|---|
| uniform 0.475 | 0.475 | 0.603 | 100.5 |
| uniform 0.05 | 0.05 | 0.568 | 103.6 |
| uniform 0.90 | 0.90 | 0.555 | 100.3 |
| split: 4 CpGs at 0.05, then 4 at 0.90 | 0.475 | 0.605 | 95.4 |
| alternating 0.05, 0.90, 0.05, ... | 0.475 | **0.895** | 86.6 |

A region split into a low half and a high half is estimated correctly, because
every CpG has its own threshold. A region whose **adjacent** CpGs alternate
between extremes is not.

`validation/test_ep_vs_exact.R` shows why. It scores the same 500 fragments with
EP and with a precise orthant probability (mvtnorm, Genz-Bretz) and reports each
kernel's log-likelihood minus that of the truth:

| layout | method | phi 0.4 | phi 0.6 (truth) | phi 0.8 | phi 1 |
|---|---|---|---|---|---|
| uniform | EP | -4.2 | 0 | -15.1 | -58.5 |
| uniform | exact | -4.5 | 0 | -14.1 | -55.2 |
| split | EP | -0.8 | 0 | -4.9 | -18.3 |
| split | exact | -1.4 | 0 | -3.8 | -15.1 |
| alternating | EP | -3.9 | 0 | **+2.0** | **+1.9** |
| alternating | exact | -3.2 | 0 | +0.5 | -2.2 |

Two things combine in the alternating layout:

- **The likelihood is nearly flat.** A pair of CpGs at 0.05 and 0.90 says little
  about their correlation: across `phi` from 0.4 to 1 the exact log-likelihood
  moves by about 5 units, against 55 for a uniform region.
- **EP tilts that flat likelihood toward high `phi`.** Its per-fragment error is
  small, but on a flat likelihood it decides the answer.

The region's mean methylation also places it in the 0.4 to 0.6 bin, although
none of its CpGs is intermediate. Such regions should be flagged rather than
interpreted.

## Rerunning

```bash
Rscript validation/test_equivalence.R
Rscript validation/test_bias_sources.R        # a few minutes
Rscript validation/test_mixed_layouts.R
Rscript validation/test_ep_vs_exact.R         # needs mvtnorm; a few minutes
Rscript validation/baseR_reference.R
Rscript -e 'testthat::test_dir("tests/testthat")'
```

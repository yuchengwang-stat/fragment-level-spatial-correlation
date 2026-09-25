# Fragment-level spatial correlation of CpG methylation

`fragcorr` is an R package. It estimates how strongly the methylation states of
CpGs **on the same sequenced DNA fragment** move together, and how that
dependence decays with the distance between them. It works region by region and
cell type by cell type, from bisulfite sequencing reads in the
[wgbstools](https://github.com/nloyfer/wgbs_tools) `.pat` / `.beta` formats.

For every region and cell type, the result is a posterior distribution over an
exponential correlation kernel

```
rho(d) = phi * exp(-d / s)
```

- `d` is the distance between two CpGs in base pairs.
- `phi` is the correlation at distance zero: how strongly neighbouring CpGs are
  co-methylated.
- `s` is the decay length in bp: how far the dependence reaches. The correlation
  falls by a factor of `e` every `s` bp.

## What kind of correlation this is

- **Within a molecule, not across samples.** Each read comes from one DNA
  molecule. The question is whether, on that molecule, a methylated CpG makes its
  neighbours more likely to be methylated. Nothing here compares methylation
  levels between samples or individuals, and different fragments are treated as
  independent.
- **On the latent (probit) scale.** Behind every CpG on a fragment sits a
  standard normal variable, and the CpG reads as methylated when that variable
  falls below a threshold. `rho` is the correlation of these latent variables,
  a tetrachoric correlation. The Pearson correlation of the 0/1 calls is smaller,
  and depends on the methylation level as well.
- **Net of each CpG's own methylation rate.** Each CpG's threshold is fixed by
  its methylation rate in that sample, read from the `.beta` file. `phi` and `s`
  describe only how the calls co-occur beyond what those rates imply. A high
  `phi` means reads tend to be methylated or unmethylated as blocks: the
  signature of heterogeneity between molecules, cells or alleles.

## How it works, in brief

1. **Regions.** A list of CpGs of interest defines the regions. Each CpG is
   padded (default +/-250 bp), overlapping pads are merged, clusters wider than
   a cap (default 1000 bp) are cut between CpGs, and regions holding fewer than
   two listed CpGs are dropped.
2. **Fragment likelihood.** Every read overlapping a region is trimmed to the
   region. Its C/T string is then scored under each of 21 candidate
   `(phi, s)` kernels. The score is the probability of that exact string under
   the latent Gaussian model, a multivariate normal orthant probability, computed
   by expectation propagation (EP).
3. **Aggregation.** Log-likelihoods are summed over the fragments of a region
   (weighted by read counts) and then over the samples of a cell type.
4. **Empirical Bayes.** Regions are grouped into five bins by methylation level.
   Within each bin, a prior over the 21 kernels is fitted by maximising the
   marginal likelihood of that bin's regions, using `mixsqp`.
5. **Posterior.** Each region and cell type gets a posterior over the 21 kernels,
   reported as the winning `(phi, s)`, the winner's posterior probability, and
   posterior means of `phi` and `s`.

The full model is in [docs/METHOD.md](docs/METHOD.md). The checks behind it are in
[docs/VALIDATION.md](docs/VALIDATION.md).

## Installation

Requirements:

- R >= 4.1 and a C++ compiler;
- the R packages `data.table`, `mixsqp`, `Rcpp` and `RcppArmadillo`;
- `awk` and `gzip` for the read prefilter.

```bash
git clone https://github.com/yuchengwang-stat/fragment-level-spatial-correlation.git
R CMD INSTALL fragment-level-spatial-correlation
```

or, from R:

```r
remotes::install_github("yuchengwang-stat/fragment-level-spatial-correlation")
```

While the repository is private, `install_github` needs a GitHub token in
`GITHUB_PAT`.

On a cluster, compile inside a job rather than on a login node.

## Quick start: a simulated example

The package ships a small simulator and an end-to-end script. The simulated
data have two cell types with known kernels:

- A: `phi = 0.8`, `s = 100`;
- B: `phi = 0.3`, `s = 300`.

Each cell type has three samples, and the whole run takes about ten seconds.

```bash
FC=$(Rscript -e 'cat(system.file("examples", package = "fragcorr"))')
bash "$FC/run_example.sh" my_example
```

It ends by printing what the pipeline recovered:

```
   celltype regions phi_mean s_mean
1:        A     111    0.733  114.5
2:        B     111    0.296  256.0
truth: A phi 0.8 s 100 | B phi 0.3 s 300
```

`run_example.sh` is also the shortest complete template for a real run.

## Inputs

| Input | What it is |
|---|---|
| CpG list | The CpGs that define the regions, as a CSV or TSV. It can have `chr` and `pos` columns, or one coordinate column such as `hg38` with values `chr1:99851871` or `chr1:99851871-99851871` (pick it with `--col`). An optional `Category` column can be filtered with `--category`. |
| CpG index | The genome-wide CpG table in `.beta` row order, where **row k is row k of every `.beta` file**. It can be an `.rds` data frame with `chr` and `pos`, or a text or BED file such as the wgbstools `CpG.bed.gz`. If it carries its own index column, that column must equal the row number. |
| `.pat.gz` | Per sample, wgbstools read patterns: chromosome, start CpG index, a C/T/. string, and read count. The file is sorted by start. |
| `.beta` | Per sample, wgbstools counts: a uint8 pair (methylated, total) for every CpG of the index. It must come from the same reads as the `.pat`. |
| Sample sheet | A TSV with a header and the columns `sample`, `celltype`, `pat` and `beta`. It is used by the SLURM templates. |

Positions are 1-based coordinates of the C of each CpG, on the same genome build
as the CpG index.

## Running the pipeline

The scripts are installed with the package:

```bash
FC=$(Rscript -e 'cat(system.file("scripts", package = "fragcorr"))')
```

**Step 1: regions.** Run this once.

```bash
Rscript $FC/01_regions.R --cpg cpgs.csv --col hg38 --ref cpg_index.rds --out regions.rds \
  [--category HighValue,EPIC] [--pad 250] [--max-width 1000] [--min-cpgs 2]
```

**Step 2: log-likelihoods.** Run this once per sample. The wrapper first keeps
only the reads that can touch a region (an `awk` pass over the `.pat.gz`), then
scores them.

```bash
bash $FC/loglik_one_sample.sh regions.rds sample.pat.gz sample.beta lik/<celltype>/<sample>.rds 8
```

The last argument is the number of worker processes. Options of the scorer
itself (`02_loglik.R`) are:

- `--min-n` (default 5): minimum `.beta` depth for a CpG to be used;
- `--use-bvn`: exact bivariate normal for two-site fragments; see Limitations;
- `--chrom`: restrict to some chromosomes.

**Step 3: prior and posterior.** Run this once per cell type, or any group of
samples.

```bash
Rscript $FC/03_posterior.R --lik lik/<celltype> --out post/<celltype>.rds --label <celltype> \
  [--prior prior/<celltype>.rds]
```

If the `--prior` file exists, it is used instead of fitting. This lets several
groups share one prior. If it does not exist, the fitted prior is written there.

**Step 4: result tables.**

```bash
Rscript $FC/04_region_table.R --regions regions.rds --cpg cpgs.csv --col hg38 \
  --post post --out results/region_ct
```

`--cpg`, `--col` and `--category` must be the same as in step 1.

All scripts accept both `--key value` and `--key=value`. An unknown option is an
error, so typos fail before a long job starts.

### On a SLURM cluster

`inst/slurm/` has two array templates:

- `loglik_array.sbatch`: one task per sample, step 2;
- `posterior_array.sbatch`: one task per cell type, step 3.

Set the partition, account and paths at the top of each. Step 2 is the expensive
part: its cost scales with the number of reads in the regions, and it
parallelises over fragments within a task. Both templates export
`OMP_NUM_THREADS=1`. SLURM does not set it, and without it BLAS and OpenMP see
every core on the node.

`loglik_one_sample.sh` writes a prefiltered `.pat` next to the output. When
scoring succeeds, it removes that file and its region index file by name. Set
`KEEP_PAT=1` to keep them.

## Outputs

**`<prefix>_long.csv`** has one row per region and cell type:

| Column | Meaning |
|---|---|
| `region_id`, `region`, `chr` | Region number and range (`chr:start-end`) |
| `n_list_cpg`, `list_cpgs` | The listed CpGs inside the region, separated by `;` |
| `celltype` | Cell type (the step-3 file name) |
| `phi`, `s` | The kernel with the largest posterior probability |
| `post_winner` | That largest posterior probability |
| `phi_mean`, `s_mean` | Posterior expectations over all 21 kernels |
| `meth` | The region's methylation in that cell type: the quantity that chose its prior bin |

**`<prefix>_wide.csv`** has one row per region. It holds the coordinates, the
width, all CpGs in the region, the listed CpGs and their categories, and then
`<celltype>_phi`, `<celltype>_s` and `<celltype>_meth` for each cell type. A
cell type with no fragment in a region has blank cells there.

The step-3 `.rds` files hold the full posterior for each region (one `_post`
column per kernel), the fitted prior, and the sample list.

**Diagnostics:**

- `05_compare.R`: how informative the data are by methylation bin, and how
  often cell types agree on a region's winning kernel.
- `06_prior_check.R`: whether disagreement between cell types survives a single
  shared prior.
- `07_decay_curves.R` and `08_plot_decay.R`: the posterior-averaged curve
  `rho(d)` per methylation bin, with a bootstrap band across cell types.

## Validation

In short (details in [docs/VALIDATION.md](docs/VALIDATION.md)):

- **The likelihood engine reproduces the original implementation exactly.**
  On 240 synthetic fragments that exercise every code path, the maximum
  difference over the 21 kernels is 0.
- **The estimator recovers known kernels on simulated data.** The simulated data
  follow the model, with `.beta` counts tallied from the same reads.
- **The estimate does not depend on the methylation level.** With a fixed truth
  of `phi = 0.6` and `s = 100`, and the methylation rate set between 0.05 and
  0.95, the recovered `phi` is 0.57 to 0.61 and `s` is 99 to 116 bp.
- **Regions whose adjacent CpGs alternate between extremes are the exception**:
  there the likelihood is flat and EP biases `phi` upward (section 5 of VALIDATION.md).
- **Unit tests** cover region building, indexing and the likelihood's limiting
  cases (`tests/testthat`).

## Reading the results, and limitations

- **The kernel grid is coarse.** `phi` takes values in
  {1, 0.8, 0.6, 0.4, 0.2, 0.1} and `s` in {30, 100, 300, 1200} bp (1200 only with
  the three largest `phi`). `phi_mean` and `s_mean` are smoother summaries than
  the argmax. Where the winning posterior is low, prefer them.
- **`phi` and `s` trade off when data are thin.** A weak long-range kernel and a
  stronger short-range one predict similar reads at typical distances. Regions
  that are narrow, CpG-poor or shallowly covered pin down `s` poorly.
- **Two-site fragments use EP, like the original.** EP is approximate there: it
  is off by up to about 0.03 in probability. `--use-bvn` switches to the exact
  closed form, which is more accurate but changes the numbers relative to the
  original method.
- **Regions whose adjacent CpGs alternate between low and high methylation are
  unreliable.** Their reads carry little information about the kernel, EP then
  tilts `phi` upward, and their average methylation puts them in the middle bin.
  A region split into a low half and a high half is fine. See section 5 of
  [docs/VALIDATION.md](docs/VALIDATION.md).
- **Groups with few samples are noisier.** Each group fits its own prior, and a
  group with one or two samples fits it from little data. `06_prior_check.R`
  shows how much this matters for your data.
- **Methylation rates are plug-in estimates.** They come from each sample's
  `.beta`, capped at 255 reads per CpG by the format. This is accurate as long
  as the `.beta` and `.pat` come from the same reads, as wgbstools produces them.

## Repository layout

```
R/                  package code: regions, indexing, reads, likelihood, prior/posterior, tables, curves, simulator
src/                C++ likelihood engine (EP orthant probabilities)
inst/scripts/       command-line steps 01-08 and the per-sample wrapper
inst/slurm/         SLURM array templates
inst/examples/      run_example.sh: the full pipeline on simulated data
tests/testthat/     unit and end-to-end tests
validation/         equivalence with the original engine, bias simulations, a base-R reference
docs/               METHOD.md, VALIDATION.md
```

## Credits

The expectation-propagation routine for multivariate normal orthant
probabilities (`src/ep_core.h`, `src/A_aux*.h`) is derived from the EPmvnCDF
package by Francesco Denti and Augusto Fasano. That code is under the MIT
licence; see `inst/COPYRIGHTS`. The rest of the package is by Yucheng Wang; see
`LICENSE`.

# Demonstration and example usage of bmfaToolkits

### 0. Install / load

Many methods are implemented in **optional backends**
(GitHub/Bioconductor packages, or bundled scripts).

``` r
# install.packages("remotes")
# remotes::install_local("path/to/bmfaToolkits")  # or install from your repo
library(bmfaToolkits)
```

Check what you can run on your machine:

``` r
available_backends()
```

If you need an optional backend (and it is unavailable above), install
it via:

``` r
install_backend("bmsfa")        # installs MSFA from GitHub
install_backend("momss")        # installs BFR.BE + Bioc deps
install_backend("sufa")         # installs SUFA from GitHub (build_vignettes = FALSE)
install_backend("curated_ovarian")
```

and no need t `library` it after installing.

> **Note on MOM-SS/BFR.BE.** We only support MOM-SS in “BFR.BE” backend
> right now. If you want to run other method mentioned in their paper,
> eg. Laplace-SS, please refer to their original repository.

> **Note on SUFA.**  
> You need to install PROJ, sqlite3 and GDAL onto PATH. Several updates
> should be done, for example, terra. These might be difficult.

> **Note on BMSFA** install_backend(“bmsfa”) installs the MSFA package
> from Mavis’s GitHub, which allows user-input scaling and centering
> options. The DeVito MSFA package does not have these options and does
> internal standardizations. If your backend is DeVito’s MSFA, do not
> specify scaling/centering.

> **Note on script backends (PFA / Tetris).**  
> In this package, PFA and Tetris are bundled under `inst/extdata/pfa/`
> and `inst/extdata/tetris/`.  
> They are loaded on demand via `load_backend("pfa")` or
> `load_backend("tetris")`. PFA may compile `PFA.cpp`, which requires a
> working C++ toolchain.

**If install_backend function does not work, you can always install them
yourself then run the bmfaToolkits package.**

### 1. Load toy data (shipped in the package)

We demonstrate the workflows using a toy simulated data with `S = 3`
studies and `P = 20` variables.

The toy data live in `inst/extdata/toy.rds`. We can load it using
[`system.file()`](https://rdrr.io/r/base/system.file.html).

``` r
toy <- readRDS(system.file("extdata", "toy.rds", package = "bmfaToolkits", mustWork = TRUE))
str(toy, max.level = 2)
```

The package wrappers expect `Y_list`: a list of study matrices with the
same number of columns (variables).

``` r
Y_list <- toy$Y_list
S <- length(Y_list)
P <- ncol(Y_list[[1]])
c(S = S, P = P)
```

Throughout this vignette: - `S` = number of studies - `P` = number of
variables

### 2. Common post-processing outputs

Each `postprocess_*()` function returns a list of covariance components
used in the tutorial:

- `Phi` / `SigmaPhi`: shared loadings / shared covariance
- `LambdaList` / `SigmaLambdaList`: study-specific loadings /
  covariances (if the method has them)
- `Psi` or `PsiList`: diagonal residual variances
- `SigmaMarginal`: per-study marginal covariance (shared + specific +
  residual)

The selection helpers choose dimensions using an **eigenvalue proportion
rule** (see
[`select_k_from_sigma()`](https://mavis-liang.github.io/bmfaToolkits/reference/select_k_from_sigma.md)),
and the `refit_*()` helpers wrap the “initial fit → select dims → refit”
workflow.

------------------------------------------------------------------------

## Method A: Stack FA (MSFA backend)

“Stack FA” stacks all studies and fits a single factor model.

### A1. Fit

``` r
fit0 <- fit_stack_fa(
  Y_list = Y_list,
  k = 5, # Overspecified
  centering = TRUE,
  scaling = FALSE,
  control = list(nrun = 1000, burn = 200)
)
```

### A2. Post-process

[`postprocess_stack_fa()`](https://mavis-liang.github.io/bmfaToolkits/reference/postprocess_stack_fa.md)
takes `S` so it can return a per-study list for `SigmaMarginal`.

``` r
post0 <- postprocess_stack_fa(fit0, S = S)
names(post0)
```

### A3. Select K and refit

``` r
select_k_stack_fa(post0, cutoff = 0.05)

out <- fit_stack_fa_2step(
  Y_list = Y_list,
  post_fit0 = post0,
  cutoff = 0.05,
  control = list(nrun = 1000, burn = 200)
)
```

### A4. Fit two-step Stack FA in one call

``` r
out <- fit_stack_fa_2step(
  Y_list = Y_list,
  k = 4,# Overspecified
  control = list(nrun = 1000, burn = 200)  
)
```

------------------------------------------------------------------------

## Method B: Ind FA (MSFA backend)

“Ind FA” fits separate factor models per study.

### B1. Fit

You specify `j_s` as either: - a single integer (recycled to all
studies), or - a length-`S` integer vector.

``` r
fit0 <- fit_ind_fa(
  Y_list = Y_list,
  j_s = rep(3, S),
  centering = TRUE,
  scaling = FALSE,
  control = list(nrun = 1000, burn = 200)
)
```

### B2. Post-process, select J_s, refit

``` r
post0 <- postprocess_ind_fa(fit0)
select_js_ind_fa(post0, cutoff = 0.05)

out <- fit_ind_fa_2step(
  Y_list = Y_list,
  post_fit0 = post0,
  cutoff = 0.05,
  control = list(nrun = 1000, burn = 200)
)
```

### B3. Fit two-step Ind FA in one call

``` r
out <- fit_ind_fa_2step(
  Y_list = Y_list,
  j_s = rep(3, S),
  centering = TRUE,
  scaling = FALSE,
  control = list(nrun = 1000, burn = 200)
)
```

------------------------------------------------------------------------

## Method C: PFA (script backend)

PFA scripts are bundled in `inst/extdata/pfa/`. The wrapper: 1) loads
the backend via `load_backend("pfa")`  
2) centers/scales each study (defaults `center=TRUE`, `scale=FALSE`)  
3) stacks into one matrix `X` and builds the study id vector `b`  
4) calls the backend `PFA(X = X, b = b, k = ..., ...)`

### C1. Fit

``` r
fit <- fit_pfa(
  Y_list = Y_list,
  k = 6, # Overspecified
  center = TRUE,
  scale = FALSE,
  nrun = 1000,
  burn = 500
  # other backend parameters via ...
)
```

### C2. Post-process and select K

``` r
post <- postprocess_pfa(fit)
select_k_pfa(post)
```

------------------------------------------------------------------------

## Method D: MOM-SS (BFR.BE backend)

MOM-SS uses the `BFR.BE` backend and requires: - `q` (aka `k`): number
of shared factors - optional covariates `X` (passed to the backend as
`v`)

**In this package wrapper**, if you pass `Y_list`, the wrapper stacks it
into `x`, builds the **membership matrix** `b`, and forwards `q` and
`v`.

### D1. Prepare covariates (optional)

`X` can be either: - a single matrix with `sum(nrow(Y_list))` rows, or -
a list of matrices aligned to `Y_list` (which will be stacked).

Example: two synthetic covariates per subject:

``` r
X_list <- lapply(Y_list, function(Y) {
  n <- nrow(Y)
  cbind(x1 = rnorm(n), x2 = rnorm(n))
})
```

### D2. Fit and post-process

``` r
fit <- fit_momss(
  Y_list = Y_list,
  k = 6,      
  X = X_list, # optional
  scaling = FALSE
)
post <- postprocess_momss(fit)

names(post)
# post$SigmaPhi, post$Psi, post$alpha, post$B, post$SigmaMarginal
```

------------------------------------------------------------------------

## Method E: SUFA (SUFA backend)

In the tutorial workflow, we typically: - center within each study
(`center = TRUE`) - do not rescale variances (`scale = FALSE`) - specify
`qmax` and `nrun`

### E1. Fit and post-process

``` r
#install_backend("sufa")
fit <- fit_sufa(
  Y_list = Y_list,
  k = 6,      # or kmax = 6 depending on your wrapper
  nrun = 1000,
  center = TRUE,
  scale = FALSE
)
post <- postprocess_sufa(fit) 

names(post)
```

------------------------------------------------------------------------

## Method F: BMSFA (MSFA backend)

### F1. Fit

You specify: - `k`: number of shared factors - `j_s`: number of
study-specific factors (length `S`)

You can control MCMC via `control = list(nrun = ..., burn = ...)`.

``` r
fit0 <- fit_bmsfa(
  Y_list = Y_list,
  k = 5, # Overspecified
  j_s = rep(2, S),
  centering = TRUE,
  scaling = FALSE,
  control = list(nrun = 2000, burn = 1000)
)
```

### F2. Post-process

``` r
post0 <- postprocess_bmsfa(fit0)

names(post0)
```

### F3. Select dimensions and refit

``` r
select_k_bmsfa(post0, cutoff = 0.05)
select_js_bmsfa(post0, cutoff = 0.05)

# refit
out <- fit_bmsfa_2step(
  Y_list = Y_list,
  post_fit0 = post0,
  cutoff = 0.05,
  control = list(nrun = 1000, burn = 200)
)
# Or `fit_bmsfa` with the selected k and j_s, and `postprocess` again.
```

### F4. Fit two-step BMSFA in one call

``` r
out <- fit_bmsfa_2step(
  Y_list = Y_list,
  k = 5, # Overspecified
  j_s = rep(2, S),
  centering = TRUE,
  scaling = FALSE,
  control = list(nrun = 1000, burn = 200)
)
```

------------------------------------------------------------------------

## Method G: CAVI (VI-MSFA backend)

CAVI is a variational inference algorithm for multi-study factor
analysis implemented in the GitHub repo `blhansen/VI-MSFA` (R package
name: `VIMSFA`).

**What you need** - `Y_list`: a list of S matrices, each `n_s × P`,
where all studies share the same variables (same `P` and column
order). - `k` (`K` in VI-MSFA): number of common factors. - `j_s` (`J_s`
in VI-MSFA): study-specific factor counts (either a scalar or a length-S
vector). - Preprocessing: this project typically **centers only**
(`center=TRUE, scale=FALSE`) to match other methods.

``` r
# install the backend once (downloads / installs VIMSFA from GitHub)
bmfaToolkits::install_backend("cavi")

Y_list <- toy$Y_list

fit_cavi <- bmfaToolkits::fit_cavi(
  Y_list = Y_list,
  k      = 3,                 # common factors
  j_s    = 2,                 # study-specific factors (scalar -> recycled to length S)
  centering = TRUE,
  scaling   = FALSE
)

post_cavi <- bmfaToolkits::postprocess_cavi(fit_cavi)
str(post_cavi, max.level = 1)
```

**What
[`postprocess_cavi()`](https://mavis-liang.github.io/bmfaToolkits/reference/post_CAVI.md)
does behind the scenes**

[`VIMSFA::cavi_msfa()`](https://rdrr.io/pkg/VIMSFA/man/cavi_msfa.html)
returns posterior means: - `mean_phi` (P × K): common loading matrix Φ -
`mean_lambda_s` (list of P × J_s): study-specific loading matrices Λ_s -
`mean_psi_s` (list of length-P): residual variances ψ_s

[`bmfaToolkits::post_CAVI()`](https://mavis-liang.github.io/bmfaToolkits/reference/post_CAVI.md)/[`postprocess_cavi()`](https://mavis-liang.github.io/bmfaToolkits/reference/post_CAVI.md)
converts these into the standardized outputs used throughout this
package: - `SigmaPhi = Φ Φᵀ` - `SigmaLambdaList[[s]] = Λ_s Λ_sᵀ` -
`SigmaMarginal[[s]] = SigmaPhi + SigmaLambdaList[[s]] + diag(ψ_s)`

#### Automatic factor selection and refit (2-step)

Like BMSFA in this package, CAVI can be run in a **2-step** workflow:

1.  Fit an initial (often over-specified) model with `(k, j_s)`.
2.  Post-process to get `SigmaPhi` and each `SigmaLambda_s`.
3.  Select `K` and `J_s` using the eigen-proportion rule (`cutoff`).
4.  Refit with the selected factor counts.

``` r
# Step 1: initial fit (choose generous k / j_s)
fit0  <- bmfaToolkits::fit_cavi(Y_list, k = 8, j_s = 6, centering = TRUE, scaling = FALSE)
post0 <- bmfaToolkits::postprocess_cavi(fit0)

# Step 2: refit using selected K and J_s
post2 <- bmfaToolkits::refit_cavi(Y_list, post_fit0 = post0, cutoff = 0.05,
                                  centering = TRUE, scaling = FALSE)

# Or: do both in one call
post2_alt <- bmfaToolkits::fit_cavi_2step(Y_list, k = 8, j_s = 6, cutoff = 0.05,
                                         centering = TRUE, scaling = FALSE)
```

## Method H. BLAST (script backend)

BLAST is implemented as a research-code repository (not a
CRAN/Bioconductor package) in: `maurilorenzo/BLAST`.

In this package we treat BLAST as a **script backend** (like
PFA/Tetris):

- When developing from source: put the BLAST scripts under
  `inst/extdata/blast/`.
- After installation: they are available under
  `system.file("extdata","blast", package="bmfaToolkits")`.

[`fit_blast()`](https://mavis-liang.github.io/bmfaToolkits/reference/fit_blast.md)
loads the scripts via `load_backend("blast")`.

**What you need** - `Y_list`: a list of S matrices, each `n_s × P` (same
P across studies). - `k`: number of common factors. - `q_s`:
study-specific factor counts (scalar or length S). - Important BLAST
args in this project: - `n_MC`: Monte Carlo iterations (e.g., 10000). -
`sample_outer_product`: if `FALSE`, BLAST requires `subsample_index` (we
default to `1:min(100, P)`).

``` r
Y_list <- toy$Y_list

fit_blast <- bmfaToolkits::fit_blast(
  Y_list = Y_list,
  k      = 3,
  q_s    = 1,
  n_MC   = 10000,
  center = TRUE,
  scale  = FALSE,
  sample_outer_product = FALSE
  # subsample_index defaults to 1:min(100, P) when sample_outer_product = FALSE
)

post_blast <- bmfaToolkits::postprocess_blast(fit_blast)
str(post_blast, max.level = 1)
```

**What
[`postprocess_blast()`](https://mavis-liang.github.io/bmfaToolkits/reference/post_BLAST.md)
does behind the scenes**

BLAST’s fit object is a list whose fields may vary depending on
settings.
[`bmfaToolkits::post_BLAST()`](https://mavis-liang.github.io/bmfaToolkits/reference/post_BLAST.md)
attempts to extract: - `Lambda_mean` (P × K) as `Phi` -
`Lambda_outer_mean` (P × P) as `SigmaPhi` if present; otherwise compute
`Phi Phiᵀ` - `Gammas_mean` (P × q × S) and/or `Gammas_outer_mean` (P × P
× S) for study-specific components - `Sigma_2s_samples` (N_mc × P) for
residual variances (ψ); if absent, residual variances are returned as
`NA`

It then returns the standardized pieces (`Phi`, `SigmaPhi`,
`LambdaList`, `SigmaLambdaList`, `PsiList`, `SigmaMarginal`) so BLAST
can plug into the same downstream code as other methods.

------------------------------------------------------------------------

## Method I: Tetris (script backend)

Tetris scripts are bundled in `inst/extdata/tetris/`.

The “3-step” workflow from the tutorial is:

1.  initial run: `tetris(..., fixed_bigT = FALSE)`  
2.  select `big_T`: `choose.A(fit, alpha_IBP = alpha, S = S)`  
3.  fixed run: `tetris(..., fixed_bigT = TRUE, bigT = big_T)`

### I1. Manual 3-step pipeline

``` r
fit1 <- fit_tetris(Y_list, alpha = "auto", beta = 1, fixed_bigT = FALSE, nrun = 500, burn = 200, nprint = 100) # Reduced iterations a faster runtime
big_T <- select_T_tetris(fit1)  # passes alpha_IBP and S from fit1 metadata by default
fit2 <- fit_tetris(Y_list, alpha = fit1$meta$alpha, beta = 1, fixed_bigT = TRUE, bigT = big_T,
                   nrun = 500, burn = 200, nprint = 100)
out <- postprocess_tetris(fit2)
```

### I2. One-shot pipeline (recommended)

``` r
out <- fit_tetris_2step(
  Y_list = Y_list,
  alpha = "auto",   # ceiling(1.25 * S)
  beta = 1,
  nrun = 500,
  burn = 200,
  nprint = 100
)
```

------------------------------------------------------------------------

## Invoke any method: fit_integrative_fa()

- `fit_integrative_fa(method = "bmsfa"|"momss"|..., ...)`
- `postprocess_integrative_fa(method = ..., fit)`

Example:

``` r
fit <- fit_integrative_fa("bmsfa", Y_list = Y_list, k = 4, j_s = rep(2, S))
post <- postprocess_integrative_fa("bmsfa", fit)
```

## Visuallization

This function give a nice heatmap plot for the loadings

``` r
plot_single_loadings(mat = out$Phi, fill_limits = c(-2, 2)) # Displays grey if values are out of range
```

------------------------------------------------------------------------

### Practical tips

- Start with **short runs** (`nrun`/`burn`) to validate your pipeline,
  then increase.
- Always check that each study has the same `P` (number of columns).
- For methods with rotational ambiguity, rely on the package
  post-processing (e.g., OP alignment) before comparing loadings.

### Session info

``` r
sessionInfo()
```

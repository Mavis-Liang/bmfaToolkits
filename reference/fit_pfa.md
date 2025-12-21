# Fit PFA (Perturbed Factor Analysis; script backend)

PFA in the tutorial is provided as scripts (R + C++ via Rcpp). This
wrapper expects the PFA scripts to be bundled under `inst/extdata/pfa/`
in this package source. At runtime, the backend is loaded via
[`load_backend()`](https://mavis-liang.github.io/bifaToolkits/reference/load_backend.md),
which sources all `.R` files in that directory and compiles `PFA.cpp`
(if present) via
[`Rcpp::sourceCpp()`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html).

## Usage

``` r
fit_pfa(
  Y_list = NULL,
  k = NULL,
  center = TRUE,
  scale = FALSE,
  nrun = 3000,
  burn = 2000,
  ...
)
```

## Arguments

- Y_list:

  List of study matrices.

- k:

  Positive integer number of shared factors.

- center:

  Logical; whether to center each variable within each study.

- scale:

  Logical; whether to scale each variable within each study.

- nrun:

  Number of MCMC iterations.

- burn:

  Number of burn-in iterations.

- ...:

  Passed to `PFA()` in the backend script. User can set hyper-parameters
  and numbers of iterations here.

## Value

A `bifa_fit` object.

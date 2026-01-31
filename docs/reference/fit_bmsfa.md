# Fit BMSFA via MSFA::sp_msfa

Fit BMSFA via MSFA::sp_msfa

## Usage

``` r
fit_bmsfa(
  Y_list,
  k,
  j_s,
  outputlevel = 1,
  scaling = FALSE,
  centering = TRUE,
  control = list(nrun = 3000, burn = 2000),
  ...
)
```

## Arguments

- Y_list:

  List of study matrices.

- k:

  Number of common factors.

- j_s:

  Study-specific factors (length S).

- outputlevel, scaling, centering, control:

  Passed to
  [`MSFA::sp_msfa()`](https://rdrr.io/pkg/MSFA/man/sp_msfa.html).

- ...:

  Passed to
  [`MSFA::sp_msfa()`](https://rdrr.io/pkg/MSFA/man/sp_msfa.html).

## Value

A `bifa_fit` object.

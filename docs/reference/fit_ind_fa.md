# Fit independent factor analysis (Ind FA) via MSFA::sp_fa (per study)

Fit independent factor analysis (Ind FA) via MSFA::sp_fa (per study)

## Usage

``` r
fit_ind_fa(
  Y_list,
  j_s,
  scaling = FALSE,
  centering = TRUE,
  control = list(nrun = 3000, burn = 2000),
  ...
)
```

## Arguments

- Y_list:

  List of study matrices.

- j_s:

  Either a single integer or a length-S integer vector of factors.

- scaling, centering, control:

  Passed to [`MSFA::sp_fa()`](https://rdrr.io/pkg/MSFA/man/sp_fa.html).

- ...:

  Passed to [`MSFA::sp_fa()`](https://rdrr.io/pkg/MSFA/man/sp_fa.html).

## Value

A `bifa_fit` object whose `$fit` is a list of fits (one per study).

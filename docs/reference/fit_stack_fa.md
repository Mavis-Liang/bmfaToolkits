# Fit stacked factor analysis (Stack FA) via MSFA::sp_fa

This is the "Stack FA" baseline used in the tutorial: stack all studies
and fit a single factor model.

## Usage

``` r
fit_stack_fa(
  Y_list,
  k,
  scaling = FALSE,
  centering = TRUE,
  control = list(nrun = 3000, burn = 2000),
  ...
)
```

## Arguments

- Y_list:

  List of study matrices (N_s x P).

- k:

  Number of factors.

- scaling, centering:

  Passed to [`MSFA::sp_fa()`](https://rdrr.io/pkg/MSFA/man/sp_fa.html).

- control:

  Passed to [`MSFA::sp_fa()`](https://rdrr.io/pkg/MSFA/man/sp_fa.html).

- ...:

  Passed to [`MSFA::sp_fa()`](https://rdrr.io/pkg/MSFA/man/sp_fa.html).

## Value

A `bifa_fit` object.

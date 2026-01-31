# Fit SUFA

Thin wrapper around
[`SUFA::fit_SUFA()`](https://rdrr.io/pkg/SUFA/man/fit_SUFA.html)

## Usage

``` r
fit_sufa(Y_list, k = NULL, nrun = 3000, center = TRUE, scale = FALSE, ...)
```

## Arguments

- Y_list:

  List of study matrices.

- k:

  Positive integer for maximum number of factors shared (`qmax` in
  SUFA).

- nrun:

  Number of MCMC iterations.

- center:

  Whether to center each study matrix (default TRUE).

- scale:

  Whether to scale each study matrix (default FALSE).

- ...:

  Passed to
  [`SUFA::fit_SUFA()`](https://rdrr.io/pkg/SUFA/man/fit_SUFA.html).

## Value

A `bifa_fit` object.

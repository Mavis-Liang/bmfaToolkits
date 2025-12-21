# Fit MOM-SS (MOM-SS / BFR.BE backend)

This is a thin wrapper around `BFR.BE::BFR.BE.EM.CV()` (or similar entry
points). Because the upstream API may evolve, we keep this wrapper
flexible:

- If you pass `...` with `X` and `M` (the membership matrix), they are
  forwarded.

- If you pass `Y_list`, we will stack it and construct `X` and `M`
  accordingly.

## Usage

``` r
fit_momss(
  Y_list,
  k,
  X = NULL,
  scaling = FALSE,
  center_within_study = TRUE,
  folds = 10,
  ...
)
```

## Arguments

- Y_list:

  Optional list of study matrices. If provided, we create a stacked
  matrix.

- k:

  Positive integer number of shared factors.

- X:

  Optional covariate matrix or list of covariate matrices aligned to
  `Y_list`.

- scaling:

  Logical; whether to scale the data.

- center_within_study:

  Logical; whether to center each study's data before fitting.
  Recommended unless you have strong reasons not to.

- folds:

  Number of folds for cross-validation within MOM-SS.

- ...:

  Parameters that can be pass directly to `BFR.BE::BFR.BE.EM.CV()`.

## Value

A `bifa_fit` object.

# Fit VI-MSFA CAVI (variational MSFA)

Wrapper around `VIMSFA::cavi_msfa()` that standardizes preprocessing and
return type.

## Usage

``` r
fit_cavi(
  Y_list,
  k,
  j_s,
  scaling = FALSE,
  centering = TRUE,
  center = NULL,
  scale = NULL,
  ...
)
```

## Arguments

- Y_list:

  List of S study matrices (each `n_s x P`).

- k:

  Integer; number of common factors (`K` in VI-MSFA).

- j_s:

  Integer; study-specific factor counts (`J_s` in VI-MSFA), length 1 or
  length S.

- scaling:

  Logical; scale each study matrix before fitting (default `FALSE`).

- centering:

  Logical; center each study matrix before fitting (default `TRUE`).

- center, scale:

  Deprecated aliases for `centering` / `scaling` (kept for backward
  compatibility).

- ...:

  Passed to `VIMSFA::cavi_msfa()`.

## Value

A `bifa_fit` object with backend `"cavi"`.

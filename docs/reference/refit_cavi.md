# Refit CAVI after selecting the number of factors

Implements a 2-step workflow: fit an initial CAVI model, post-process to
obtain covariance components, select factor counts using an
eigen-proportion rule, and refit CAVI with the selected `K` and `J_s`.

## Usage

``` r
refit_cavi(Y_list, fit0 = NULL, post_fit0 = NULL, cutoff = 0.05, ...)
```

## Arguments

- Y_list:

  List of S study matrices (each `n_s x P`).

- fit0:

  Optional initial fit from
  [`fit_cavi()`](https://mavis-liang.github.io/bmfaToolkits/reference/fit_cavi.md).

- post_fit0:

  Optional post-processed initial fit from
  [`postprocess_cavi()`](https://mavis-liang.github.io/bmfaToolkits/reference/post_CAVI.md).
  If both `fit0` and `post_fit0` are supplied, `post_fit0` is used.

- cutoff:

  Eigen-proportion cutoff passed to
  [`select_k_from_sigma()`](https://mavis-liang.github.io/bmfaToolkits/reference/select_k_from_sigma.md).

- ...:

  Passed to
  [`fit_cavi()`](https://mavis-liang.github.io/bmfaToolkits/reference/fit_cavi.md)
  for the refit.

## Value

Post-processed output of the refit (same structure as
[`postprocess_cavi()`](https://mavis-liang.github.io/bmfaToolkits/reference/post_CAVI.md)).

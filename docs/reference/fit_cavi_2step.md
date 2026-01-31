# Fit CAVI with automatic factor selection (2-step)

Convenience wrapper that (i) fits an initial CAVI model, (ii) selects
factor counts from covariance components, and (iii) refits with the
selected counts.

## Usage

``` r
fit_cavi_2step(Y_list, post_fit0 = NULL, k = NULL, j_s = NULL, cutoff = 0.05, ...)
```

## Arguments

- Y_list:

  List of S study matrices (each `n_s x P`).

- post_fit0:

  Optional post-processed initial fit (from
  [`postprocess_cavi()`](https://mavis-liang.github.io/bmfaToolkits/reference/post_CAVI.md)).
  If provided, the initial fit stage is skipped.

- k:

  Initial number of common factors (required if `post_fit0` is `NULL`).

- j_s:

  Initial study-specific factor counts (required if `post_fit0` is
  `NULL`).

- cutoff:

  Eigen-proportion cutoff passed to
  [`select_k_from_sigma()`](https://mavis-liang.github.io/bmfaToolkits/reference/select_k_from_sigma.md).

- ...:

  Passed to
  [`fit_cavi()`](https://mavis-liang.github.io/bmfaToolkits/reference/fit_cavi.md)
  for both stages.

## Value

Post-processed output of the refit (same structure as
[`postprocess_cavi()`](https://mavis-liang.github.io/bmfaToolkits/reference/post_CAVI.md)).

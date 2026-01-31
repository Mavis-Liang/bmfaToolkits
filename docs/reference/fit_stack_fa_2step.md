# Fit Stack FA all-at-once If `post_fit0` is NULL, fit an initial StackFA with given `k`, then automatically select the number of factors and refit.

Fit Stack FA all-at-once If `post_fit0` is NULL, fit an initial StackFA
with given `k`, then automatically select the number of factors and
refit.

## Usage

``` r
fit_stack_fa_2step(Y_list, post_fit0 = NULL, k = NULL, cutoff = 0.05, ...)
```

## Arguments

- Y_list:

  Data list.

- post_fit0:

  Initial post-processed result (or NULL to compute internally).

- k:

  Initial number of factors (used only if `post_fit0` is NULL).

- cutoff:

  Eigen proportion cutoff.

- ...:

  Passed to
  [`fit_stack_fa()`](https://mavis-liang.github.io/bmfaToolkits/reference/fit_stack_fa.md).

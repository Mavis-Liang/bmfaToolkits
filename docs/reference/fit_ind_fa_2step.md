# Fit Ind FA all at once If `post_fit0` is NULL, fit an initial IndFA with given `j_s`, then automatically select the number of factors and refit.

Fit Ind FA all at once If `post_fit0` is NULL, fit an initial IndFA with
given `j_s`, then automatically select the number of factors and refit.

## Usage

``` r
fit_ind_fa_2step(Y_list, post_fit0 = NULL, j_s = NULL, cutoff = 0.05, ...)
```

## Arguments

- Y_list:

  Data list.

- post_fit0:

  Initial post-processing result (or NULL to compute internally).

- j_s:

  Initial number of factors (or NULL to use all).

- cutoff:

  Eigen proportion cutoff.

- ...:

  Passed to
  [`fit_ind_fa()`](https://mavis-liang.github.io/bmfaToolkits/reference/fit_ind_fa.md).

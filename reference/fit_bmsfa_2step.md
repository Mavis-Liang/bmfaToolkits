# Fit BMSFA all-at-once If `post_fit0` is NULL, fit an initial BMSFA with given `k` and `j_s`, then automatically select the number of factors and refit. Otherwise, use the `post_fit0` then select the number of factors and refit.

Fit BMSFA all-at-once If `post_fit0` is NULL, fit an initial BMSFA with
given `k` and `j_s`, then automatically select the number of factors and
refit. Otherwise, use the `post_fit0` then select the number of factors
and refit.

## Usage

``` r
fit_bmsfa_2step(
  Y_list,
  post_fit0 = NULL,
  k = NULL,
  j_s = NULL,
  cutoff = 0.05,
  ...
)
```

## Arguments

- Y_list:

  Data list.

- post_fit0:

  Post-processed initial fit (or NULL).

- k:

  Initial number of common factors (if `post_fit0` is NULL).

- j_s:

  Initial study-specific factors (if `post_fit0` is NULL).

- cutoff:

  Eigen proportion cutoff.

- ...:

  Passed to
  [`fit_bmsfa()`](https://mavis-liang.github.io/bifaToolkits/reference/fit_bmsfa.md).

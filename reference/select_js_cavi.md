# Select J_s for CAVI from study-specific covariances

Applies
[`select_k_from_sigma()`](https://mavis-liang.github.io/bmfaToolkits/reference/select_k_from_sigma.md)
to each element of `post$SigmaLambdaList`.

## Usage

``` r
select_js_cavi(post, cutoff = 0.05)
```

## Arguments

- post:

  Output of
  [`postprocess_cavi()`](https://mavis-liang.github.io/bmfaToolkits/reference/post_CAVI.md)
  (or
  [`post_CAVI()`](https://mavis-liang.github.io/bmfaToolkits/reference/post_CAVI.md)).

- cutoff:

  Eigen-proportion cutoff.

## Value

A list of length S with selected `J_s`.

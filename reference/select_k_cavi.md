# Select K for CAVI from the common covariance

Applies
[`select_k_from_sigma()`](https://mavis-liang.github.io/bmfaToolkits/reference/select_k_from_sigma.md)
to `post$SigmaPhi`.

## Usage

``` r
select_k_cavi(post, cutoff = 0.05)
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

An integer selected `K`.

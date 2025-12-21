# Select K for PFA via nontrivial-column rule

In the nutrition case study, a factor is considered "kept" if at least
one loading exceeds a small threshold in absolute value.

## Usage

``` r
select_k_pfa(post, threshold = 0.001)
```

## Arguments

- post:

  Output of
  [`postprocess_pfa()`](https://mavis-liang.github.io/bifaToolkits/reference/postprocess_pfa.md).

- threshold:

  Threshold for considering a column nontrivial (default 1e-3).

## Value

Integer K.

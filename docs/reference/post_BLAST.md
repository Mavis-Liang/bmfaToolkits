# Post-process BLAST output into a standard format

Extracts posterior mean loadings/covariances from BLAST output and
constructs common/study-specific covariance pieces.

## Usage

``` r
post_BLAST(fit) 
postprocess_blast(fit, ...)
```

## Arguments

- fit:

  A `bifa_fit` returned by
  [`fit_blast()`](https://mavis-liang.github.io/bmfaToolkits/reference/fit_blast.md)
  or a raw BLAST fit list.

- ...:

  Ignored (kept for interface consistency).

## Value

A list with covariance components and marginal covariance per study.

## See also

[`fit_blast`](https://mavis-liang.github.io/bmfaToolkits/reference/fit_blast.md)

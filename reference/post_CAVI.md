# Post-process CAVI output into a standard format

Extracts posterior means from VI-MSFA CAVI output and constructs
common/study-specific covariance pieces.

## Usage

``` r
post_CAVI(fit)
postprocess_cavi(fit, ...)
```

## Arguments

- fit:

  A `bifa_fit` returned by
  [`fit_cavi()`](https://mavis-liang.github.io/bmfaToolkits/reference/fit_cavi.md)
  or a raw CAVI fit list.

- ...:

  Ignored (kept for interface consistency).

## Value

A list with covariance components and marginal covariance per study.

## See also

[`fit_cavi`](https://mavis-liang.github.io/bmfaToolkits/reference/fit_cavi.md)

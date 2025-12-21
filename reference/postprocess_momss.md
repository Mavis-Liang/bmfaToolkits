# Post-process MOM-SS output

The tutorial uses `fit$M` as the loading matrix and `fit$Sigma` as the
list of study-specific marginal covariances.

## Usage

``` r
postprocess_momss(fit)
```

## Arguments

- fit:

  A `bifa_fit` from
  [`fit_momss()`](https://mavis-liang.github.io/bifaToolkits/reference/fit_momss.md)
  or a raw MOM-SS fit.

## Value

A list with `Phi`, `SigmaPhi`, `SigmaMarginal`, `PsiList`.

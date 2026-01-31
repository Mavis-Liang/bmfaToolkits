# Post-process Ind FA output (nutrition tutorial)

Post-process Ind FA output (nutrition tutorial)

## Usage

``` r
postprocess_ind_fa(fit)
```

## Arguments

- fit:

  A `bifa_fit` from
  [`fit_ind_fa()`](https://mavis-liang.github.io/bmfaToolkits/reference/fit_ind_fa.md)
  or a raw list of `MSFA::sp_fa()` fits.

## Value

A list with `LambdaList`, `SigmaLambdaList`, `Psi` (list),
`SigmaMarginal` (list).

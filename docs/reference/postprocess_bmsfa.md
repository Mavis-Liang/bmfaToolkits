# Post-process BMSFA output

OP on Phi and each Lambda_s, form SigmaPhi, SigmaLambdaList, estimate
Psi_s by averaging `fit$psi[[s]]` over iterations, and build marginal
covariance as SigmaPhi + SigmaLambda_s + diag(Psi_s).

## Usage

``` r
postprocess_bmsfa(fit)
```

## Arguments

- fit:

  A `bifa_fit` from
  [`fit_bmsfa()`](https://mavis-liang.github.io/bmfaToolkits/reference/fit_bmsfa.md)
  or a raw
  [`MSFA::sp_msfa()`](https://rdrr.io/pkg/MSFA/man/sp_msfa.html) fit.

## Value

A list with `Phi`, `SigmaPhi`, `LambdaList`, `SigmaLambdaList`,
`PsiList`, `SigmaMarginal`.

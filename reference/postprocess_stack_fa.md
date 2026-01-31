# Post-process Stack FA output (nutrition tutorial)

Implements the post-processing in the nutrition case study: apply
orthogonal Procrustes (OP) to posterior loadings, form SigmaPhi = Phi
Phi^T, estimate marginal covariance by averaging Sigma draws, and
estimate Psi by averaging Sigma - Lambda Lambda^T over MCMC iterations.

## Usage

``` r
postprocess_stack_fa(fit, S)
```

## Arguments

- fit:

  A `bifa_fit` object from
  [`fit_stack_fa()`](https://mavis-liang.github.io/bmfaToolkits/reference/fit_stack_fa.md)
  or a raw `MSFA::sp_fa()` fit.

- S:

  Number of studies (used only to replicate marginal covariance into a
  list).

## Value

A list with `Phi`, `SigmaPhi`, `Psi`, `SigmaMarginal`.

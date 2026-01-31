# Post-process Tetris output

The tutorial extracts study-specific loadings via `getLambda(fit)` from
the Tetris script and builds marginal covariances as Lambda Lambda^T +
Psi.

## Usage

``` r
postprocess_tetris(fit, num_samps = NULL)
```

## Arguments

- fit:

  A `bifa_fit` from
  [`fit_tetris()`](https://mavis-liang.github.io/bmfaToolkits/reference/fit_tetris.md)
  or a raw Tetris fit.

- num_samps:

  Number of posterior samples to use. If NULL, uses min(2000,
  available).

## Value

A list with `Phi`, `SigmaPhi`, `LambdaList`, `SigmaLambdaList`,
`PsiList`, `SigmaMarginal`.

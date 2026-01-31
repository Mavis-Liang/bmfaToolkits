# Select big_T (A matrix) for Tetris

Uses choose.A() on a Tetris fit to select the A matrix.

## Usage

``` r
select_T_tetris(fit, alpha_IBP = NULL, S = NULL, ...)
```

## Arguments

- fit:

  A `bifa_fit` returned by
  [`fit_tetris()`](https://mavis-liang.github.io/bmfaToolkits/reference/fit_tetris.md)
  (typically the non-fixed run).

- alpha_IBP:

  Alpha used by choose.A(). If NULL, uses fit metadata alpha.

- S:

  Number of studies. If NULL, uses fit metadata S.

- ...:

  Passed to choose.A().

## Value

The selected A matrix (big_T).

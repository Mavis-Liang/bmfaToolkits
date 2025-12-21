# Fit Tetris (single run)

Runs one Tetris fit. Optionally uses a fixed A matrix (selected by
choose.A()).

## Usage

``` r
fit_tetris(
  Y_list,
  alpha = c("auto"),
  beta = 1,
  fixed_bigT = FALSE,
  bigT = NULL,
  nrun = 3000,
  burn = 2000,
  nprint = 200,
  ...
)
```

## Arguments

- Y_list:

  List of study matrices.

- alpha:

  Either numeric alpha, or "auto" to use ceiling(1.25 \* S).

- beta:

  Beta hyperparameter for Tetris.

- fixed_bigT:

  Logical; whether to run with fixed big T.

- bigT:

  Fixed T matrix used when fixed=TRUE.

- nrun:

  Number of iterations.

- burn:

  Burn-in iterations.

- nprint:

  Print frequency.

- ...:

  Passed to backend tetris().

## Value

A `bifa_fit` object.

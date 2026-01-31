# Run the full Tetris pipeline: initial fit -\> choose T -\> fixed refit

Mirrors the common workflow:

1.  tetris(..., fixed_bigT=FALSE)

2.  big_T \<- choose.A(fit, ...)

3.  tetris(..., fixed_bigT=TRUE, A_fixed=big_T)

## Usage

``` r
fit_tetris_2step(
  Y_list,
  alpha = "auto",
  beta = 1,
  nrun = 3000,
  burn = 2000,
  nprint = 200,
  chooseA_args = list(),
  ...
)
```

## Arguments

- Y_list:

  List of study matrices.

- alpha:

  Either numeric alpha, or "auto" to use ceiling(1.25 \* S).

- beta:

  Beta hyperparameter.

- nrun, burn, nprint:

  MCMC controls for both runs.

- chooseA_args:

  List of extra arguments passed to choose.A().

- ...:

  Extra args passed to tetris() in both runs.

## Value

A `bifa_fit` (fixed refit) with `stage1` and `big_T` attached.

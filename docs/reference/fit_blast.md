# Fit BLAST (script backend)

Wrapper that loads BLAST scripts from the package's bundled extdata (or
the source tree under `inst/extdata/blast/`) and calls BLAST's
`fit_blast()`.

## Usage

``` r
fit_blast(Y_list, k, q_s, n_MC = 10000, center = TRUE, scale = FALSE, subsample_index = NULL, sample_outer_product = FALSE, ...)
```

## Arguments

- Y_list:

  List of S study matrices (each n_s x P).

- k:

  Number of common factors.

- q_s:

  Study-specific factor counts, length 1 or length S.

- n_MC:

  Number of Monte Carlo iterations used by BLAST.

- center:

  Logical; center each study matrix before fitting.

- scale:

  Logical; scale each study matrix before fitting.

- subsample_index:

  Integer indices used by BLAST when `sample_outer_product = FALSE`.

- sample_outer_product:

  Logical; forwarded to BLAST.

- ...:

  Passed to the BLAST implementation of `fit_blast()`.

## Value

A `bifa_fit` object with backend `"blast"`.

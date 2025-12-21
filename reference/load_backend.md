# Load a bundled script backend into a private environment

For PFA/Tetris, this sources all `.R` files under
`inst/extdata/<backend>/`. For PFA, it also compiles `PFA.cpp` (if
present) via
[`Rcpp::sourceCpp()`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html).

## Usage

``` r
load_backend(backend)
```

## Arguments

- backend:

  "pfa" or "tetris".

## Value

An environment containing the sourced functions.

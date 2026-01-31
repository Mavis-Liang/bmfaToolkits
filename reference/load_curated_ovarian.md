# Load curated ovarian data (Bioconductor)

This does not ship the data inside `bmfaToolkits`. Instead, it provides
a thin wrapper that loads from the Bioconductor package
`curatedOvarianData`.

## Usage

``` r
load_curated_ovarian(...)
```

## Arguments

- ...:

  Passed to the underlying data accessor.

## Value

An object from `curatedOvarianData`.

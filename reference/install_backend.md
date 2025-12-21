# Install optional backends used in the tutorial

This is a convenience wrapper for non-CRAN dependencies.

## Usage

``` r
install_backend(backend, ...)
```

## Arguments

- backend:

  One of: "bmsfa", "momss", "sufa", "curated_ovarian".

- ...:

  Passed to
  [`remotes::install_github()`](https://remotes.r-lib.org/reference/install_github.html)
  or
  [`BiocManager::install()`](https://bioconductor.github.io/BiocManager/reference/install.html)
  when applicable.

## Details

- For GitHub packages, we call
  [`remotes::install_github()`](https://remotes.r-lib.org/reference/install_github.html).

- For Bioconductor packages, we call
  [`BiocManager::install()`](https://bioconductor.github.io/BiocManager/reference/install.html).

- For script-based methods (PFA / Tetris), scripts are bundled under
  `inst/extdata/`.

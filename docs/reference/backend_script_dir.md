# Get the installed directory for a script-based backend

Script backends (PFA/Tetris) are bundled inside this package under
`inst/extdata/`. Use this helper to locate them after installation.

## Usage

``` r
backend_script_dir(backend)
```

## Arguments

- backend:

  "pfa" or "tetris".

## Value

Full path to the backend directory inside the installed package.

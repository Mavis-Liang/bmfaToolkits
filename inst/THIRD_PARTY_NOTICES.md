# Third-party notices (scripts/backends)

This package is intended as a lightweight wrapper layer for a tutorial paper.
Some methods referenced in the tutorial are implemented externally.

## GitHub / Bioconductor packages (not redistributed)

These are installed by users and are not bundled inside this package:
- MSFA (BMSFA / Stack FA / Ind FA) via `remotes::install_github("Mavis-Liang/MSFA")`
- BFR.BE (MOM-SS) via `remotes::install_github("AleAviP/BFR.BE")`
- SUFA via `remotes::install_github("noirritchandra/SUFA", build_vignettes = FALSE)`
- curatedOvarianData (Bioconductor) via `BiocManager::install("curatedOvarianData")`

## Script backends (downloaded into user cache)

For convenience and to avoid redistributing third-party code, `bifa::install_backend()`
can download scripts used in the tutorial into a user cache directory (see `?backend_path`).

- PFA script URL used by default:
  - https://github.com/Mavis-Liang/Bayesian_integrative_FA_tutorial/raw/refs/heads/main/FBPFA-PFA.R

- Tetris script URL used by default:
  - https://github.com/Mavis-Liang/Bayesian_integrative_FA_tutorial/raw/refs/heads/main/Tetris.R

If you plan to redistribute those scripts inside this package, you should:
1) verify the original repositories' licensing terms,
2) include their LICENSE text under `inst/`,
3) cite them in `inst/CITATION` and in your vignette(s),
4) record the exact commit hash you vendored.

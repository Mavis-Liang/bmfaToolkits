# bifa

A lightweight R package skeleton to support the *Bayesian integrative factor analysis* tutorial/paper.

See its demo usage in the accompanying vignette: `vignette("bifaToolkits")`, or better, 
at https://mavis-liang.github.io/Bayesian_integrative_FA_tutorial_book/.


## Install (local)

```r
# From a folder that contains this package:
# install.packages("devtools")
devtools::install(".", dependencies = TRUE, build_vignettes = FALSE),
# Or install from GitHub:
devtools::install_github("your_github_username/bifaToolkits", dependencies = TRUE, build_vignettes = FALSE)
# The vignette may fail to build if optional backends are not installed.
```

## Quick start

```r
library(bifaToolkits)

# See which optional backends are available on your machine
available_backends()

# Example: quick fit BMSFA (requires MSFA)
data <- bifaToolkits::load_toy() #n_s=c(52,58,40), P=20
fit  <- fit_bmsfa(data$Y_list, k = 4, j_s = rep(2, length(data$Y_list)))
post <- postprocess_bmsfa(fit)
select_k_bmsfa(post)
```

## Optional backends

For GitHub / Bioconductor backends, see `?install_backend`.

### Script backends (PFA / Tetris)

This package expects the PFA/Tetris scripts to be bundled under:

- `inst/extdata/pfa/` (one or more `.R` files; optional `PFA.cpp`)
- `inst/extdata/tetris/` (one or more `.R` files)

At runtime, the package sources all `.R` files in those folders via `load_backend()`.

# bifaToolkits

<!-- badges: start -->
[![R-CMD-check](https://github.com/Mavis-Liang/bifaToolkits/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Mavis-Liang/bifaToolkits/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Toolkits and reproducible workflows for Bayesian integrative factor analysis methods, with a lightweight “works out of the box” interface and optional backends for additional models.

## Install

### Development version (GitHub)

```r
install.packages("remotes")
remotes::install_github("Mavis-Liang/bifaToolkits")
```

## Optional backends

Some methods depend on optional packages and/or vendored scripts. Use available_backends() to see what is usable on your system, and install extras as needed.

```
library(bifaToolkits)
available_backends()
```

If a backend is missing, the corresponding fit_*() function will give a clear message telling you what to install.

Detail usage see https://mavis-liang.github.io/bifaToolkits/. Navigate to the "Articles" tab for vignettes.

## Quick start (toy data)

```{r}
library(bifaToolkits)

toy <- readRDS(system.file("extdata", "toy.rds", package = "bifaToolkits", mustWork = TRUE))
Y_list <- toy$Y_list

# Example: MOM-SS (requires BFR.BE + its runtime deps)
fit <- fit_momss(Y_list, k = 6, scaling = FALSE)
pp  <- postprocess_momss(fit)

# Plot an estimated shared covariance (example)
plot_single_loadings(pp$Phi, show_colorbar = TRUE)
```




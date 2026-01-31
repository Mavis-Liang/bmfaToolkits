# bmfaToolkits

Toolkits and reproducible workflows for Bayesian integrative factor analysis methods, with a lightweight “works out of the box” interface and optional backends for additional models.

Detail usage see https://mavis-liang.github.io/bmfaToolkits/. Navigate to the "Articles" tab for vignettes.

Can also refer to https://mavis-liang.github.io/Bayesian_integrative_FA_tutorial_book/. (Archived 1/31/2026)


## Install

### Development version (GitHub)

```r
install.packages("remotes")
remotes::install_github("Mavis-Liang/bmfaToolkits")
```

## Optional backends

Some methods depend on optional packages and/or vendored scripts. Use available_backends() to see what is usable on your system, and install extras as needed.

```
library(bmfaToolkits)
available_backends()
```

If a backend is missing, the corresponding fit_*() function will give a clear message telling you what to install.


## Quick start

```
library(bmfaToolkits)

# Example: MOM-SS (requires BFR.BE)
fit <- fit_momss(Y_list, k = 6, scaling = FALSE) # Y_list is a list of N_s x P data matrices
pp  <- postprocess_momss(fit)

# Plot an estimated shared covariance (example)
plot_single_loadings(pp$Phi, show_colorbar = TRUE)
```




#' Load the simulated nutrition dataset from this package
#'
#' The full simulated nutrition dataset used in the tutorial is shipped as an RDS
#' under `inst/extdata/` (to avoid loading it automatically).
#'
#' @return An object read by `readRDS()`.
#' @export
load_simulated_nutrition_data <- function() {
  path <- system.file("extdata", "simulated_nutrition_data.rds", package = "bifaToolkits")
  if (path == "") stop("extdata file not found (package not installed correctly).", call. = FALSE)
  readRDS(path)
}

#' Load the toy dataset from this package
#'
#' The toy datasetis shipped as an RDS
#' under `inst/extdata/`.
#'
#' @return An object read by `readRDS()`.
#' @export
load_toy_data <- function() {
  path <- system.file("extdata", "toy.rds", package = "bifaToolkits")
  if (path == "") stop("extdata file not found (package not installed correctly).", call. = FALSE)
  readRDS(path)
}

#' Load curated ovarian data (Bioconductor)
#'
#' This does not ship the data inside `bifaToolkits`. Instead, it provides a thin wrapper
#' that loads from the Bioconductor package `curatedOvarianData`.
#'
#' @param ... Passed to the underlying data accessor.
#' @return An object from `curatedOvarianData`.
#' @export
load_curated_ovarian <- function(...) {
  .require_pkg("curatedOvarianData", "loading curated ovarian datasets")
  # Users can call curatedOvarianData::loadOvarianEsets() etc.
  curatedOvarianData::loadOvarianEsets(...)
}

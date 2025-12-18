#' Example data shipped in extdata
#'
#' This package ships example data files under \code{inst/extdata/}.
#' Use the loader functions below to access them.
#'
#' Because these are stored in \code{extdata/}, they are not available via \code{data()}.
#'
#' @name bifa_extdata
NULL

#' Load the toy example object
#'
#' Loads the file \code{inst/extdata/toy.rds}.
#'
#' @return The R object stored in \code{toy.rds}. Often a list containing toy datasets.
#' @export
#' @examples
#' toy <- load_toy()
#' str(toy)
load_toy <- function() {
  f <- system.file("extdata", "toy.rds", package = "bifaToolkits", mustWork = TRUE)
  readRDS(f)
}

#' Load the simulated nutrition dataset
#'
#' Loads the file \code{inst/extdata/simulated_nutrition_data.rds}.
#'
#' @return The R object stored in \code{simulated_nutrition_data.rds}.
#' @export
#' @examples
#' nut <- load_simulated_nutrition_data()
#' str(nut)
load_simulated_nutrition_data <- function() {
  f <- system.file("extdata", "simulated_nutrition_data.rds", package = "bifaToolkits", mustWork = TRUE)
  readRDS(f)
}

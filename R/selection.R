#' Select number of factors via eigenvalue proportion rule
#'
#' This follows the simple rule used in the nutrition case study: choose K as the
#' number of eigenvalues explaining > `cutoff` proportion of total variance.
#'
#' @param Sigma A symmetric covariance matrix.
#' @param cutoff Proportion threshold (default 0.05).
#' @return Integer K.
#' @export
select_k_from_sigma <- function(Sigma, cutoff = 0.05) {
  val <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  prop <- val / sum(val)
  sum(prop > cutoff)
}

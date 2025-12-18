#' Fit stacked factor analysis (Stack FA) via MSFA::sp_fa
#'
#' This is the "Stack FA" baseline used in the tutorial: stack all studies and fit
#' a single factor model.
#'
#' @param Y_list List of study matrices (N_s x P).
#' @param k Number of factors.
#' @param scaling,centering Passed to `MSFA::sp_fa()`.
#' @param control Passed to `MSFA::sp_fa()`.
#' @param ... Passed to `MSFA::sp_fa()`.
#' @return A `bifa_fit` object.
#' @export
fit_stack_fa <- function(Y_list, k, scaling = FALSE, centering = TRUE,
                         control = list(nrun = 3000, burn = 2000), ...) {
  .require_pkg("MSFA", "Stack FA backend")
  Y_list <- .as_matrix_list(Y_list)
  .check_same_p(Y_list)
  Y_mat <- do.call(rbind, Y_list)
  fit <- MSFA::sp_fa(Y_mat, k = k, scaling = scaling, centering = centering,
                     control = control, ...)
  .new_bifa_fit("stack_fa", fit, meta = list(S = length(Y_list), P = ncol(Y_mat), k = k))
}

#' Post-process Stack FA output (nutrition tutorial)
#'
#' Implements the post-processing in the nutrition case study:
#' apply orthogonal Procrustes (OP) to posterior loadings, form SigmaPhi = Phi Phi^T,
#' estimate marginal covariance by averaging Sigma draws, and estimate Psi by averaging
#' Sigma - Lambda Lambda^T over MCMC iterations.
#'
#' @param fit A `bifa_fit` object from `fit_stack_fa()` or a raw `MSFA::sp_fa()` fit.
#' @param S Number of studies (used only to replicate marginal covariance into a list).
#' @return A list with `Phi`, `SigmaPhi`, `Psi`, `SigmaMarginal`.
#' @export
postprocess_stack_fa <- function(fit, S) {
  .require_pkg("MSFA", "Stack FA post-processing")
  fit <- .unwrap_fit(fit)
  est_Phi <- MSFA::sp_OP(fit$Lambda, trace = FALSE)$Phi
  est_SigmaPhi <- tcrossprod(est_Phi)
  est_SigmaMarginal <- lapply(seq_len(S), function(s) apply(fit$Sigma, c(1, 2), mean))
  Psi_chain <- vector("list", dim(fit$Sigma)[3])
  for (i in seq_len(dim(fit$Sigma)[3])) {
    Psi_chain[[i]] <- fit$Sigma[, , i] - tcrossprod(fit$Lambda[, , i])
  }
  est_Psi <- Reduce("+", Psi_chain) / length(Psi_chain)
  list(Phi = est_Phi, SigmaPhi = est_SigmaPhi, Psi = est_Psi, SigmaMarginal = est_SigmaMarginal)
}

#' Select K for Stack FA using eigen proportion rule
#' @param post Output of `postprocess_stack_fa()`.
#' @param cutoff Proportion threshold (default 0.05).
#' @export
select_k_stack_fa <- function(post, cutoff = 0.05) {
  select_k_from_sigma(post$SigmaPhi, cutoff = cutoff)
}

#' Fit Stack FA all-at-once
#' If `post_fit0` is NULL, fit an initial StackFA with given `k`, then automatically select
#' the number of factors and refit.
#'
#' @param Y_list Data list.
#' @param post_fit0 Initial post-processed result (or NULL to compute internally).
#' @param k Initial number of factors (used only if `post_fit0` is NULL).
#' @param cutoff Eigen proportion cutoff.
#' @param ... Passed to `fit_stack_fa()`.
#' @export
fit_stack_fa_2step <- function(Y_list, post_fit0 = NULL, k=NULL, cutoff = 0.05, ...) {
  if (is.null(post_fit0)) {
    fit0 <- fit_stack_fa(Y_list, k = k, ...)
    post_fit0 <- postprocess_stack_fa(fit0, S = length(.as_matrix_list(Y_list)))
  }
  K <- select_k_stack_fa(post_fit0, cutoff = cutoff)
  message(sprintf("Selected K = %d\n", K))
  fit2 <- fit_stack_fa(Y_list, k = K, ...)
  postprocess_stack_fa(fit2, S = length(.as_matrix_list(Y_list)))
}
